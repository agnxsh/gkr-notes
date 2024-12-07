use std::vec;

use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};
use rand_chacha::rand_core::le;

use super::UniformCircuit;
use crate::helpers::SumCheckProof;
use crate::{
    oracle::GKROracle,
    side::{Prover as SumCheckProver, Verifier as SumCheckVerifier},
};

#[derive(Clone, Debug, Default)]
pub struct GKRProofEnd<F: PrimeField> {
    pub items: Vec<Vec<F>>,
    pub sumcheck_proofs: Vec<SumCheckProof<F>>,
}

// Each end, be it the prover or the verifier's end will instantiate a `Transcript`
// The prover end will write the GKR layer claims into the transcript in order to 
// make the overall the proving system non-interactive in nature, the prover will do 
// so entirely separately as sumcheck messages, thereby constructing the `GKRProofEnd`
// within the transcript.

// The verifier end in turn will feed into the transcript the values received in the
// `GKRProof` itself.

#[derive(Clone, Debug)]
pub struct Transcript<F: PrimeField + Absorb> {
    pub proof_end: GKRProofEnd<F>,
}

impl<F: PrimeField + Absorb> Transcript <F> {
   fn new() -> Self {
    Self {
        proof_end: GKRProofEnd::default(),
    }
   }

   fn update(&mut self, elements: &[F], sponge: &mut PoseidonSponge<F>, n: usize) -> Vec<F> {
    for elem in elements.iter(){
        sponge.absorb(elem);
    }
    self.proof_end.items.push(elements.to_vec());
    sponge.squeeze_field_elements::<F>(n)
   }
}

pub struct Prover<F: PrimeField + Absorb, const s: usize> {
    circuit: UniformCircuit<F, s>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

impl<F: PrimeField + Absorb, const s: usize> Prover<F, s> {
    pub fn new(circuit: UniformCircuit<F, s>, sponge: PoseidonSponge<F>) -> Self {
        Self {
            circuit,
            sponge,
            transcript: Transcript::new(),
        }
    }
}

impl<F: PrimeField + Absorb, const s: usize> Prover <F, s> {
    pub fn run(&mut self, v: Vec<F>) -> Transcript<F> {
        // d defines the number of layers of the proving system
        let d = self.circuit.layers.len();

        let distinct_one =
            DenseMultilinearExtension::from_evaluations_vec(s, vec![F::one(); 1 << s]);
        
        // compute all the layers of the proving system
        let w = self.circuit.evaluate(v);
        // iterate over the circuit layers in reverse indices order for w,
        // since the input layer is at the 0th index
        let r0 = self.transcript.update(&w[d-1], &mut self.sponge, s);

        let mut alpha = F::one();
        let mut beta = F::zero();
        let mut ui = r0.clone();
        let mut vi = r0;

        for i in 0..d {
            let wi_plus_mle =
                DenseMultilinearExtension::from_evaluations_vec(s, w[d - 1 - i].clone());
            let [add_i_to_mle, mul_i_to_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i]).into();

            let f_i_1 = (
                add_i_to_mle.clone(),
                wi_plus_mle.clone(),
                distinct_one.clone(),
            );
            let f_i_2 = (add_i_to_mle, distinct_one.clone(), wi_plus_mle.clone());
            let f_i_3 = (mul_i_to_mle, wi_plus_mle.clone(), wi_plus_mle.clone());

            let mut sumcheck_prover = SumCheckProver::new(
                vec![f_i_1.into(), f_i_2.into(), f_i_3.into()].into(),
                self.sponge.clone(),
            );
            
            let (sumcheck_proof, random_sumcheck_challenges) =
                sumcheck_prover.run(&ui, &vi, alpha, beta);

            // first half of the transcript is b*, the second is c*
            let (b_star, c_star) = random_sumcheck_challenges.split_at((1 << s) / 2);
            let wb_star = wi_plus_mle.evaluate(b_star).unwrap();
            let wc_star = wi_plus_mle.evaluate(c_star).unwrap();

            ui = b_star.to_vec();
            vi = c_star.to_vec();

            // get the random scalars from the verifier to compute the combination
            let temp_scalars = self
                .transcript
                .update(&[wb_star, wc_star], &mut self.sponge, 2);
            
            alpha = temp_scalars[0];
            beta = temp_scalars[1];

            // add the sumcheck transcript to the GKR transcript
            self.transcript.proof_end.sumcheck_proofs.push(sumcheck_proof);
        }

        Transcript {
            proof_end: self.transcript.proof_end.clone(),
        }
    }
}

pub struct Verifer<F: PrimeField + Absorb, const s: usize> {
    circuit: UniformCircuit<F, s>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

impl <F: PrimeField + Absorb, const s: usize> Verifer<F, s> {
    pub fn new(
        circuit: UniformCircuit<F, s>,
        sponge: PoseidonSponge<F>,
        proof_end: GKRProofEnd<F>,
    ) -> Self {
        Self {
            circuit,
            sponge,
            transcript: Transcript { proof_end },
        }
    }

    pub fn run(&mut self) -> bool {
        let d = self.circuit.layers.len();
        let claim = self.transcript.proof_end.items[0].clone();
        let r0 = self.transcript.update(&claim, &mut self.sponge, s);

        // verifying the transcript correctly
        let mut ui = r0.clone();
        let mut vi = r0;

        let mut alpha = F::one();
        let mut beta = F::zero();

        for i in 0..d {
            let [add_i_mle, mul_i_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i]).into();

            let wui = self.transcript.proof_end.items[i+1][0];
            let wvi = self.transcript.proof_end.items[i+1][1];

            let fi1 = (add_i_mle.clone(), wui, F::one());
            let fi2 = (add_i_mle, F::one(), wvi);
            let fi3 = (mul_i_mle, wui, wvi);
            let sum_of_prods = vec![fi1.into(), fi2.into(), fi3.into()].into();

            let mut sumcheck_verifier = SumCheckVerifier::new(
                GKROracle::new(sum_of_prods, ui, vi, alpha, beta),
                self.sponge.clone(),
                self.transcript.proof_end.sumcheck_proofs[i].clone(),
            );

            let (result, random_challenges) = sumcheck_verifier.run();
            assert!(result, "sumcheck has failed at round {}", i);

            let (tv1, tv2) = random_challenges.split_at((1 << s) / 2);
            ui = tv1.to_vec();
            vi = tv2.to_vec();

            let temp_scalars = self.transcript.update(&[wui, wvi], &mut self.sponge, 2);

            alpha = temp_scalars[0];
            beta = temp_scalars[1];
        }

        true
    }
}