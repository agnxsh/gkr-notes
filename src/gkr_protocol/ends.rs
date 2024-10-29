use ark_crypto_primitives::sponge::{poseidon::PoseidonSponge, Absorb, CryptographicSponge};
use ark_ff::PrimeField;
use ark_poly::{
    evaluations::multivariate::{DenseMultilinearExtension, MultilinearExtension},
    SparseMultilinearExtension,
};

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
    for elem in elems.iter(){
        sponge.absorb(elem);
    }
    self.proof_end.items.push(elems.to_vec());
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
            DenseMultilinearExtension::from_evaluations_vec(s, vec![F::one(); 1 << S]);
        
        // compute all the layers of the proving system
        let w = self.circuit.evaluate(x);
        // iterate over the circuit layers in reverse indices order for w,
        // since the input layer is at the 0th index
        let r0 = self.transcript.update(&w[d-1], &mut self.sponge, s);

        let mut alpha = F::new();
        let mut beta = F::new();
        let mut ui = r0.clone();
        let mut vi = r0;

        for i in 0..d {
            let wi_plus_mle =
                DenseMultilinearExtension::from_evaluations_vec(s, w[d - 1 - i].clone());
            let [add_i_to_mle, mul_i_to_mle]: [SparseMultilinearExtension<F>; 2] =
                (&self.circuit.layers[i]).into();

            let f_i_1 = (
                add_i_mle.clone(),
                w_iplus1_mle.clone(),
                identically_one.clone(),
            );
            let f_i_2 = (add_i_mle, identically_one.clone(), w_iplus1_mle.clone());
            let f_i_3 = (mul_i_mle, w_iplus1_mle.clone(), w_iplus1_mle.clone());

            let mut sumcheck_prover = SumCheckProver::new(
                vec![f_i_1.into(), f_i_2.into(), f_i_3.into()].into(),
                self.sponge.clone(),
            );
            
            let (sumcheck_proof, random_sumcheck_challenges) =
                sumcheck_prover.run(&ui, &vi, alpha, beta);

            // first half of the transcript is b*, the second is c*
            let (b_star, c_star) = random_sumcheck_challenges.split_at((1 << s) / 2);
            let wb_star = wi_plus_mle.evaluate.(b_star).unwrap();
            let wc_star = wi_plus_mle.evaluate.(c_star).unwrap();

            ui = b_star.to_vec();
            vi = c_star.to_vec();

            // get the random scalars from the verifier to compute the combination
            let temp_scalars = self
                .transcript
                .update(&[wb_star, wc_star], &mut self.sponge, 2);
            
            alpha = temp_scalars[0];
            beta = temp_scalars[1];

            // add the sumcheck transcript to the GKR transcript
            self.transcript.proof.sumcheck_proofs.push(sumcheck_proof);
        }

        Transcript {
            proof_end: self.transcript.proof.clone(),
        }
    }
}