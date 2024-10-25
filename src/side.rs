use std::{collections::VecDeque, iter::Product, result};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};
use rand::random;

use crate:: {
    views::{GKRVerifSumOfProduct, GKRVerifProduct, MLEProduct, GKRProverSumOfProduct},
    oracle::{CombinedPolynomialOracle, GKROracle, Oracle, PolynomialOracle},
    helpers::{univariate_polynomial_interpolation, test_sponge, to_zxy_form, SumCheckProof, Transcript},
};

pub struct Prover<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> {
    sum_of_products: GKRProverSumOfProduct<F, MLE>,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

pub struct Verifier<F: PrimeField + Absorb, O: Oracle<F>> {
    oracle: O,
    sponge: PoseidonSponge<F>,
    transcript: Transcript<F>,
}

impl<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> Prover<F, MLE> {
    pub fn new(sum_of_products: GKRProverSumOfProduct<F, MLE>, sponge: PoseidonSponge<F>) -> Self {
        Self {
            sum_of_products,
            sponge,
            transcript: Transcript::new(),
        }
    }

    #[allow(non_snake_case)]
    fn sumcheck_product(
        &mut self,
        A_fs: &mut Vec<Vec<F>>,
        A_gs: &mut Vec<Vec<F>>,
        v: usize,
    ) -> (F, Vec<F>) {
        // in the first round the claimed_value contains the prover wants to prove
        let num_summands = A_fs.len();
        assert_eq!(num_summands, A_gs.len());
        let claimed_value = (0..num_summands)
            .map(|summand| {
                (0..1 << v)
                    .map(|i| A_fs[summand][i] * A_gs[summand][i])
                    .sum::<F>()
            }).sum();

        let le_indices = to_le_indices(v);
        let mut random_challenges = Vec::new();

        for i in 1..=v {
            let mut values = vec![F::zero(); 3];

            for summand in 0..num_summands {
                let mut f_values = vec![vec![F::zero(); 1 << (v-1)]; 3];
                let mut g_values = vec![vec![F::zero(); 1 << (v-1)]; 3];
                let A_f = &A_fs[summand];
                let A_g = &A_gs[summand];

                for b in 0..(1 << (v-1)) {
                    f_values[0][b] = A_f[le_indices[b]];
                    f_values[1][b] = A_f[le_indices[b + (1 << (v - i))]];
                    f_values[2][b] = -A_f[le_indices[b]] + A_f[le_indices[b + (1 << (v - i))]] * F::from(2u64);

                    g_values[0][b] = A_g[le_indices[b]];
                    g_values[1][b] = A_g[le_indices[b + (1 << (v - i))]];
                    g_values[2][b] =  -A_g[le_indices[b]] + A_g[le_indices[b + (1 << (v - i))]] * F::from(2u64);                
                }

                let temp_vals : Vec<F> = (0..=2)
                    .map(|t| ((0..1 << (v - i)).map(|b| f_values[t][b] * g_values[t][b])).sum())
                    .collect();

                values[0] += temp_vals[0];
                values[1] += temp_vals[1];
                values[2] += temp_vals[2];
            }

            let r_i = self.transcript.update(&values, &mut self.sponge);
            random_challenges.push(r_i);

            for summand in 0..num_summands {
                let A_f = &mut A_fs[summand];
                let A_g = &mut A_gs[summand];

                for b in 0..(1 << (v - i)) {
                    A_f[le_indices[b]] = A_f[le_indices[b]] * (F::one() - r_i)
                        + A_f[le_indices[b + (1 << (v - i))]] * r_i;
                    A_g[le_indices[b]] = A_g[le_indices[b]] * (F::one() - r_i)
                        + A_g[le_indices[b + (1 << (v - i))]] * r_i;
                }
            }
        }
        (claimed_value, random_challenges)
    }

    #[allow(non_snake_case)]
    pub fn run(&mut self, g1: &[F], g2: &[F], alpha: F, beta: F) -> (SumCheckProof<F>, Vec<F>) {
        let MLEProduct(_, pol, _) = &self.sum_of_products.elements[0];
        let k: usize = pol.num_vars();

        let mut A_hs: Vec<Vec<F>> = Vec::new();
        let mut A_f2s: Vec<Vec<F>> = Vec::new();

        //phase 1
        for MLEProduct(f1, f2, f3) in self.sum_of_products.elements.iter() {
            A_hs.push(phase1_start(f1, f3, g1, g2, alpha, beta));
            A_f2s.push(f2.to_evaluations());
        }

        let (out, u) = self.sumcheck_product(&mut A_hs, &mut A_f2s, k);

        self.transcript.freeze_claim(out);

        //phase 2
        let _A_f3s: Vec<Vec<F>> = Vec::new();
        let mut A_f1s: Vec<Vec<F>> = Vec::new();

        let mut A_f3_f2_us: Vec<Vec<F>> = Vec::new();
        for MLEProduct(f1, f2, f3) in &self.sum_of_products.elements {
            let f2_u = f2.evaluate(&u).unwrap();
            let A_f3 = f3.to_evaluations();

            A_f1s.push(phase2_start::<F>(f1, g1, g2, &u, alpha, beta));
            A_f3_f2_us.push(A_f3.iter().map(|x| *x * f2_u).collect::<Vec<F>>());
        }

        let (_,v) = self.sumcheck_product(&mut A_f1s, &mut A_f3_f2_us, k);

        println!("Prover finished successfully");

        (self.transcript.proof.clone(), [u, v].concat())
    }
}

impl<F: PrimeField + Absorb, O: Oracle<F>> Verifier<F, O> {
    pub fn new(oracle: O, sponge: PoseidonSponge<F>, proof: SumCheckProof<F>) -> Self {
        Self {
            oracle,
            sponge,
            transcript: Transcript { proof },
        }
    }

    pub fn run(&mut self) -> (bool, Vec<F>) {
        let v = self.oracle.num_rounds();

        // the last polynomial is initialised to a singleton list containing the claimed
        // value in subsequent rounds, it contains the values of the expected polynomial
        // at 0, 1 and 2.
        let mut last_polynomial = vec![self.transcript.proof.claimed.unwrap()];

        // dummy scalar to be used when evaluating the constant polynomial in the first
        // iteration
        let mut last_sent_scalar = F::one();

        let mut random_challenges = Vec::new();

        for r in 0..v {
            let new_polynomial_evals = self.transcript.proof.values[r].clone();

            let claimed_sum = new_polynomial_evals[0] + new_polynomial_evals[1];
            let new_polynomial = Vec::from(&new_polynomial_evals[..]);

            if claimed_sum != univariate_polynomial_interpolation(&last_polynomial, last_sent_scalar) {
                return (false, random_challenges);
            }
            last_sent_scalar = self.transcript.update(&new_polynomial_evals, &mut self.sponge);
            random_challenges.push(last_sent_scalar);

            last_polynomial = new_polynomial;
        }

        if univariate_polynomial_interpolation(&last_polynomial, last_sent_scalar) 
            != self
                .oracle
                .ithroundops(&random_challenges[..(v / 2)], &random_challenges[(v / 2)..]) 
            
        {
            println!("verifier found inconsistent computation in the oracle call");
            return (false, random_challenges);
        }
        (true, random_challenges)
    }
}

pub(crate) fn precompute_table<F: PrimeField>(vals: &[F]) -> Vec<F> {
    let n = vals.len();

    let mut table_1 = vec![F::zero(); 1 << n];
    let mut table_2 = vec![F::zero(); 1 << n];

    let mut prev_table = &mut table_1;
    let mut next_table = &mut table_2;

    prev_table[0] = F::one();

    for i in 0..n {
        for b in 0..(1 << i) {
            next_table[2 * b] = prev_table[b] * (F::one() - vals[i]);
            next_table[2 * b + 1] = prev_table[b] * vals[i];
        }

        std::mem::swap(&mut prev_table, &mut next_table);
    }

    if n % 2 == 0 {
        table_1
    } else {
        table_2
    }
}

pub (crate) fn to_le_indices(v: usize) -> Vec<usize> {
    (0usize..(1 << v))
        .map(|i| i.reverse_bits() >> (usize::BITS as usize - v))
        .collect()
}

pub(crate) fn phase1_start <F: PrimeField, MLE: MultilinearExtension<F>> (
    f_1: &SparseMultilinearExtension<F>,
    f_3: &MLE,
    g1: &[F],
    g2: &[F],
    alpha: F,
    beta: F,
) -> Vec<F> {
    let v = f_3.num_vars();
    let le_indices_for_f1 = to_le_indices(f_1.num_vars);
    let le_indices_for_f3 = to_le_indices(v);

    let table_g1 = precompute_table(g1);
    let table_g2 = precompute_table(g2);
    let table_g = table_g1
        .iter()
        .zip(table_g2.iter())
        .map(|(x, y)| alpha * *x + beta * *y)
        .collect::<Vec<F>>();

    let mut ahg = vec![F::zero(); 1 << v];

    for(idx_le, val) in f_1.evaluations.iter() {
        let idx = le_indices_for_f1[*idx_le];
        let (z, x, y) = to_zxy_form(idx, v);
        ahg[le_indices_for_f3[x]] += table_g[z] * val * f_3.to_evaluations()[le_indices_for_f3[y]];
    }
    ahg
}

pub(crate) fn phase2_start <F: PrimeField> (
    f_1: &SparseMultilinearExtension<F>,
    g1: &[F],
    g2: &[F],
    u: &[F],
    alpha: F,
    beta: F,
) -> Vec<F> {
    let v = g1.len();
    let le_indices_f1 = to_le_indices(f_1.num_vars);
    let le_indices_g = to_le_indices(v);

    let table_g1 = precompute_table(g1);
    let table_g2 = precompute_table(g2);
    let table_u = precompute_table(u);

    let mut af1 = vec![F::zero(); 1 << v];

    for (idx_le, val) in f_1.evaluations.iter() {
        let idx = le_indices_f1[*idx_le];
        let (z, x, y) = to_zxy_form(idx, v);
        af1[le_indices_g[y]] += (alpha * table_g1[z] + beta * table_g2[z]) * table_u[x] * val;
    }

    af1
}

// Run through the protocol and return true iff the verifier does NOT reject
pub fn run_sumcheck_protocol<F: PrimeField + Absorb, MLE: MultilinearExtension<F>> (
    f1: SparseMultilinearExtension<F>,
    f2: MLE,
    f3: MLE,
    g: &[F],
) {
    let simple_sum = GKRProverSumOfProduct {
        elements: vec![MLEProduct(f1, f2, f3)],
    };

    // Need to use the sample sponge, since it's initialized the random values
    let sponge = test_sponge();
    let mut prover = Prover::new(simple_sum.clone(), sponge.clone());

    let (proof, _) = prover.run(g, g, F::one(), F::zero());

    let mut verifier = Verifier::new(
        PolynomialOracle::<F, MLE>::new(simple_sum, g.to_vec()),
        sponge,
        proof,
    );

    let (result, _) = verifier.run();
    assert!(result)
}