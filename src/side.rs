use std::{collections::VecDeque, iter::Product};

use ark_crypto_primitives::sponge::poseidon::PoseidonSponge;
use ark_crypto_primitives::sponge::Absorb;
use ark_ff::PrimeField;

use ark_poly::{MultilinearExtension, SparseMultilinearExtension};

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

pub(crate) fn phase2_start <F: PrimeField, MLE: MultilinearExtension<F>> (
    f_1: &SparseMultilinearExtension<F>,
    g1: &[F],
    g2: &[F],
    u: &[F],
    alpha: F,
    beta: F,
) -> Vec<F> {
    let v
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
    pub fn run(&mut self, g1: &[F], g2: &[F], alpha: F, beta: F) -> (SumCheckProof<F, Vec<F>) {
        let MLEProduct(_, pol, _) = &self.sum_of_products.elements[0];
        let k: usize = pol.num_vars();

        let mut A_hs: Vec<Vec<F>> = Vec::new();
        let mut A_f2s: Vec<Vec<F>> = Vec::new();

        //phase 1
        for MLEProduct(f1, f2, f3) in self.sum_of_products.iter() {
            A_hs.push()
        }
    }
}