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

pub (crate) fn to_le_indices(v: usize) -> Vec<usize> {
    (0usize..(1 << v))
        .map(|i| i.reverse_bits() >> (usize::BITS as usize - v))
        .collect()
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


            }
        }


        (claimed_value, random_challenges)
    }
}