use std::iter::Once;

use ark_crypto_primitives::sponge::{
    poseidon::{PoseidonConfig, PoseidonSponge},
    Absorb, CryptographicSponge,
};
use ark_ff::PrimeField;
use ark_std::{iterable::Iterable, test_rng};

// Define publicly usable constants
pub const FULL_ROUNDS: usize = 8;
pub const PARTIAL_ROUNDS: usize = 31;
pub const ALPHA: u64 = 17; // Assuming `alpha` is an exponent and should be a positive integer

pub fn get_mds<F: PrimeField>() -> Vec<Vec<F>> {
    vec![
        vec![F::one(), F::zero(), F::one()],
        vec![F::one(), F::one(), F::zero()],
        vec![F::zero(), F::one(), F::one()],
    ]
}

pub fn get_round_constants<F: PrimeField>() -> Vec<Vec<F>> {
    let mut v = Vec::new();
    let mut ark_rng = test_rng();
    for _ in 0..(FULL_ROUNDS + PARTIAL_ROUNDS) {
        let mut res = Vec::new();
        for _ in 0..3 {
            res.push(F::rand(&mut ark_rng));
        }
        v.push(res);
    }
    v
}

pub(crate) fn test_sponge<F: PrimeField>() -> PoseidonSponge<F> {
    let mds = get_mds::<F>();
    let round_constants = get_round_constants::<F>();

    let config = PoseidonConfig::new(FULL_ROUNDS, PARTIAL_ROUNDS, ALPHA, mds, round_constants, 2, 1);
    PoseidonSponge::new(&config)
}

/// `SumCheckProof` can be securely shared with the verifier in full
#[derive(Clone, Debug, Default)]
pub struct SumCheckProof<F: PrimeField> {
    pub values: Vec<Vec<F>>,
    pub claimed: Option<F>, 
}

// Each side will initiate a Transcript
// The prover will write into the transcript the values of the polynomial
// The verifier will write into the transcript the values received via the `SumCheckProof`
#[derive(Clone, Debug, Default)]
pub struct Transcript<F: PrimeField + Absorb> {
    pub proof: SumCheckProof<F>,
}

impl<F: PrimeField + Absorb> Transcript<F> {
    pub fn new() -> Self {
        Self {
            proof: SumCheckProof::default(),
        }
    }

    pub fn update(&mut self, elements: &[F], sponge: &mut PoseidonSponge<F>) -> F {
        for elem in elements.iter() {
            sponge.absorb(elem);

        }

        self.proof.values.push(elements.to_vec());

        sponge.squeeze_field_elements::<F>(1)[0]
    }

    pub fn freeze_claim(&mut self, claimed_value: F) {
        self.proof.claimed = Some(claimed_value);
    }
}

// Separate a bit string represented as a usize into three values
// the last of which correspond to
#[inline]
pub fn to_zxy_form(idx: usize, v: usize) -> (usize, usize, usize) {
    let vp = 1 << v;

    let y = idx % vp;
    let idx_zx = idx / vp;
    let x = idx_zx % vp;
    let z = idx_zx / vp;

    (z, x, y)
}

/// perform a polynomial interpolation of the univariate polynomial of degree
/// `atmost` pi.len()-1 passing through the y-axis in pi at x = 0,1,2,...pi.len()-1
/// and evaluate the polynomial at `eval_at`.
pub(crate) fn univariate_polynomial_interpolation<F: PrimeField>(p_i: &[F], evaluate_at: F) -> F {
    let pilen = p_i.len();
    
    let mut evals = vec![];
    let mut prod = evaluate_at;
    evals.push(evaluate_at);

    // ∏ = ∏_j (evaluate_at - j)
    let mut check = F::zero();
    for i in 1..pilen {
        if evaluate_at == check {
            return p_i[i -1];
        }
        check += F::one();

        let temp = evaluate_at - check;
        evals.push(temp);
        prod *= temp;
    }

    if evaluate_at == check {
        return p_i[pilen - 1];
    }

    let mut res = F::zero();

    res

}