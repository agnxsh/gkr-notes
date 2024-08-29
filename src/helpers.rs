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
}