pub mod views;
#[allow(non_upper_case_globals)]
pub mod oracle;
pub mod gkr_protocol;
pub mod side;
mod helpers;
mod tests;

use ark_ff::fields::{Fp64, MontBackend, MontConfig};

/// Configuration for the finite field Fq.
/// - Modulus: A prime number defining the field's order.
/// - Generator: A multiplicative generator for the field.
///
/// The modulus is chosen to be 101 (a small prime), and the generator is 2.
#[derive(MontConfig)]
#[modulus = "101"]
#[generator = "2"]
pub struct FqConfig;

/// Type alias for the finite field Fq.
/// This uses `Fp64`, which is optimized for fields with modulus less than 2^64.
pub type Fq = Fp64<MontBackend<FqConfig, 1>>;
