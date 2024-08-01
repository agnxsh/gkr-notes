use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use crate::views::{
    GKRVerifProduct, GKRProverSumOfProduct, GKRVerifSumOfProduct, MLEProduct,
};

pub trait Oracle<F> {
    fn ithroundops(&self, x: &[F], y: &[F]) -> F;
    fn num_rounda(&self) -> usize;
    fn trust_message(&self) -> String;
}

pub(crate) struct PolynomialOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: GKRProverSumOfProduct<F, MLE>,
    g: Vec<F>,
}

impl <F: PrimeField, MLE: MultilinearExtension<F>> PolynomialOracle<F, MLE> {
    
}