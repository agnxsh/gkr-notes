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
    pub(crate) fn new(sum_of_products: GKRProverSumOfProduct<F, MLE>, g: Vec<F>) -> Self {
        Self {
            sum_of_products, g
        }
    }
}

impl <F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for PolynomialOracle<F, MLE> {
    fn ithroundops(&self, x: &[F], y: &[F]) -> F {
        let mut zxy: Vec<F> = self.g.clone();
        zxy.extend(x);
        zxy.extend(y);
        let mut sum = F::zero();

        for MLEProduct(f1, f2, f3) in &self.sum_of_products.elements {
            sum += f1.evaluate(&zxy).unwrap() * f2.evaluate(x).unwrap() * f3.evaluate(y).unwrap()
        }
        sum
    }

    // return the number of variables of each of the 2 polynomials
    fn num_rounda(&self) -> usize {
        let MLEProduct(_,f2,_) = &self.sum_of_products.elements[0];
        2*f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}