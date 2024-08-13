use std::iter::Product;

use ark_ff::PrimeField;
use ark_poly::MultilinearExtension;
use crate::views::{
    GKRVerifProduct, GKRProverSumOfProduct, GKRVerifSumOfProduct, MLEProduct,
};

pub trait Oracle<F> {
    fn ithroundops(&self, x: &[F], y: &[F]) -> F;
    fn num_rounds(&self) -> usize;
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
    fn num_rounds(&self) -> usize {
        let MLEProduct(_,f2,_) = &self.sum_of_products.elements[0];
        2*f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

pub(crate) struct CombinedPolynomialOracle<F: PrimeField, MLE: MultilinearExtension<F>> {
    sum_of_products: GKRProverSumOfProduct<F, MLE>,
    g1: Vec<F>,
    g2: Vec<F>,
    alpha: F,
    beta: F,
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> CombinedPolynomialOracle<F, MLE> {
    pub(crate) fn new(
        sum_of_products: GKRProverSumOfProduct<F, MLE>,
        g1: Vec<F>,
        g2: Vec<F>,
        alpha: F,
        beta: F,
    ) -> Self {
        Self {
            sum_of_products,
            g1,
            g2,
            alpha,
            beta,
        }
    }
}

impl<F: PrimeField, MLE: MultilinearExtension<F>> Oracle<F> for CombinedPolynomialOracle<F, MLE> {
    fn ithroundops(&self, x: &[F], y: &[F]) -> F {
        let mut sums = Vec::new();

        for (coeff, g) in vec![(self.alpha, self.g1.clone()), (self.beta, self.g2.clone())] {
            let mut zxy: Vec<F> = g.clone();
            zxy.extend(x);
            zxy.extend(y);

            let mut sum = F::zero();

            for MLEProduct(f1, f2, f3) in &self.sum_of_products.elements {
                sum += f1.evaluate(&zxy).unwrap() * f2.evaluate(x).unwrap() * f3.evaluate(y).unwrap()
            }

            sums.push(coeff * sum);
        }
        
        sums[0] + sums[1]
    }

    // returns the number of variables of each of the 2 polynomials
    // the total number of variables of the product is twice that
    fn num_rounds(&self) -> usize {
        let MLEProduct(_, f2, _) = &self.sum_of_products.elements[0];
        2 * f2.num_vars()
    }

    fn trust_message(&self) -> String {
        String::from("unconditionally")
    }
}

