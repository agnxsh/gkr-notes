use ark_ff::PrimeField;
use ark_poly::{
    MultilinearExtension,
    SparseMultilinearExtension
};

#[derive(Clone, Debug)]
pub struct GKRVerifProduct<F: PrimeField> (
    pub SparseMultilinearExtension<F>,
    pub F,
    pub F
);

impl<F: PrimeField> From <(SparseMultilinearExtension<F>, F, F)> for GKRVerifProduct<F> {
    fn from((sme, f1, f2): (SparseMultilinearExtension<F>, F, F)) -> Self {
        GKRVerifProduct(sme, f1, f2)
    }
}

// Computing the product of Spare Multilinear Extension
// with 2 other Multilinear Extensions. When dealing with GKR
#[derive(Clone, Debug)]
pub struct MLEProduct<F: PrimeField, MLE: MultilinearExtension<F>> (
    pub SparseMultilinearExtension<F>,
    pub MLE,
    pub MLE,
);

impl <F: PrimeField, MLE: MultilinearExtension<F>> From<(SparseMultilinearExtension<F>, MLE, MLE)>
    for MLEProduct<F, MLE>
    {
        fn from(tuple: (SparseMultilinearExtension<F>, MLE, MLE)) -> Self {
            assert_eq!(
                tuple.1.num_vars(),
                tuple.2.num_vars(),
                "f3 and f2 should be having number of variables"
            );
            Self(tuple.0, tuple.1, tuple.2)
        }
    }


/// Sum of products as seen by the verifier
/// For eg, for GKR, we write the verifier's view
/// in a multivariate polynomial's General-Purpose
/// Form
#[derive(Clone, Debug)]
pub struct GKRVerifSumOfProduct<F: PrimeField> {
    pub elements: Vec<GKRVerifProduct<F>>
}

impl <F: PrimeField> From<Vec<GKRVerifProduct<F>>> for GKRVerifSumOfProduct<F> {
    fn from (elements: Vec<GKRVerifProduct<F>>) -> Self {
        Self {
            elements
        }
    }
}

/// Sum of products as seen by the prover's view
#[derive(Clone, Debug)]
pub struct GKRProverSumOfProduct<F: PrimeField, MLE: MultilinearExtension<F>> {
    pub elements: Vec<MLEProduct<F, MLE>>,
}

impl <F: PrimeField, MLE: MultilinearExtension<F>> From <Vec<MLEProduct<F, MLE>>>
    for GKRProverSumOfProduct<F, MLE>

{
    fn from (elements: Vec<MLEProduct<F, MLE>>) -> Self {
        Self {
            elements
        }
    }
}

