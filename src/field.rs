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

