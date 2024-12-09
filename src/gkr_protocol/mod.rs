use ark_poly::SparseMultilinearExtension;
use rand::distributions::Slice;
use rand_chacha::rand_core::le;
use std::marker::PhantomData;
use ark_ff::PrimeField;

use crate::{side::to_le_indices};

pub mod ends;
mod tests;
#[derive(Debug)]
pub struct Wiring<const s: usize> {
    current_idx: usize,
    left: usize,
    right: usize,
}

impl <const s: usize> Wiring<s> {
    fn new (current_idx: usize, left: usize, right: usize) -> Self {
        assert!(current_idx < 1 << s);
        assert!(left < 1 << s);
        assert!(right < 1 << s);
        Self {
            current_idx,
            left,
            right,
        }
    }
}

#[derive(Debug)]
pub struct Layering<const s: usize> {
    pub add: Vec<Wiring<s>>,
    pub multiply: Vec<Wiring<s>>,
}

impl <const s: usize> Layering<s> {
    pub fn new(add: Vec<Wiring<s>>, multiply: Vec<Wiring<s>>) -> Self {
        Self {
            add,
            multiply,
        }
    }
}

impl <F: PrimeField, const s: usize> From<&Layering<s>> for [SparseMultilinearExtension<F>; 2] {
    /// Assume a uniform circuit first
    fn from(val: &Layering<s>) -> Self {
        let leindices = to_le_indices(3 * s);
        let addindices: Vec<usize> = val
            .add
            .iter()
            .map(|w| {
                // cook the index from the current index, left index and right index
                let index_be: usize = (w.current_idx << (2*s)) + (w.left << (s)) + w.right;
                leindices[index_be]
            })
            .collect();
        let add_evals: Vec<(usize, F)> = addindices.iter().map(|i| (*i, F::one())).collect();
        let add_mle = SparseMultilinearExtension::from_evaluations(3 * s, &add_evals);

        let mul_indices: Vec<usize> = val
            .multiply
            .iter()
            .map(|w| {
                // cook the index from the current index, left index and right index
                let indexbe: usize = (w.current_idx << (2 * s)) + (w.left << (s)) + w.right;
                leindices[indexbe]
            }).collect();
        let mul_evals: Vec<(usize, F)> = mul_indices.iter().map(|i| (*i, F::one())).collect();
        let mul_mle = SparseMultilinearExtension::from_evaluations(3 * s, &mul_evals);
        // returning both add and mul evaluations
        [add_mle, mul_mle]
    }
}

pub struct UniformCircuit<F, const s: usize> {
    layers: Vec<Layering<s>>,
    phantom: PhantomData<F>,
}

impl <F: PrimeField, const s: usize> UniformCircuit<F, s> {
    /// Padding each layer to have 1 << s elements,
    /// in order to have a uniform circuit representation.
    /// We do this by adding padded addition gates
    pub fn new (layers: Vec<Layering<s>>) -> Self {
        let mut layers = layers;
        for layer in layers.iter_mut() {
            let num_of_gates = layer.add.len() + layer.multiply.len();
            let diff = (1 << s) - (num_of_gates);
            for i in 0..diff {
                layer.add.push(Wiring::new(i + num_of_gates, 0, 0));
            }
        }

        Self {
            layers,
            phantom: PhantomData::<F>,
        }
    }

    /// this is a plain circuit evaluation given the input x 
    /// asserts that each layer is uniform with 1 << s evaluations
    pub fn evaluate(&self, x: Vec<F>) -> Vec<Vec<F>> {
        let mut evals = Vec::new();
        let mut last_layer = x;

        for layer in self.layers.iter().rev() {
            let mut new_layer: Vec<F> = vec![F::zero(); layer.add.len() + layer.multiply.len()];
            
            // handle addition
            for Wiring {
                current_idx,
                left,
                right,
            } in layer.add.iter(){
                new_layer[*current_idx] = last_layer[*left] + last_layer[*right];
            }

            // handle multiplication
            for Wiring {
                current_idx,
                left,
                right
            } in layer.multiply.iter()
            {
                new_layer[*current_idx] = last_layer[*left] * last_layer[*right];
            }

            assert_eq!(new_layer.len(), 1 << s, "non uniform circuit");
            evals.push(new_layer.clone());

            last_layer = new_layer;
        }
        evals
    }
}