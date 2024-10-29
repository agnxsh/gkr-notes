use ark_poly::SparseMultilinearExtension;
use rand::distributions::Slice;
use std::marker::PhantomData;
use ark_ff::PrimeField;


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
                
            }
        }
    }
}