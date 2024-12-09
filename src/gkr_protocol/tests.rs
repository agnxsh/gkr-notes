#[cfg(test)]
mod tests {
    use crate::{
        gkr_protocol::{
            ends::{Prover, Verifier},
            *,
        },
        helpers::test_sponge,
    };
    use ark_bls12_381::Fq;
    use ark_crypto_primitives::sponge::Absorb;

    fn cook_test_circuit<F: PrimeField + Absorb>() -> UniformCircuit<F, 2> {
        // multtply, layer 0 (o/p)
        let mul00 = Wiring::new(0, 0, 1);
        let mul01 = Wiring::new(1, 2, 3);

        let layer0 = Layering::new(Vec::new(), vec![mul00, mul01]);

        // add, layer 0
        // empty 
        let mul10 = Wiring::new(0, 0, 0);
        let mul11 = Wiring::new(1, 1, 1);
        let mul12 = Wiring::new(2, 1, 2);

        // add, layer 1
        let add13 = Wiring::new(3, 3, 3);
        let layer1 = Layering::new(vec![add13], vec![mul10, mul11, mul12]);

        UniformCircuit::<F, 2>::new(vec![layer0, layer1])
    }

    #[test]
    fn simple_circuit_cooker() {
        let circuit = cook_test_circuit::<Fq>();
        let computed = circuit.evaluate(
            vec![3,2,3,1]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect::<Vec<Fq>>()
        );
        assert_eq!(
            *computed.last().unwrap(),
            vec![36, 12, 18, 18]
                .iter()
                .map(|x| Fq::from(*x as u64))
                .collect::<Vec<Fq>>()
        );
    }

    #[test]
    fn layer_to_multilinerar_ext() {
        let circuit = cook_test_circuit::<Fq>();
        for layer in circuit.layers.iter() {
            let _: [SparseMultilinearExtension<Fq>; 2] = layer.into();
        }
    }

    #[test]
    fn test_gkr() {
        // Need to use the same sponge, since it's initialized with random values
        let sponge = test_sponge();

        let mut gkr_prover = Prover::<Fq, 2>::new(cook_test_circuit(), sponge.clone());
        let circuit_input = vec![3, 2, 3, 1]
            .iter()
            .map(|x| Fq::from(*x as u64))
            .collect();
        let transcript = gkr_prover.run(circuit_input);

        let mut gkr_verifier =
            Verifier::<Fq, 2>::new(cook_test_circuit(), sponge, transcript.proof_end);
        assert!(gkr_verifier.run());
    }
}
