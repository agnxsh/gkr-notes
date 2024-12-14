#[cfg(test)]
mod tests {
    use core::num;

    use crate::{
        helpers::to_zxy_form, side::{
            phase1_start as p1_org,
            phase2_start as p2_org, run_sumcheck_protocol,
            run_sumcheck_protocol_combined, run_sumcheck_protocol_combined_multiproof
        }, views::{GKRProverSumOfProduct, GKRVerifSumOfProduct, MLEProduct}, Fq
    };
    use ark_ff::PrimeField;
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension, SparseMultilinearExtension};
    use ark_std::UniformRand;

    fn start_testing_phase1<F: PrimeField, MLE: MultilinearExtension<F>>(
        f1: &SparseMultilinearExtension<F>,
        f3: &MLE,
        g: &[F],
    ) -> Vec<F> {
        p1_org(f1, f3, g, g, F::one(), F::zero())
    }

    fn start_testing_phase2<F: PrimeField>(
        f1: &SparseMultilinearExtension<F>,
        g: &[F],
        u: &[F],
    ) -> Vec<F> {
        p2_org(f1, g, g, u, F::one(), F::zero())
    }

    #[test]
    fn test_int_zxy() {
        let v = 2;
        let (z0, x0, y0) = (3, 0, 1);

        let vp = 1 << v;
        let idx = z0 * (vp * vp) + x0 * vp + y0;

        assert_eq!((z0, x0, y0), to_zxy_form(idx, v));
    }

    #[test]
    fn test_phase1_expected0() {
        let mut test_rng = ark_std::test_rng();
        
        // let f(x1, x2, x3, x4, x5, x6) = f(g1, g2, g3, g4, g5, g6)
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // assume |x| = |y| = |z| =  2, so f1 would be a 6-variate polynomial
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // when the first 2 variables are 1 and 1, then f1 should be 0 regardless of what f3 is
        let g = vec![Fq::from(1u64), Fq::from(1u64)];

        // f3 is a random 2 variate MLE
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            (0..4).map(|_| Fq::rand(&mut test_rng)).collect(),
        );

        let actual = start_testing_phase1(&f1, &f3, &g);

        // expected value
        let expected = vec![
            Fq::from(0u64),
            Fq::from(0u64),
            Fq::from(0u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);

    }

    #[test]
    fn test_phase1() {
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![3, 3, 3, 3]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let actual = start_testing_phase1(&f1, &f3, &g);

        // expected values
        let expected = vec![
            Fq::from(6u64),
            Fq::from(0u64),
            Fq::from(3u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_phase1_full() {
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let actual = start_testing_phase1(&f1, &f3, &g);

        // expected value
        let expected = vec![
            Fq::from(10u64),
            Fq::from(0u64),
            Fq::from(1u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_phase1_full_with_non_zero_gval() {
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);
        let g = vec![Fq::from(2u64), Fq::from(3u64)];
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let actual = start_testing_phase1(&f1, &f3, &g);

        // expected value
        let expected = vec![
            Fq::from(20u64),
            Fq::from(0u64),
            Fq::from(2u64),
            Fq::from(0u64),
        ];
        assert_eq!(actual, expected)
    }

    #[test]
    fn tests_phase1_max_value() {
        let points = vec![
            (4, Fq::from(6_u64)),
            (5, Fq::from(6_u64)),
            (6, Fq::from(6_u64)),
            (7, Fq::from(97_u64)),
            (12, Fq::from(10_u64)),
            (13, Fq::from(10_u64)),
            (14, Fq::from(10_u64)),
            (20, Fq::from(7_u64)),
            (21, Fq::from(7_u64)),
            (22, Fq::from(7_u64)),
            (23, Fq::from(98_u64)),
            (28, Fq::from(11_u64)),
            (29, Fq::from(11_u64)),
            (30, Fq::from(11_u64)),
            (31, Fq::from(1_u64)),
            (36, Fq::from(6_u64)),
            (37, Fq::from(6_u64)),
            (38, Fq::from(7_u64)),
            (39, Fq::from(98_u64)),
            (44, Fq::from(3_u64)),
            (45, Fq::from(3_u64)),
            (46, Fq::from(4_u64)),
            (47, Fq::from(95_u64)),
            (52, Fq::from(7_u64)),
            (53, Fq::from(7_u64)),
            (54, Fq::from(8_u64)),
            (55, Fq::from(99_u64)),
            (60, Fq::from(4_u64)),
            (61, Fq::from(4_u64)),
            (62, Fq::from(5_u64)),
            (63, Fq::from(96_u64)),
        ];

        // asssume |x| = |y| = 2, so f1 would be a 6-variate
        // |z| = 2
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        // when first 2 variables are 1, 1, then f1 should be 0 no matter what f3 is
        let g = vec![Fq::from(0u64), Fq::from(0u64)];

        // eval vector in LE indexing
        // f3(y1, y2) = 1 - 10 * y1 * y2 + 5 * y1 + 4 * y2
        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let actual = start_testing_phase1(&f1, &f3, &g);

        // expected values
        let expected = vec![
            Fq::from(0u64),
            Fq::from(78u64),
            Fq::from(0u64),
            Fq::from(91u64),
        ];
        assert_eq!(actual, expected);
    }

    #[test]
    fn test_phase2_max() {
        let points = vec![
            (4, Fq::from(6_u64)),
            (5, Fq::from(6_u64)),
            (6, Fq::from(6_u64)),
            (7, Fq::from(97_u64)),
            (12, Fq::from(10_u64)),
            (13, Fq::from(10_u64)),
            (14, Fq::from(10_u64)),
            (20, Fq::from(7_u64)),
            (21, Fq::from(7_u64)),
            (22, Fq::from(7_u64)),
            (23, Fq::from(98_u64)),
            (28, Fq::from(11_u64)),
            (29, Fq::from(11_u64)),
            (30, Fq::from(11_u64)),
            (31, Fq::from(1_u64)),
            (36, Fq::from(6_u64)),
            (37, Fq::from(6_u64)),
            (38, Fq::from(7_u64)),
            (39, Fq::from(98_u64)),
            (44, Fq::from(3_u64)),
            (45, Fq::from(3_u64)),
            (46, Fq::from(4_u64)),
            (47, Fq::from(95_u64)),
            (52, Fq::from(7_u64)),
            (53, Fq::from(7_u64)),
            (54, Fq::from(8_u64)),
            (55, Fq::from(99_u64)),
            (60, Fq::from(4_u64)),
            (61, Fq::from(4_u64)),
            (62, Fq::from(5_u64)),
            (63, Fq::from(96_u64)),
        ];

        let f1 = SparseMultilinearExtension::from_evaluations(
            6,
            &points
        );

        let g = vec![Fq::from(80u64), Fq::from(6u64)];
        let u = vec![Fq::from(27u64), Fq::from(12u64)];

        let actual = start_testing_phase2::<Fq>(
            &f1, &g, &u
        );

        // expected value
        let expected = vec![
            Fq::from(27_u64),
            Fq::from(54_u64),
            Fq::from(42_u64),
            Fq::from(69_u64),
        ];

        assert_eq!(expected, actual);
    }

    #[test]
    fn test_sumcheck_basic() {
        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);

        let g = vec![Fq::from(2u64), Fq::from(3u64)];
        let f2 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        run_sumcheck_protocol(f1, f2, f3, &g);   
    }

    #[test]
    fn tests_sumcheck_random() {
        let mut test_rng = ark_std::test_rng();
        // degree_x is the degree of both f2 (in x) and f3 (in y)
        // only go up to degree 7, otherwise tests start taking longer than that
        for degree_x in 3..8 {
            // degree of f1
            let d = 3 * degree_x;
            // cook a random sparse MLE with 2^(d - 2) evaluation points (25% density)
            let num_evaluations: usize = 1 << (d - 2);
            let points = (0..num_evaluations)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f1 = SparseMultilinearExtension::from_evaluations(
                d,
                &points);
            // cook a random g
            let g = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            let f2 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            let f3 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );

            run_sumcheck_protocol(f1, f2, f3, &g);
        }
    }

    #[test]
    fn tests_sumcheck_simple_combination() {
        // same parameters as the test tests_sumcheck, using the linear
        // combination, first with the coefficient (1, 0) and (0, 1)

        let points = vec![(8, Fq::from(1_u64)), (32, Fq::from(2_u64))];

        // assume |x| = |y| = |z| = 2 so f1 would be a 6-variate
        let f1 = SparseMultilinearExtension::from_evaluations(6, &points);
        let g = vec![Fq::from(2u64), Fq::from(3u64)];
        let f2 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );

        let f3 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![1, 6, 5, 0]
                .into_iter()
                .map(|x| Fq::from(x as u64))
                .collect(),
        );


        let zero_vector = vec![Fq::from(0_u64); g.len()];
        let ones_vector = vec![Fq::from(0_u64); g.len()];

        run_sumcheck_protocol_combined(
            f1.clone(),
            f2.clone(),
            f3.clone(),
            &g,
            &zero_vector,
            Fq::from(1u64),
            Fq::from(0_u64),
        );

        run_sumcheck_protocol_combined(
            f1,
            f2,
            f3,
            &ones_vector,
            &g,
            Fq::from(0_u64),
            Fq::from(1u64),
        );
    }

    #[test]
    fn tests_sumcheck_randomized_combination(){
        let mut test_rng = ark_std::test_rng();
        // degree x is the degree of both f2 and f3
        for degree_x in 3..8{
            // total degree in f1
            let d = 3 * degree_x;
            // cook a random sparse MLE with 2^(d-2) evaluation points (25% density)
            let num_evals: usize = 1 << (d - 2);
            let points = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f1 = SparseMultilinearExtension::from_evaluations(
                d,
                &points
            );
            // generate a random g1
            let g1 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate a random g2
            let g2 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate a random f2
            let f2 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),   
            );
            // generate a random f3
            let f3 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),   
            );

            run_sumcheck_protocol_combined(
                f1.clone(),
                f2.clone(),
                f3.clone(),
                &g1,
                &g2,
                Fq::rand(&mut test_rng),
                Fq::rand(&mut test_rng),
            );
        }
    }

    #[test]
    fn tests_sumcheck_randomized_combination_multiproof() {
        let mut test_rng = ark_std::test_rng();
        // degree_x is the degree of both f2 and f3
        // only climbing upto degree 7
        for degree_x in 3..8 {
            // total degree of f1
            let d = 3 * degree_x;
            // cook a random sparse MLE with 2^(d-2) evaluation points (25% density)
            let num_evals: usize = 1 << (d - 2);
            let points = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f11 = SparseMultilinearExtension::from_evaluations(
                d,
                &points,
            );
            let points2 = (0..num_evals)
                .map(|i| (i, Fq::rand(&mut test_rng)))
                .collect::<Vec<_>>();
            let f12 = SparseMultilinearExtension::from_evaluations(
                d,
                &points,
            );
            // generate random g1
            let g1 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generate random g2
            let g2 = (0..degree_x)
                .map(|_| Fq::rand(&mut test_rng))
                .collect::<Vec<_>>();
            // generating first f2 fron the random coefficients
            let f21 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            // generating first f2 fron the random coefficients
            let f31 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            let f22 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );
            let _f32 = DenseMultilinearExtension::from_evaluations_vec(
                degree_x,
                (0..(1 << degree_x))
                    .map(|_| Fq::rand(&mut test_rng))
                    .collect::<Vec<_>>(),
            );

            let sum_of_products = GKRProverSumOfProduct {
                elements: vec![
                    MLEProduct(f11.clone(), f21.clone(), f31.clone()),
                    MLEProduct(f12.clone(), f22.clone(), _f32.clone()),
                ],
            };

            run_sumcheck_protocol_combined_multiproof(
                sum_of_products,
                &g1,
                &g2,
                Fq::rand(&mut test_rng),
                Fq::rand(&mut test_rng),
            );

        }
    }
}