#[cfg(test)]
mod tests {
    use crate::{
        views::{MLEProduct, GKRProverSumOfProduct},
        side::{
            phase1_start as p1_org,
            phase2_start as p2_org, run_sumcheck_protocol,
            run_sumcheck_protocol_combined, run_sumcheck_protocol_combined_multiproof
        },
        helpers::to_zxy_form,
        Fq,
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
}