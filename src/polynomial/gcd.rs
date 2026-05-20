use std::{cmp::Reverse, collections::BTreeSet};

use crate::{
    domain::{euclidean_algorithm, NormalForm, Ring},
    integer::Int64,
    modular::{chinese_remainder, Int64ModP, Modular},
    polynomial::{
        poly::{Poly, PolyDomain},
        primes::PRIMES,
        unipoly::{UniPoly, UniPolyDomain},
    },
    Atom,
};

// Algorithms for Computer Algebra (page 318 of pdf)
// Adapted for univariate case
// TODO: Would be really cool if we could combine with pgcd, or at least break out shared logic
pub fn mgcd(mut a: Poly<i64>, mut b: Poly<i64>) -> Poly<i64> {
    let mut vars = BTreeSet::new();
    vars.extend(a.vars());
    vars.extend(b.vars());
    let vars: Vec<_> = vars.into_iter().collect();

    if a.is_zero() {
        return b;
    } else if b.is_zero() {
        return a;
    }

    // Handle this case directly
    if vars.len() == 0 {
        let a = a.into_terms().next().map(|term| term.coef).unwrap_or(0);
        let b = b.into_terms().next().map(|term| term.coef).unwrap_or(0);

        return Poly::dense("", vec![euclidean_algorithm(&Int64, a, b)]);
    }

    // TODO: Better name? Should this be an input parameter?
    let base_ring = PolyDomain { coef_domain: Int64 };

    let a_cont = a.remove_content(&base_ring.coef_domain);
    let b_cont = b.remove_content(&base_ring.coef_domain);

    let lc_a = *a.leading_term().unwrap().coef;
    let lc_b = *b.leading_term().unwrap().coef;

    let c = euclidean_algorithm(&Int64, a_cont, b_cont);
    let g = euclidean_algorithm(&Int64, lc_a, lc_b);

    let mut q = Int64ModP::new(0);
    let mut h = Poly::zero();

    let mut n = None;

    // Use precomputed primes
    let mut prime_idx = 0;
    let mut get_prime = || {
        let p = PRIMES[prime_idx];
        prime_idx += 1;
        p as i64
    };

    loop {
        // Get a prime that doesn't divide g
        let mut p = get_prime();
        while g % p == 0 || lc_a % p == 0 {
            p = get_prime();
        }

        let ring = PolyDomain {
            coef_domain: Int64ModP::new(p),
        };

        // Construct image of a,b, and g
        let a_p = a.clone().to_modp(p);
        let b_p = b.clone().to_modp(p);
        let g_p = g % p;

        let Some(mut c_p) = pgcd(&ring, &vars, a_p, b_p) else {
            // Restart
            q = Int64ModP::new(0);
            h = Poly::zero();
            n = None;
            continue;
        };

        // Normalize so that g_p = lc(Cp)
        c_p = ring.normal(c_p);
        c_p.ring_mul_by(&ring.coef_domain, g_p);

        // TODO: We need to do reverse because our monomial ordering is backwards
        let m = Some(Reverse(c_p.leading_term().unwrap().mono.clone()));

        // if h.is_zero() || m < n {
        if h.is_zero() || m < n {
            // Unlocky homomorphisms, restart
            q = Int64ModP::new(p);
            h = c_p;
            n = m;
        } else if m == n {
            (q, h) = chinese_remainder(vec![(q, h), (Int64ModP::new(p), c_p)]);
            // Symmetric range!
            h.mods(q.p);
        }

        let lc_h = h.leading_term().unwrap().coef;
        if lc_h == &g {
            // Update coefficients of gcd candidate
            let mut pp_h = h.clone();
            pp_h.remove_content(&base_ring.coef_domain);

            // TODO: We could directly return a/g and b/g here to save work for consumer
            if base_ring.divides(a.clone(), pp_h.clone())
                && base_ring.divides(b.clone(), pp_h.clone())
            {
                return Poly::from_coef(c) * pp_h;
            }
        }
    }
}

// Useful pseudocode: https://www.cecm.sfu.ca/CAG/theses/suling.pdf
//
// NOTE: We don't directly return the recursive representation, because it has to be collapsed anyway
// Base case Z_p[x]
// First interpolations Z_p[y][x]
// Collase to Z_p[y,x]
// Second interpolations Z_p[z][y,x]
// Etc.
fn pgcd(
    ring: &PolyDomain<Int64ModP>,
    vars: &[Atom],
    a: Poly<i64>,
    b: Poly<i64>,
) -> Option<Poly<i64>> {
    let k = vars.len();

    // Recursive polynomial ring Z_p[x][y,z,...]
    let rec_ring = PolyDomain {
        coef_domain: UniPolyDomain {
            coef_domain: ring.coef_domain.clone(),
        },
    };

    // Base case, polynomials are univariate
    if k <= 1 {
        // Run pgcd, assumed  univariate case
        let c = euclidean_algorithm(
            &rec_ring.coef_domain,
            a.try_into().expect("k <= 1"),
            b.try_into().expect("k <= 1"),
        );

        if c.is_zero() {
            return Some(ring.one());
        } else {
            return Some(c.to_poly(&rec_ring.coef_domain.coef_domain, vars[0].clone()));
        }
    }

    // Variable we interpolate over
    let var = vars[k - 1].clone();

    // Convert a and b to recursive representation over var
    let mut a = a.clone().move_var_to_coefs(&rec_ring.coef_domain, &var);
    let mut b = b.clone().move_var_to_coefs(&rec_ring.coef_domain, &var);

    let a_cont = a.remove_content(&rec_ring.coef_domain);
    let b_cont = b.remove_content(&rec_ring.coef_domain);

    // Compute gcd using euclidean algorithm
    let c = euclidean_algorithm(&rec_ring.coef_domain, a_cont.clone(), b_cont.clone());

    let lc_a = a.leading_term().unwrap().coef.clone();
    let lc_b = b.leading_term().unwrap().coef.clone();

    let g = euclidean_algorithm(&rec_ring.coef_domain, lc_a.clone(), lc_b.clone());

    let mut q = Modular {
        domain: rec_ring.coef_domain.clone(),
        p: rec_ring.coef_domain.one(),
    };
    let mut h = rec_ring.zero();

    let mut n = None;

    // NOTE: DO NOT USE BETA = 0! It's reducible?
    let mut beta = 0; // Equal to 1 due to first statement in loop
    while beta < ring.coef_domain.p - 1 {
        beta += 1;
        // Get a new member of Z[p] (where member is x - beta) where g(beta) != 0
        let g_beta = g.clone().eval(&rec_ring.coef_domain.coef_domain, beta);
        if g_beta == 0 {
            continue;
        }

        let lc_a_beta = lc_a.eval(&rec_ring.coef_domain.coef_domain, beta);
        if lc_a_beta == 0 {
            continue;
        }

        let a_beta = a.clone().coef_eval(&rec_ring.coef_domain, beta);
        let b_beta = b.clone().coef_eval(&rec_ring.coef_domain, beta);

        let Some(c_beta) = pgcd(ring, &vars[0..(k - 1)], a_beta, b_beta) else {
            // Restart
            q = Modular {
                domain: rec_ring.coef_domain.clone(),
                p: rec_ring.coef_domain.one(),
            };
            h = rec_ring.zero();
            n = None;
            continue;
        };

        let mut c_beta = ring.normal(c_beta);
        c_beta.ring_mul_by(&ring.coef_domain, g_beta);

        let c_beta = c_beta.move_var_to_coefs(&rec_ring.coef_domain, &var);

        let beta_poly: UniPoly<_> = Poly::dense(var.clone(), vec![1, -beta])
            .to_modp(ring.coef_domain.p)
            .try_into()
            .expect("constructed from a dense univariate polynomial");

        let new = Modular {
            domain: rec_ring.coef_domain.clone(),
            p: beta_poly.clone(),
        };

        // TODO: We need to do reverse because our monomial ordering is backwards
        let m = Some(Reverse(c_beta.leading_term().unwrap().mono.clone()));

        // Test for unlucky homomorphism
        if h.is_zero() || m < n {
            n = m;
            q = new;
            h = c_beta;
        } else if m == n {
            // Build h via polynomial interpolation
            (q, h) = chinese_remainder(vec![(q, h), (new, c_beta)]);
        }

        let lc_h = h.leading_term().unwrap().coef;
        if lc_h == &g {
            let mut pp_h = h.clone();
            pp_h.remove_content(&rec_ring.coef_domain);

            // Check divisibility
            if rec_ring.divides(a.clone(), pp_h.clone())
                && rec_ring.divides(b.clone(), pp_h.clone())
            {
                let mut result = pp_h;
                result.ring_mul_by(&rec_ring.coef_domain, c);

                // Collapse result back into a dense polynomial
                return Some(result.collapse(&ring, var));
            }
        }
    }

    // Exhausted possible beta values, ring is too small
    return None;
}

#[cfg(test)]
#[allow(clippy::zero_prefixed_literal)]
// TODO: We can probably remove the non proptest tests, they were mostly useful for debugging during development
mod tests {

    fn pgcd_test_case() -> (Poly<i64>, Poly<i64>) {
        let a = Poly::sparse(
            ["x", "y", "z"],
            [
                //
                (0009, [5, 0, 0]),
                (0002, [4, 1, 1]),
                (-189, [3, 3, 1]),
                (0117, [3, 1, 2]),
                (0003, [3, 0, 0]),
                (-042, [2, 4, 2]),
                (0026, [2, 2, 3]),
                //
                (0018, [2, 0, 0]),
                (-063, [1, 3, 1]),
                (0039, [1, 1, 2]),
                (0004, [1, 1, 1]),
                (0006, [0, 0, 0]),
            ],
        );
        let b = Poly::sparse(
            ["x", "y", "z"],
            [
                //
                (0006, [6, 0, 0]),
                (-126, [4, 3, 1]),
                (0078, [4, 1, 2]),
                (0001, [4, 1, 0]),
                (0001, [4, 0, 1]),
                (0013, [3, 0, 0]),
                (-021, [2, 4, 1]),
                (-021, [2, 3, 2]),
                //
                (0013, [2, 2, 2]),
                (0013, [2, 1, 3]),
                (-021, [1, 3, 1]),
                (0013, [1, 1, 2]),
                (0002, [1, 1, 0]),
                (0002, [1, 0, 1]),
                (0002, [0, 0, 0]),
            ],
        );

        (a, b)
    }

    #[test]
    fn pgcd_mod11() {
        let ring = PolyDomain {
            coef_domain: Int64ModP::new(11),
        };

        let (a, b) = pgcd_test_case();
        let expected = Poly::sparse(
            ["x", "y", "z"],
            [
                (1, [3, 0, 0]),
                (1, [1, 3, 1]),
                (2, [1, 1, 2]),
                (2, [0, 0, 0]),
            ],
        )
        .to_modp(11);
        assert!(
            expected
                == dbg!(pgcd(
                    &ring,
                    &[Atom::from("x"), Atom::from("y"), Atom::from("z")],
                    a.clone().to_modp(11),
                    b.clone().to_modp(11)
                )
                .unwrap())
        );
    }

    #[test]
    fn pgcd_mod13() {
        let p = 13;
        let ring = PolyDomain {
            coef_domain: Int64ModP::new(p),
        };

        let (a, b) = pgcd_test_case();
        let mut res = pgcd(
            &ring,
            &[Atom::from("x"), Atom::from("y"), Atom::from("z")],
            a.clone().to_modp(p),
            b.clone().to_modp(p),
        )
        .unwrap();
        res.ring_mul_by(&ring.coef_domain, 3); // Book only gives result after normalizing

        let expected = Poly::sparse(
            ["x", "y", "z"],
            [(3, [3, 0, 0]), (2, [1, 3, 1]), (6, [0, 0, 0])],
        )
        .to_modp(p);
        assert!(expected == res);
    }

    #[test]
    fn pgcd_mod17() {
        let p = 17;
        let ring = PolyDomain {
            coef_domain: Int64ModP::new(p),
        };

        let (a, b) = pgcd_test_case();
        let mut res = pgcd(
            &ring,
            &[Atom::from("x"), Atom::from("y"), Atom::from("z")],
            a.clone().to_modp(p),
            b.clone().to_modp(p),
        )
        .unwrap();
        res.ring_mul_by(&ring.coef_domain, 3); // Book only gives result after normalizing

        let expected = Poly::sparse(
            ["x", "y", "z"],
            [
                (3, [3, 0, 0]),
                (5, [1, 3, 1]),
                (5, [1, 1, 2]),
                (6, [0, 0, 0]),
            ],
        )
        .to_modp(p);
        assert!(expected == res);
    }

    #[test]
    fn mgcd_multi() {
        let (a, b) = pgcd_test_case();
        let expected = Poly::sparse(
            ["x", "y", "z"],
            [
                (1, [3, 0, 0]),
                (-21, [1, 3, 1]),
                (13, [1, 1, 2]),
                (2, [0, 0, 0]),
            ],
        );
        let res = mgcd(a.clone(), b.clone());

        assert!(expected == res);
    }

    #[test]
    fn uni_mgcd_monagan() {
        // based on https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/lec8/Lec8A.pdf
        fn check_gcd(g: Poly<i64>) {
            let a = g.clone() * Poly::dense("x", vec![5, 18]);
            let b = g.clone() * Poly::dense("x", vec![5, 1]);
            assert!(mgcd(a, b) == g);
        }
        check_gcd(Poly::dense("x", vec![13, -11]));
        check_gcd(Poly::dense("x", vec![13, 11]));
    }

    #[test]
    fn uni_mgcd_winkler() {
        // Example 4.2.1 on page 94 of Polynomial Algorithms in Computer Algebra

        let a = Poly::dense("x", vec![2, -13, 20, 12, -20, -15, -18]);
        let b = Poly::dense("x", vec![2, 1, -14, -11, 22, 28, 8]);

        assert!(mgcd(a, b) == Poly::dense("x", vec![1, -1, -2]));
    }

    use proptest::prelude::*;

    use crate::{
        domain::Ring,
        integer::Int64,
        modular::Int64ModP,
        polynomial::{gcd::pgcd, mgcd, PolyDomain},
        Atom, Poly,
    };

    fn uni_poly_strategy() -> impl Strategy<Value = Poly<i64>> {
        let coeffs = prop::collection::vec(-50i64..=50, 1..10);
        coeffs.prop_map(|mut c| {
            c[0] = 1;
            Poly::dense("x", c)
        })
    }

    proptest! {
        #[test]
        fn uni_proptest(p1 in uni_poly_strategy(), p2 in uni_poly_strategy(), g in uni_poly_strategy()) {
            let domain = PolyDomain{coef_domain: Int64};
            if mgcd(p1.clone(), p2.clone()) != domain.one() {
                return Err(TestCaseError::Reject("g will not be gcd".into()));
            }

            let a = g.clone() * p1;
            let b = g.clone() * p2;

            let d = mgcd(a, b);

            assert_eq!(g, d);
        }
    }

    fn multi_poly_strategy() -> impl Strategy<Value = Poly<i64>> {
        let vars = ["x", "y", "z"];
        // Each exponent vector has length 3 (for x, y, z)
        let term_strategy = (-20i64..=20, prop::array::uniform3(0..4u64));
        let terms_strategy = prop::collection::vec(term_strategy, 1..5);

        terms_strategy.prop_map(move |terms| Poly::sparse(vars, terms))
    }

    proptest! {
        #[test]
        fn multi_proptest(p1 in multi_poly_strategy(), p2 in multi_poly_strategy(), mut g in multi_poly_strategy()) {


            let domain = PolyDomain{coef_domain: Int64};
            if mgcd(p1.clone(), p2.clone()) != domain.one() || g.is_zero() {
                return Err(TestCaseError::Reject("g will not be gcd".into()));
            }

            g.force_normalize_leading_coef();

            let a = g.clone() * p1;
            let b = g.clone() * p2;

            let d = mgcd(a, b);

            assert_eq!(g, d);
        }
    }
}
