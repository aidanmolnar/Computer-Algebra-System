use std::{cmp::Reverse, collections::BTreeMap};

use super::generated_group;
use crate::{
    domain::{Domain, Monoid},
    integer::UInt64,
    Atom,
};

#[derive(Clone, Eq, PartialEq, Debug)]
// TODO: Should be generic over exponent domain
pub struct Mono(BTreeMap<Atom, u64>);

#[derive(PartialEq, Clone, Debug, Default)]
pub struct MonoDomain<E> {
    pub exp_domain: E,
}

// TODO: This should be generic over exponent group
impl Domain for MonoDomain<UInt64> {
    type Element = Mono;
}

impl Monoid for MonoDomain<UInt64> {
    fn identity(&self) -> Self::Element {
        Mono(BTreeMap::new())
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        generated_group::op_assign(&self.exp_domain, &mut lhs.0, rhs.0);
    }
}

// TODO: Allow monomials to be group if exponent domain is also a group?
// impl Group for MonoDomain<Int64> {
//     fn invert(&self, element: &mut Self::Element) {
//         generated_group_ops::invert(&self.exponent_group, &mut element.0);
//     }
// }

// Order monomials using lexographical order
impl PartialOrd for Mono {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(
            self.0
                .iter()
                .map(|(a, b)| (Reverse(a), b))
                .cmp(other.0.iter().map(|(a, b)| (Reverse(a), b)))
                .reverse(),
        )
    }
}

impl Ord for Mono {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl Mono {
    pub fn one() -> Self {
        MonoDomain::<UInt64>::default().identity()
    }

    pub fn is_one(&self) -> bool {
        self.0.is_empty()
    }

    pub fn total_deg(&self) -> u64 {
        self.0.values().sum()
    }

    pub fn power_of(&self, var: &Atom) -> u64 {
        self.0.get(var).cloned().unwrap_or(0)
    }

    pub fn vars(&self) -> impl Iterator<Item = &Atom> {
        self.0.keys()
    }

    pub fn remove_var(&mut self, var: &Atom) -> u64 {
        self.0.remove(var).unwrap_or(0)
    }

    // Constructs a monomial where var is raised to the power provided
    // TODO: This shares logic with Term for Polynomial...
    pub fn power(var: Atom, pow: u64) -> Mono {
        if pow == 0 {
            // Protects invariant
            Mono(BTreeMap::new())
        } else {
            Mono(BTreeMap::from([(var.clone(), pow)]))
        }
    }

    pub fn powers<const N: usize>(vars: [&'static str; N], pows: [u64; N]) -> Self {
        Mono(
            vars.into_iter()
                .zip(pows)
                .filter_map(|(var, pow)| {
                    if pow == 0 {
                        None
                    } else {
                        Some((Atom::from(var), pow))
                    }
                })
                .collect(),
        )
    }

    pub fn try_divide(self, other: Self) -> Option<Self> {
        let mut q = self;

        for (atom, divisor_pow) in other.0 {
            if let Some(dividend_pow) = q.0.get_mut(&atom) {
                if divisor_pow > *dividend_pow {
                    // To many of symbol in divisor: indivisible
                    return None;
                } else if divisor_pow == *dividend_pow {
                    // Exactly as many in divisor: remove symbol
                    q.0.remove(&atom);
                } else {
                    // More in dividend then divisor: subtract any from divisor
                    *dividend_pow -= divisor_pow;
                }
            } else {
                // Divisor has atoms that dividend doesn't
                return None;
            }
        }

        Some(q)
    }

    // TODO: Maybe this should not consume self...
    pub fn try_into_univariate(self) -> Result<Option<(Atom, u64)>, ()> {
        let mut mono = self.0.into_iter();
        let res = mono.next();

        if mono.next().is_some() {
            Err(())
        } else {
            Ok(res)
        }
    }
}

#[test]
fn test_mono_ordering() {
    use crate::polynomial::poly::Poly;

    // https://en.wikipedia.org/wiki/Monomial_order#Lexicographic_order
    let p = Poly::sparse(
        ["x2", "x1"],
        [
            (1, [0, 2]),
            (2, [1, 1]),
            (3, [0, 1]),
            (4, [2, 0]),
            (5, [1, 0]),
            (6, [0, 0]),
        ],
    );

    assert!(p.into_terms().map(|x| x.coef).eq(1..=6));
}
