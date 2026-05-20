use crate::{
    domain::{Domain, Group, Monoid, NormalForm, Ring},
    integer::Int64,
    polynomial::{mgcd, Poly, PolyDomain, Term},
};

use super::Atom;

// TODO: We may want to store C instead of Poly<C> and then we also generalize to other domains
//       We would need a gcd trait (UFD?) for that.
#[derive(Debug, Clone, PartialEq)]
pub struct RationalFunc<C> {
    num: Poly<C>,
    den: Poly<C>,
}

pub struct RationalFuncDomain<C> {
    pub poly_domain: PolyDomain<C>,
}

impl RationalFunc<i64> {
    pub fn one() -> Self {
        Self {
            num: Poly::one(),
            den: Poly::one(),
        }
    }

    pub fn zero() -> Self {
        Self {
            num: Poly::zero(),
            den: Poly::one(),
        }
    }

    pub fn dense(var: impl Into<Atom>, coefs: Vec<i64>) -> Self {
        Self {
            num: Poly::dense(var, coefs),
            den: Poly::one(),
        }
    }

    // TODO: Rename?
    // pub fn power(coef: i64, var: Atom, pow: u64) -> Self {
    //     let mono = Mono::power(var, pow);
    //     let num = PolyDomain { coef_domain: Int64 }.term_to_element(Term { coef, mono });
    //     Self {
    //         num,
    //         den: Poly::one(),
    //     }
    // }
    pub fn from_term(term: Term<i64>) -> Self {
        Self {
            num: PolyDomain { coef_domain: Int64 }.term_to_element(term),
            den: Poly::one(),
        }
    }

    pub fn into_parts(self) -> (Poly<i64>, Poly<i64>) {
        (self.num, self.den)
    }

    pub fn from_poly(poly: Poly<i64>) -> Self {
        Self {
            num: poly,
            den: Poly::one(),
        }
    }
}

// TODO: Generalize pgcd over arbitrary field so that this can have various coefficients
impl RationalFuncDomain<Int64> {
    // TODO: jenk
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            poly_domain: PolyDomain { coef_domain: Int64 },
        }
    }

    // TODO: Simplify by returning cofactors directly from mgcd
    // TODO: if the denominator is equal to zero, this divides out by a pole.  We may want to not allow that?
    pub fn from_frac(&self, num: Poly<i64>, den: Poly<i64>) -> RationalFunc<i64> {
        let g = mgcd(num.clone(), den.clone());

        let (mut num, r) = self.poly_domain.multi_poly_div(num, g.clone());
        assert!(r.is_zero());
        let (mut den, r) = self.poly_domain.multi_poly_div(den, g);
        assert!(r.is_zero());

        // Normalize by making numerator carry negatives
        let u = self
            .poly_domain
            .coef_domain
            .inverse_leading_unit(den.leading_term().unwrap().coef);

        // Multiplying num and den by a value is the same as multiplying by one
        num.ring_mul_by(&self.poly_domain.coef_domain, u);
        den.ring_mul_by(&self.poly_domain.coef_domain, u);

        RationalFunc { num, den }
    }
}

impl<C: Domain> Domain for RationalFuncDomain<C> {
    type Element = RationalFunc<C::Element>;
}

impl Monoid for RationalFuncDomain<Int64> {
    fn identity(&self) -> Self::Element {
        RationalFunc::one()
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        // a/b + c/d = (a*d + c*b)/(b*d)
        let num = lhs.num.clone() * rhs.den.clone() + rhs.num * lhs.den.clone();
        let den = lhs.den.clone() * rhs.den;
        *lhs = self.from_frac(num, den)
    }
}

impl Group for RationalFuncDomain<Int64> {
    fn invert(&self, element: &mut Self::Element) {
        self.poly_domain.invert(&mut element.num)
    }
}

impl Ring for RationalFuncDomain<Int64> {
    fn one(&self) -> Self::Element {
        Self::Element::one()
    }

    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        self.from_frac(lhs.num * rhs.num, lhs.den * rhs.den)
    }
}

// TODO: impl Field for RationalFuncDomain?

// TODO: We probably don't want these implemented directly on Rational...
mod ops {
    use std::ops::{Add, Div, Mul, Sub};

    use crate::domain::Ring;

    use super::{RationalFunc, RationalFuncDomain};

    impl Add for RationalFunc<i64> {
        type Output = Self;

        fn add(self, rhs: Self) -> Self::Output {
            RationalFuncDomain::new().add(self, rhs)
        }
    }

    impl Sub for RationalFunc<i64> {
        type Output = Self;

        fn sub(self, rhs: Self) -> Self::Output {
            RationalFuncDomain::new().sub(self, rhs)
        }
    }

    impl Mul for RationalFunc<i64> {
        type Output = Self;

        fn mul(self, rhs: Self) -> Self::Output {
            RationalFuncDomain::new().mul(self, rhs)
        }
    }

    impl Div for RationalFunc<i64> {
        type Output = Self;

        fn div(self, rhs: Self) -> Self::Output {
            RationalFuncDomain::new().from_frac(self.num * rhs.den, self.den * rhs.num)
        }
    }
}

#[test]
fn rational_test() {
    let a = (
        //
        RationalFunc::dense("x", vec![1, 2, 0])
            + RationalFunc::dense("y", vec![2, 0])
            + RationalFunc::dense("y", vec![1, 0]) * RationalFunc::dense("x", vec![1, 0])
    ) / RationalFunc::dense("x", vec![1, 2])
        - RationalFunc::dense("y", vec![1, 0]);
    assert!(dbg!(a) == RationalFunc::dense("x", vec![1, 0]));
}
