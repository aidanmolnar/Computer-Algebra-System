use std::collections::BTreeMap;

use crate::{
    domain::{Domain, EuclideanDomain, Field, Group, Monoid, NormalForm, Ring},
    integer::UInt64,
    polynomial::{generated_group, Mono, Poly, Term},
    Atom,
};

#[derive(Clone, Debug, PartialEq)]
pub struct UniPolyDomain<C> {
    pub coef_domain: C,
}

#[derive(Clone, Debug, PartialEq, Default)]
pub struct UniPoly<C>(BTreeMap<u64, C>);

impl<C> TryFrom<Poly<C>> for UniPoly<C> {
    type Error = ();

    fn try_from(value: Poly<C>) -> Result<Self, ()> {
        let mut new = BTreeMap::new();
        let mut var = None;

        for Term { mono, coef } in value.into_terms() {
            let pow = if let Some((atom, pow)) = mono.try_into_univariate()? {
                if let Some(exist) = &var {
                    if &atom != exist {
                        return Err(());
                    }
                } else {
                    var = Some(atom);
                }
                pow
            } else {
                0
            };
            new.insert(pow, coef);
        }

        Ok(Self(new))
    }
}

impl<C> UniPoly<C> {
    // TODO: Maybe move this function to the domain?
    pub fn eval<R>(&self, domain: &R, alpha: C) -> C
    where
        C: Clone,
        R: Ring<Element = C>,
    {
        // It's a bit of a bummer we need to reallocate here, but we need to modify all the keys of the underlying btreemap
        let mut new = domain.zero();

        for (exp, mut coef) in self.0.iter().map(|(a, b)| (a.clone(), b.clone())) {
            // TODO: Simplify manually rolled exponentiation, maybe this could be a function on ring?
            for _ in 0..exp {
                coef = domain.mul(coef, alpha.clone());
            }

            new = domain.add(new, coef);
        }
        new
    }

    pub fn deg(&self) -> u32 {
        self.0.iter().map(|(a, _)| *a as u32).max().unwrap_or(0)
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub fn leading_coef(&self) -> Option<&C> {
        self.0.last_key_value().map(|term| term.1)
    }

    pub fn to_poly<R>(self, ring: &R, var: Atom) -> Poly<C>
    where
        R: Ring<Element = C>,
        C: PartialEq,
    {
        Poly::from_map(
            ring,
            self.0
                .into_iter()
                .map(|(pow, coef)| (Mono::power(var.clone(), pow), coef))
                .collect(),
        )
    }

    pub fn from_power<R>(ring: &R, power: u64, coef: C) -> Self
    where
        R: Ring<Element = C>,
        C: PartialEq,
    {
        if coef != ring.zero() {
            Self(BTreeMap::from([(power, coef)]))
        } else {
            Self(BTreeMap::new())
        }
    }
}

impl<C: Domain> Domain for UniPolyDomain<C> {
    type Element = UniPoly<C::Element>;
}

impl<C: Group> Monoid for UniPolyDomain<C> {
    fn identity(&self) -> Self::Element {
        UniPoly(BTreeMap::new())
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        // TODO: We need to check that we're adding univariate polynomials of the same variable!!
        generated_group::op_assign(&self.coef_domain, &mut lhs.0, rhs.0);
    }
}

impl<C: Group> Group for UniPolyDomain<C> {
    fn invert(&self, element: &mut Self::Element) {
        generated_group::invert(&self.coef_domain, &mut element.0);
    }
}

impl<C> Ring for UniPolyDomain<C>
where
    C: Ring,
    C::Element: Clone,
{
    fn one(&self) -> Self::Element {
        UniPoly(BTreeMap::from([(0, self.coef_domain.one())]))
    }

    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        UniPoly(generated_group::mul(
            &self.coef_domain,
            &UInt64,
            lhs.0,
            rhs.0,
        ))
    }
}

impl<C> EuclideanDomain for UniPolyDomain<C>
where
    C: Field,
    C::Element: Clone,
{
    // Algorithm 2.5 on page 39 of Modern Computer Algebra
    // Univariate polynomial division with remainder
    fn div(
        &self,
        dividend: Self::Element,
        divisor: Self::Element,
    ) -> (Self::Element, Self::Element) {
        assert!(!divisor.is_zero(), "Can't divide by zero");

        // Leading coefficient of divisor
        let c = divisor.leading_coef().unwrap();
        let u = self.coef_domain.mul_inverse(c.clone());

        // Input degrees
        let n = dividend.deg();
        let m = divisor.deg();

        let mut q = self.zero(); // TODO: Make a UniPoly::zero function
        let mut r = dividend;

        if m > n {
            return (q, r);
        }

        for i in (0..=n - m).rev() {
            if r.deg() == m + i {
                let lc_r = r.leading_coef().unwrap();
                let coef = self.coef_domain.mul(lc_r.clone(), u.clone());

                let q_i = UniPoly(BTreeMap::from([(i as u64, coef)]));

                q = self.add(q, q_i.clone());
                r = self.sub(r, self.mul(q_i, divisor.clone()));

                // TODO: Shouldn't be necessary? is there a bug
                if r.is_zero() {
                    break;
                }
            }
        }

        return (q, r);
    }

    fn euclidean_function(&self, element: &Self::Element) -> u32 {
        element.deg()
    }
}

impl<C> NormalForm for UniPolyDomain<C>
where
    C: Field,
    C::Element: Clone,
{
    fn normal(&self, element: Self::Element) -> Self::Element {
        self.mul(self.inverse_leading_unit(&element), element)
    }

    fn inverse_leading_unit(&self, element: &Self::Element) -> Self::Element {
        let c = if let Some(coef) = element.leading_coef() {
            self.coef_domain.mul_inverse(coef.clone())
        } else {
            self.coef_domain.one()
        };
        assert!(self.coef_domain.zero() != c);
        UniPoly(BTreeMap::from([(0, c)]))
    }
}

#[test]
fn test_leading_coef() {
    let poly: UniPoly<_> = Poly::dense("x", vec![1, 2, 3, 4, 5]).try_into().unwrap();
    assert!(*poly.leading_coef().unwrap() == 1);
}

#[test]
fn test_eval() {
    use crate::integer::Int64;

    let poly: UniPoly<_> = Poly::one().try_into().unwrap();
    let evaled = poly.eval(&Int64, 10);
    assert!(evaled == 1);
}

// TODO: Add test for univariate division
// Compare to multi_poly_div.
// Would be easier to write with a generic field type, or we assume that leading coefficient is a unit?
