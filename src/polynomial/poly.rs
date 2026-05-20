use std::collections::{BTreeMap, HashSet};

use crate::{
    domain::{
        euclidean_algorithm, Domain, EuclideanDomain, Field, Group, Monoid, NormalForm, Ring,
    },
    integer::UInt64,
    polynomial::mono::{Mono, MonoDomain},
    Atom,
};

use super::generated_group;

#[derive(PartialEq, Clone, Debug, Default)]
pub struct PolyDomain<C> {
    pub coef_domain: C,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Poly<C>(BTreeMap<Mono, C>);

// TODO: Would be nice to make these more general where possible
impl Poly<i64> {
    // TODO: Maybe a better way to do this?
    pub fn to_modp(mut self, p: i64) -> Poly<i64> {
        self.modp(p);
        Poly(self.0)
    }

    pub fn modp(&mut self, p: i64) {
        self.0.retain(|_, coef| {
            *coef = coef.rem_euclid(p);
            *coef != 0
        });
    }

    pub fn mods(&mut self, p: i64) {
        self.0.retain(|_, coef| {
            *coef = coef.rem_euclid(p);
            // Put in symmetric range
            if *coef > p / 2 {
                *coef -= p;
            }
            *coef != 0
        });
    }

    // Used in tests!
    pub fn dense(var: impl Into<Atom>, coefs: Vec<i64>) -> Self {
        let var = var.into();
        Poly(
            coefs
                .into_iter()
                .rev()
                .enumerate()
                .filter_map(|(pow, coef)| {
                    if coef == 0 {
                        None
                    } else {
                        Some((Mono::power(var.clone(), pow as u64), coef))
                    }
                })
                .collect(),
        )
    }

    pub fn sparse<const N: usize>(
        vars: [&'static str; N],
        terms: impl IntoIterator<Item = (i64, [u64; N])>,
    ) -> Self {
        Poly(
            terms
                .into_iter()
                .filter_map(|(coef, pows)| {
                    if coef == 0 {
                        None
                    } else {
                        Some((Mono::powers(vars, pows), coef))
                    }
                })
                .collect(),
        )
    }

    pub fn from_coef(coef: i64) -> Self {
        if coef == 0 {
            Poly::zero()
        } else {
            Poly(BTreeMap::from([(Mono::one(), coef)]))
        }
    }

    pub fn one() -> Self {
        Self::from_coef(1)
    }

    // Modern Computer Algebra (page 158)
    // ||f||_inf
    pub fn max_norm(&self) -> i64 {
        self.0.values().map(|x| x.abs()).max().unwrap_or(0)
    }

    // Modern Computer Algebra (page 165)
    // ||f||_1
    pub fn one_norm(&self) -> i64 {
        self.0.values().map(|x| x.abs()).sum()
    }

    // Only used during testing.  Sets the coefficient of leading term to 1
    pub fn force_normalize_leading_coef(&mut self) {
        *self.0.first_entry().unwrap().get_mut() = 1;
    }
}

impl<C> Poly<C> {
    pub fn iter_terms(&self) -> impl Iterator<Item = TermRef<'_, C>> {
        self.0.iter().map(|(mono, coef)| TermRef { mono, coef })
    }
    pub fn into_terms(self) -> impl Iterator<Item = Term<C>> {
        self.0.into_iter().map(|(mono, coef)| Term { mono, coef })
    }

    pub fn from_map<R>(coef_ring: &R, mut map: BTreeMap<Mono, C>) -> Self
    where
        R: Ring<Element = C>,
        C: PartialEq,
    {
        map.retain(|_, coef| *coef != coef_ring.zero());
        Self(map)
    }

    pub fn leading_term(&self) -> Option<TermRef<'_, C>> {
        self.0
            .first_key_value()
            .map(|(mono, coef)| TermRef { mono, coef })
    }

    pub fn deg_of(&self, var: &Atom) -> u32 {
        self.0.keys().map(|m| m.power_of(var)).max().unwrap_or(0) as u32
    }

    pub fn total_deg(&self) -> u32 {
        self.0.keys().map(|m| m.total_deg()).max().unwrap_or(0) as u32
    }

    pub fn vars(&self) -> HashSet<Atom> {
        let mut vars = HashSet::new();
        for m in self.0.keys() {
            vars.extend(m.vars().cloned())
        }

        vars
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    pub fn zero() -> Self {
        Poly(BTreeMap::new())
    }

    pub fn ring_mul_by<R>(&mut self, ring: &R, a: C)
    where
        R: Ring<Element = C>,
        C: Clone + PartialEq,
    {
        if a == ring.zero() {
            self.0 = BTreeMap::new();
            return;
        }

        for (_, existing) in self.0.iter_mut() {
            *existing = ring.mul(existing.clone(), a.clone());
        }
    }

    pub fn ring_modp<R>(&mut self, ring: &R, p: C)
    where
        R: EuclideanDomain<Element = C>,
        C: Clone + PartialEq,
    {
        self.0.retain(|_, coef| {
            // TODO: div_assign for EuclidaenDomain?
            *coef = ring.div(coef.clone(), p.clone()).1;
            // * coef = coef.rem_euclid(p);
            *coef != ring.zero()
        });
    }

    // Modern Computer Algebra (page 147)
    pub fn content<R>(&self, coefficient_domain: &R) -> C
    where
        R: EuclideanDomain<Element = C> + NormalForm,
        C: Clone + PartialEq,
    {
        let mut coefs = self.0.values();
        let Some(mut g) = coefs.next().cloned() else {
            return coefficient_domain.zero();
        };
        for a in coefs {
            g = euclidean_algorithm(coefficient_domain, g, a.clone());
        }
        g
    }

    // Returns the content and leaves the polynomial as just the primitive part
    pub fn remove_content<R>(&mut self, coefficient_domain: &R) -> C
    where
        R: EuclideanDomain<Element = C> + NormalForm,
        C: Clone + PartialEq + Default,
    {
        let content = self.content(coefficient_domain);

        for coef in self.0.values_mut() {
            let (q, r) = coefficient_domain.div(std::mem::take(coef), content.clone());
            assert!(r == coefficient_domain.zero());
            *coef = q;
        }

        content
    }
}

impl<C> Domain for PolyDomain<C>
where
    C: Domain,
{
    type Element = Poly<C::Element>;
}

impl<C> Monoid for PolyDomain<C>
where
    C: Group,
{
    fn identity(&self) -> Self::Element {
        Poly(BTreeMap::new())
    }
    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        generated_group::op_assign(&self.coef_domain, &mut lhs.0, rhs.0);
    }
}

impl<C> Group for PolyDomain<C>
where
    C: Group,
{
    fn invert(&self, element: &mut Self::Element) {
        generated_group::invert(&self.coef_domain, &mut element.0);
    }
}

impl<C> Ring for PolyDomain<C>
where
    C: Ring,
    C::Element: Clone,
{
    fn one(&self) -> Self::Element {
        Poly(BTreeMap::from([(Mono::one(), self.coef_domain.one())]))
    }

    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        Poly(generated_group::mul(
            &self.coef_domain,
            &MonoDomain::<UInt64>::default(),
            lhs.0,
            rhs.0,
        ))
    }
}

impl<C> NormalForm for PolyDomain<C>
where
    C: Field,
    C::Element: Clone,
{
    fn normal(&self, element: Self::Element) -> Self::Element {
        let inverse_leading = self.inverse_leading_unit(&element);
        self.mul(element, inverse_leading)
    }

    fn inverse_leading_unit(&self, element: &Self::Element) -> Self::Element {
        let coef = if let Some(term) = element.leading_term() {
            self.coef_domain.mul_inverse(term.coef.clone())
        } else {
            self.coef_domain.one()
        };

        self.term_to_element(Term {
            coef,
            mono: Mono::one(),
        })
    }
}

impl<C: Domain> PolyDomain<C> {
    // TODO: Rename
    pub fn term_to_element(&self, term: Term<C::Element>) -> Poly<C::Element>
    where
        C: Ring,
    {
        if term.coef == self.coef_domain.identity() {
            self.identity()
        } else {
            Poly(BTreeMap::from([(term.mono, term.coef)]))
        }
    }

    // TODO: Should this be a function on Poly?
    pub fn eval(&self, poly: Poly<C::Element>, var: &Atom, alpha: C::Element) -> Poly<C::Element>
    where
        C: Ring,
        C::Element: Clone,
    {
        // It's a bit of a bummer we need to reallocate here, but we need to modify all the keys of the underlying btreemap
        let mut new = Poly::zero();

        for (mut mono, mut coef) in poly.0 {
            let exp = mono.remove_var(var);
            // TODO: Simplify manually rolled exponentiation, maybe this could be a function on ring?
            for _ in 0..exp {
                coef = self.coef_domain.mul(coef, alpha.clone());
            }
            new = self.add(new, self.term_to_element(Term { mono, coef }));
        }
        new
    }

    // Modern Computer Algebra page 599
    // Simplified for case where s = 1 (TODO: Add general case?)
    // TODO: Make a special return type that holds quotient and remainder?
    pub fn multi_poly_div(
        &self,
        dividend: Poly<C::Element>,
        divisor: Poly<C::Element>,
    ) -> (Poly<C::Element>, Poly<C::Element>)
    where
        C: EuclideanDomain, // We need to be able to check divisibility of coefficients
        C::Element: Clone,
    {
        let mut r = self.zero();
        let mut p = dividend;
        let mut q = self.zero();

        let lt_f = divisor.leading_term().expect("not zero").to_owned();

        while !p.is_zero() {
            let lt_p = p.leading_term().expect("not zero").to_owned();
            // TODO: Should divide take a reference to the divisor?
            if let Some(lt_q) = lt_p.clone().try_divide(lt_f.clone(), &self.coef_domain) {
                q = self.add(q, self.term_to_element(lt_q.clone()));
                p = self.sub(
                    p,
                    self.mul(self.term_to_element(lt_q.clone()), divisor.clone()),
                );
            } else {
                r = self.add(r, self.term_to_element(lt_p.clone()));
                p = self.sub(p, self.term_to_element(lt_p.clone()));
            }
        }

        (q, r)
    }

    pub fn divides(&self, dividend: Poly<C::Element>, divisor: Poly<C::Element>) -> bool
    where
        C: EuclideanDomain,
        C::Element: Clone,
    {
        // Check that remainder is zero
        self.multi_poly_div(dividend, divisor).1.is_zero()
    }
}

pub struct TermRef<'a, C> {
    pub mono: &'a Mono,
    pub coef: &'a C,
}

impl<'a, C: Clone> TermRef<'a, C> {
    #[allow(clippy::wrong_self_convention)] // TODO
    pub fn to_owned(self) -> Term<C> {
        Term {
            mono: self.mono.clone(),
            coef: self.coef.clone(),
        }
    }
}

#[derive(Debug, Clone, Eq, PartialEq, Ord, PartialOrd)]
pub struct Term<C> {
    pub coef: C,
    pub mono: Mono,
}

// TODO: Having to implement these manually is annoying.
// Maybe we do want Term to be generic over an element of the ring?

impl<C> Term<C> {
    pub fn try_divide<R>(self, other: Self, domain: &R) -> Option<Self>
    where
        R: EuclideanDomain<Element = C>,
        C: Clone + PartialEq,
    {
        let (coef_q, coef_r) = domain.div(self.coef.clone(), other.coef.clone());

        if coef_r != domain.zero() {
            // Coefficients don't divide
            return None;
        }

        let mono_q = self.mono.try_divide(other.mono)?;

        Some(Term {
            coef: coef_q,
            mono: mono_q,
        })
    }
}

// TODO: Ultimately we don't want these implemented directly on Poly
mod op_impls {
    use std::ops::{Add, Mul, Sub};

    use crate::{
        domain::{Group, Monoid, Ring},
        integer::Int64,
    };

    use super::{Poly, PolyDomain};

    impl Add for Poly<i64> {
        type Output = Poly<i64>;

        fn add(self, rhs: Self) -> Self::Output {
            PolyDomain::<Int64>::default().op(self, rhs)
        }
    }

    impl Sub for Poly<i64> {
        type Output = Poly<i64>;

        fn sub(self, mut rhs: Self) -> Self::Output {
            PolyDomain::<Int64>::default().invert(&mut rhs);
            PolyDomain::<Int64>::default().op(self, rhs)
        }
    }

    impl Mul for Poly<i64> {
        type Output = Poly<i64>;

        fn mul(self, rhs: Self) -> Self::Output {
            PolyDomain::<Int64>::default().mul(self, rhs)
        }
    }
}

#[test]
fn test_uni_eval() {
    let p = Poly::dense("x", vec![1, 2, 3, -1]);
    let alpha = 7i64;
    let res = alpha.pow(3) + 2 * alpha.pow(2) + 3 * alpha.pow(1) - 1;
    let evaled = PolyDomain {
        coef_domain: crate::integer::Int64,
    }
    .eval(p, &Atom::from("x"), alpha);
    assert!(evaled == Poly::from_coef(res))
}

#[test]
fn test_multi_eval() {
    let c1 = Poly::dense("x", vec![-2, 3, 0]);
    let c2 = Poly::dense("x", vec![1, -1]);
    let y = Poly::dense("y", vec![1, 0]);

    let p = c1.clone() * y.clone() + c2.clone();

    let alpha = 11i64;
    let c1_eval = Poly::from_coef(-2 * alpha.pow(2) + 3 * alpha.pow(1));
    let c2_eval = Poly::from_coef(alpha.pow(1) - 1);
    let p_eval = c1_eval * y + c2_eval;

    let evaled = PolyDomain {
        coef_domain: crate::integer::Int64,
    }
    .eval(p, &Atom::from("x"), alpha);
    assert!(evaled == p_eval)
}

#[test]
fn simple_multi_poly_div() {
    let a = Poly::dense("x", vec![1, 2, 3, 4, 2]) * Poly::dense("y", vec![-2, 4]);
    let b = Poly::sparse(["x", "y"], [(1, [1, 0]), (3, [0, 1])]);
    // let b = "x".poly() + 3.poly() * "y".poly();
    let ab = a.clone() * b.clone();

    let (q, r) = PolyDomain {
        coef_domain: crate::integer::Int64,
    }
    .multi_poly_div(ab, a);

    assert!(r.is_zero());
    assert!(q == b);
}
