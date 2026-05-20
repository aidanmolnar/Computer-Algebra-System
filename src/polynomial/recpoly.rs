use std::collections::BTreeMap;

use crate::{
    domain::{Monoid, Ring},
    polynomial::unipoly::UniPolyDomain,
    Atom,
};

use super::{
    poly::{PolyDomain, Term},
    Poly, UniPoly,
};

// A polynomial where the coefficients are univariate polynomials
pub type RecPoly<C> = Poly<UniPoly<C>>;

impl<C> Poly<C> {
    pub fn move_var_to_coefs<R>(self, ring: &UniPolyDomain<R>, var: &Atom) -> RecPoly<C>
    where
        R: Ring<Element = C>,
        C: Clone + PartialEq,
    {
        // Maps from monomials in [x1, x2, x3] to polys with coefs in Z[var]
        let mut exploded = BTreeMap::new();
        for Term { mut mono, coef } in self.into_terms() {
            let pow = mono.remove_var(var);
            let extracted = UniPoly::from_power(&ring.coef_domain, pow, coef);

            let existing = exploded.entry(mono.clone()).or_insert(ring.zero());
            // TODO: Add add_assign, sub_assign, and mul_assign to Ring. op_assign here is addition
            ring.op_assign(existing, extracted);
        }

        Poly::from_map(ring, exploded)
    }
}

impl<C> RecPoly<C> {
    // Largest power of the univariate polynomial coefficients
    pub fn coef_deg(&self) -> u32 {
        self.iter_terms()
            .map(|term| term.coef.deg())
            .max()
            .unwrap_or(0)
    }
}

impl<C> RecPoly<C> {
    // Converts back to a dense polynomial by evaluating the coefficient polynomials at the supplied alpha
    pub fn coef_eval<R>(self, ring: &UniPolyDomain<R>, alpha: C) -> Poly<C>
    where
        C: Clone + PartialEq,
        // UniPolyDomain<R>: Ring<Element = UniPoly<C>> + Clone,
        R: Ring<Element = C>,
    {
        Poly::from_map(
            &ring.coef_domain,
            self.into_terms()
                .map(|term| (term.mono, term.coef.eval(&ring.coef_domain, alpha.clone())))
                .collect(),
        )
    }

    // Converts back to a dense polynomial by evaluating the coefficient polynomials at the supplied alpha
    pub fn collapse<R>(self, ring: &PolyDomain<R>, var: Atom) -> Poly<C>
    where
        R: Ring<Element = C>,
        C: Clone + PartialEq,
    {
        let mut result = Poly::zero();

        for Term { mono, coef } in self.into_terms() {
            let term = ring.mul(
                ring.term_to_element(Term {
                    coef: ring.coef_domain.one(),
                    mono,
                }),
                coef.to_poly(&ring.coef_domain, var.clone()),
            );
            result = ring.add(result, term);
        }

        result
    }
}
