use std::collections::{BTreeMap, HashMap, HashSet};

use super::{
    domain::Domain,
    integer::Int64,
    polynomial::{Mono, Term},
    rational::RationalFuncDomain,
    Atom, Id, Poly, RationalFunc,
};

pub trait DifAlgebra: Domain {
    fn deriv(&self, element: Self::Element, by: &Id) -> Self::Element;
}

pub struct RationalDifDomain {
    pub domain: RationalFuncDomain<Int64>,
    // Describes what atoms are functions of what other atoms
    pub functions: HashMap<Id, HashSet<Id>>,
}

impl RationalDifDomain {
    fn deriv_atom(&self, atom: Atom, by: &Id) -> RationalFunc<i64> {
        // TODO: Would be simplified by representing all atoms as partials
        let (var, mut wrt) = match atom {
            Atom::Id(id) => {
                if id == *by {
                    return RationalFunc::one();
                }

                (id, BTreeMap::new())
            }
            Atom::Partial { var, wrt } => {
                // The derivative of var by itself would be 1, which would then be destroyed by other partials
                if var == *by {
                    return RationalFunc::zero();
                }
                (var, wrt)
            }
        };

        let Some(dependents) = self.functions.get(&var) else {
            // atom dependes on no variables --> atom is a constant
            return RationalFunc::zero();
        };

        // var doesn't depend on 'by'
        // TODO: What if var indirectly depends on 'by'?
        //  --> Would need to do a recursive tree search or we need to forbid
        //      this situation when constructing domain
        if !dependents.contains(&by) {
            return RationalFunc::zero();
        }

        // Update degree of partial derivative wrt to by
        let num = wrt.entry(by.clone()).or_insert(0);
        *num += 1;

        RationalFunc::dense(Atom::Partial { var, wrt }, vec![1, 0])
    }

    fn deriv_mono(&self, mut mono: Mono, by: &Id) -> RationalFunc<i64> {
        let Some(var) = mono.vars().next().cloned() else {
            // mono is a constant
            return RationalFunc::zero();
        };
        let pow = mono.remove_var(&var);

        let f = RationalFunc::from_term(Term {
            coef: 1,
            mono: Mono::power(var.clone(), pow),
        });
        let fp = RationalFunc::from_term(Term {
            coef: pow as i64,
            mono: Mono::power(var.clone(), pow - 1),
        }) * self.deriv_atom(var, by);

        let gp = self.deriv_mono(mono.clone(), by);
        let g = RationalFunc::from_term(Term { coef: 1, mono });

        f * gp + fp * g
    }

    fn deriv_poly(&self, poly: Poly<i64>, by: &Id) -> RationalFunc<i64> {
        let mut result = RationalFunc::zero();

        for term in poly.into_terms() {
            // TODO: This is a cursed way to make a constant
            result =
                result + RationalFunc::dense("", vec![term.coef]) * self.deriv_mono(term.mono, by);
        }

        result
    }
}

impl Domain for RationalDifDomain {
    type Element = RationalFunc<i64>;
}

impl DifAlgebra for RationalDifDomain {
    fn deriv(&self, element: Self::Element, by: &Id) -> Self::Element {
        let (hi, lo) = element.into_parts();

        let dhi = self.deriv_poly(hi.clone(), by);
        let dlo = self.deriv_poly(lo.clone(), by);
        let hi = RationalFunc::from_poly(hi);
        let lo = RationalFunc::from_poly(lo);

        (lo.clone() * dhi - hi * dlo) / (lo.clone() * lo)
    }
}

#[test]
fn simple_test() {
    use crate::polynomial::PolyDomain;

    let a = RationalFunc::dense("x", vec![1, 2, 3, 2]);

    // TODO: This is rediculous
    let domain = RationalDifDomain {
        domain: RationalFuncDomain {
            poly_domain: PolyDomain { coef_domain: Int64 },
        },
        functions: HashMap::new(),
    };

    let da = domain.deriv(a, &Id("x"));

    assert!(da == RationalFunc::dense("x", vec![3, 4, 3]));
}
