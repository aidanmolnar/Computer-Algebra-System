use crate::{domain::Ring, integer::Int64, polynomial::PolyDomain};

use super::{rational::RationalFuncDomain, Atom, RationalFunc};

// TODO: Could be replaced in future by gaussian elimination or other more sophisticated method
// TODO: Maybe equation should just store poylnomials because we flatten them anyway
#[derive(Debug, Clone, PartialEq)]
pub struct Equation {
    pub lhs: RationalFunc<i64>,
    pub rhs: RationalFunc<i64>,
}

impl Equation {
    // Support functions we need:
    // 1. try_into_poly on rational
    // 2. Contains var on poly

    pub fn implicit(expr: RationalFunc<i64>) -> Self {
        Self {
            lhs: RationalFunc::zero(),
            rhs: expr,
        }
    }

    // TODO: How does this handle poles / zeros in denominator
    pub fn try_linear_solve_for(self, var: &Atom) -> Result<RationalFunc<i64>, ()> {
        // 1. Multiply by the denominators on both sides so that we have a flat equation
        let (lhs_num, lhs_den) = self.lhs.into_parts();
        let (rhs_num, rhs_den) = self.rhs.into_parts();
        let combined = lhs_num * rhs_den - rhs_num * lhs_den;
        // TODO: Could divide by gcd / try to remove content?

        // 2. Move all terms that contain var to the lhs and all those that don't to the rhs
        let domain = PolyDomain { coef_domain: Int64 };
        let mut lhs = domain.zero();
        let mut rhs = domain.zero();

        for mut term in combined.into_terms() {
            let pow = term.mono.remove_var(&var);
            if pow == 0 {
                rhs = rhs - domain.term_to_element(term);
            } else if pow == 1 {
                lhs = lhs + domain.term_to_element(term);
            } else {
                // TODO: Try to handle this case better.
                //  This is like looking at it as a univariate polynomial where
                //  the coefficients are multivariate in other vars.
                return Err(());
            }
        }

        //  b. Return rhs / factored out lhs if it was linear in target
        let rat_domain = RationalFuncDomain {
            poly_domain: domain,
        };
        Ok(rat_domain.from_frac(rhs, lhs))
    }
}

#[test]
fn test_linsolve() {
    use crate::Poly;

    let domain = RationalFuncDomain {
        poly_domain: PolyDomain { coef_domain: Int64 },
    };

    let expr = RationalFunc::dense("x", vec![2, 1]);
    let x = Equation::implicit(expr)
        .try_linear_solve_for(&Atom::from("x"))
        .unwrap();
    assert!(x == domain.from_frac(Poly::dense("", vec![-1]), Poly::dense("", vec![2])))
}
