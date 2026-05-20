use crate::DifAlgebra;

use super::{Id, RationalDifDomain, RationalFunc};

// Vector type to simplify taking gradient and divergence
// TODO: Generalize this to a tensor type?
pub struct VectorExpr {
    // TODO: Do we want this to be const generic with vars being implicit?
    //       Dependes closely on the RationalDifDomain...
    elements: Vec<RationalFunc<i64>>,
}

impl VectorExpr {
    // Takes the divergence wrt to the provided variables
    pub fn div(self, domain: &RationalDifDomain, vars: &[Id]) -> RationalFunc<i64> {
        assert!(self.elements.len() == vars.len());
        let mut result = RationalFunc::zero();
        for (element, by) in self.elements.into_iter().zip(vars) {
            result = result + domain.deriv(element, by);
        }
        result
    }
}

impl RationalFunc<i64> {
    pub fn grad(self, domain: &RationalDifDomain, vars: &[Id]) -> VectorExpr {
        VectorExpr {
            elements: vars
                .iter()
                .map(|by| domain.deriv(self.clone(), by))
                .collect(),
        }
    }
}
