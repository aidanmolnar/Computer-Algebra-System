mod derivative;
pub mod domain;
mod equation;
pub mod integer;
mod modular;
pub mod polynomial;
mod rational;
mod subsitute;
mod vector;

pub use derivative::{DifAlgebra, RationalDifDomain};
pub use equation::Equation;
pub use modular::newton_interp;
pub use polynomial::{Poly, PolyDomain};
pub use rational::{RationalFunc, RationalFuncDomain};

use std::collections::BTreeMap;

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Copy)]
pub struct Id(pub &'static str);

// TODO: Should we get rid of Id field and just use partial where wrt is empty?
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub enum Atom {
    Id(Id),
    // TODO: This is a monoid?
    Partial { var: Id, wrt: BTreeMap<Id, u8> },
}

impl From<&'static str> for Atom {
    fn from(value: &'static str) -> Self {
        Self::Id(Id(value))
    }
}
