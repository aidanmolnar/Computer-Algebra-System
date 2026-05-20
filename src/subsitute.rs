use super::{
    polynomial::{Mono, Term},
    Atom, Equation, Poly, RationalFunc,
};

// TODO: Generalize this to allow for other superdomains (ex. Poly)

impl Atom {
    pub fn substitute(self, replace: &Atom, by: &RationalFunc<i64>) -> RationalFunc<i64> {
        if self == *replace {
            by.clone()
        } else {
            RationalFunc::dense(self, vec![1, 0])
        }
    }
}

impl Mono {
    pub fn substitute(mut self, replace: &Atom, by: &RationalFunc<i64>) -> RationalFunc<i64> {
        let pow = self.remove_var(replace);
        let mut result = RationalFunc::from_term(Term {
            coef: 1,
            mono: self,
        });
        // TODO: Probably a more efficient way to do this
        //       For one, this probably triggers a lot of redundant gcd calls
        for _ in 0..pow {
            result = result * by.clone();
        }
        result
    }
}

impl Poly<i64> {
    pub fn substitute(self, replace: &Atom, by: &RationalFunc<i64>) -> RationalFunc<i64> {
        let mut result = RationalFunc::zero();

        for term in self.into_terms() {
            result = result
                // TODO better helper for constructing a rationalfunc from a coefficient
                + RationalFunc::from_term(Term {
                    coef: term.coef,
                    mono: Mono::one(),
                }) * term.mono.substitute(replace, by)
        }

        result
    }
}

impl RationalFunc<i64> {
    pub fn substitute(self, replace: &Atom, by: &RationalFunc<i64>) -> RationalFunc<i64> {
        let (num, den) = self.into_parts();
        let num = num.substitute(replace, by);
        let den = den.substitute(replace, by);
        num / den
    }
}

impl Equation {
    pub fn substitute(mut self, replace: &Atom, by: &RationalFunc<i64>) -> Self {
        self.lhs = self.lhs.substitute(replace, by);
        self.rhs = self.rhs.substitute(replace, by);
        self
    }
}

#[test]
fn test_substitute() {
    let a = RationalFunc::dense("x", vec![1, 2, 3, 4]) * RationalFunc::dense("y", vec![-2, 1]);
    let b = RationalFunc::dense("y", vec![1, 2, 3, 4]) * RationalFunc::dense("y", vec![-2, 1]);
    assert!(a.substitute(&"x".into(), &RationalFunc::dense("y", vec![1, 0])) == b);
}
