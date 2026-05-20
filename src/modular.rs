use crate::{
    domain::{
        extended_euclidean_algorithm, Domain, EuclideanDomain, Field, Group, Monoid, NormalForm,
        Ring,
    },
    integer::Int64,
    polynomial::{Mono, Poly, PolyDomain, Term},
};

use super::Atom;

pub type Int64ModP = Modular<Int64>;

impl Int64ModP {
    pub fn new(p: i64) -> Self {
        //TODO: check that p is prime and then protect it
        Self { domain: Int64, p }
    }
}

impl Field for Int64ModP {
    fn mul_inverse(&self, element: Self::Element) -> Self::Element {
        let euclid = extended_euclidean_algorithm(&Int64, self.p, element);
        assert!(euclid.r == 1);

        // TODO: Is this normalization sufficient?
        let mut i = euclid.t;
        if i < 0 {
            i += self.p;
        }
        i
    }
}

// TODO: Move this to it's own file
#[derive(Clone, Debug, PartialEq)]
pub struct Modular<D: Domain> {
    pub domain: D,
    pub p: D::Element,
}

impl<D: Domain> Domain for Modular<D> {
    type Element = D::Element;
}

impl<D: EuclideanDomain> Monoid for Modular<D>
where
    D::Element: Clone,
{
    fn identity(&self) -> Self::Element {
        self.domain.identity()
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        self.domain.op_assign(lhs, rhs);
        *lhs = self.domain.div(lhs.clone(), self.p.clone()).1;
    }
}

impl<D: EuclideanDomain> Group for Modular<D>
where
    D::Element: Clone,
{
    fn invert(&self, element: &mut Self::Element) {
        self.domain.invert(element);
        *element = self.domain.div(element.clone(), self.p.clone()).1;
    }
}

impl<D: EuclideanDomain> Ring for Modular<D>
where
    D::Element: Clone,
{
    fn one(&self) -> Self::Element {
        self.domain.one()
    }

    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        let result = self.domain.mul(lhs, rhs);
        self.domain.div(result, self.p.clone()).1
    }
}

impl<D: EuclideanDomain> EuclideanDomain for Modular<D>
where
    D::Element: Clone,
{
    fn div(
        &self,
        dividend: Self::Element,
        divisor: Self::Element,
    ) -> (Self::Element, Self::Element) {
        self.domain.div(dividend, divisor)
    }

    fn euclidean_function(&self, element: &Self::Element) -> u32 {
        self.domain.euclidean_function(element)
    }
}

impl<D: NormalForm + EuclideanDomain> NormalForm for Modular<D>
where
    D::Element: Clone,
{
    fn normal(&self, element: Self::Element) -> Self::Element {
        self.domain.normal(element)
    }

    fn inverse_leading_unit(&self, element: &Self::Element) -> Self::Element {
        self.domain.inverse_leading_unit(element)
    }
}

// Algorithm 5.2 on page 188 of Algorithms for Computer Algebra
// TODO: Delete if unused
pub fn newton_interp(
    ring: &PolyDomain<Int64ModP>,
    atom: Atom,
    alpha: &[i64],
    u: Vec<Poly<i64>>,
) -> Poly<i64> {
    let x = ring.term_to_element(Term {
        coef: 1,
        mono: Mono::power(atom, 1),
    });

    let n = alpha.len();
    assert!(u.len() == n, "not enough evals");
    assert!(n > 1);

    let phi = &ring.coef_domain;

    // Step 1: Compute gammas
    let mut gamma = vec![0]; // TODO: Placeholder value!

    for k in 1..n {
        let mut product = ring.coef_domain.sub(alpha[k], alpha[0]);

        for i in 1..(k - 1) {
            product = phi.mul(product, phi.sub(alpha[k], alpha[i]));
        }

        gamma.push(phi.mul_inverse(product));
    }

    // Step 2: Compute newton coefficients
    let mut nu = vec![u[0].clone()];

    for k in 1..n {
        let mut temp = nu[k - 1].clone();
        for j in (0..(k - 1)).rev() {
            temp.ring_mul_by(phi, phi.sub(alpha[k], alpha[j]));
            temp = ring.add(temp, nu[j].clone());
        }
        nu.push(ring.sub(u[k].clone(), temp));
        nu[k].ring_mul_by(phi, gamma[k].clone());
    }

    // Step 3: Convert to standard form
    let mut result = nu[n - 1].clone(); // Should not need to clone nu
    for k in (0..(n - 1)).rev() {
        let a = ring.term_to_element(Term {
            coef: alpha[k],
            mono: Mono::one(),
        });
        result = ring.add(ring.mul(result, ring.sub(x.clone(), a)), nu[k].clone());
    }
    result
}

// Algorithm 5.4 of Modern Computer Algebra page 106
pub fn chinese_remainder<C>(
    values: Vec<(Modular<C>, Poly<C::Element>)>,
) -> (Modular<C>, Poly<C::Element>)
where
    C: EuclideanDomain + NormalForm + PartialEq,
    C::Element: Clone,
    Modular<C>: Clone,
{
    let mut c = Poly::zero();

    // Build the output ring
    // TODO: Maybe make this a helper function
    let mut modus = values.iter().map(|(a, _)| a.clone());
    let mut m = modus.next().expect("at least one input value");
    while let Some(modu) = modus.next() {
        assert!(modu.domain == m.domain);
        m.p = m.domain.mul(m.p, modu.p);
    }

    let ring = PolyDomain { coef_domain: m };

    for (modu, mut v) in values.into_iter() {
        let q = ring
            .coef_domain
            .div(ring.coef_domain.p.clone(), modu.p.clone())
            .0;

        let euclid = extended_euclidean_algorithm(&ring.coef_domain, q.clone(), modu.p.clone());
        assert!(euclid.r == ring.coef_domain.one());

        // ci = vi * si rem mi
        v.ring_mul_by(&ring.coef_domain, euclid.s);
        v.ring_modp(&ring.coef_domain, modu.p);

        v.ring_mul_by(&ring.coef_domain, q);
        c = ring.add(c, v);

        c.ring_modp(&ring.coef_domain, ring.coef_domain.p.clone());
    }

    (ring.coef_domain, c)
}

#[cfg(test)]
mod tests {
    use super::{chinese_remainder, newton_interp, Int64ModP};
    use crate::{domain::Field, polynomial::PolyDomain, Atom, Poly};

    #[test]
    fn test_mul_inverse() {
        let m = Int64ModP::new(1113);
        let a = 2;
        let i = m.mul_inverse(a);

        assert!((a * i) % m.p == 1);
    }

    #[test]
    fn simple_chinese_remainder() {
        let a = Poly::dense("x", vec![3, 1]);
        let b = Poly::dense("x", vec![5, 2]);

        let mut result = chinese_remainder(vec![(Int64ModP::new(5), a), (Int64ModP::new(7), b)]).1;
        result.modp(5 * 7);
        assert!(result == Poly::dense("x", vec![33, 16]));
    }

    #[test]
    fn test_newton_interp() {
        let p = 97;
        let image_x_is_0 = newton_interp(
            &PolyDomain {
                coef_domain: Int64ModP::new(p),
            },
            Atom::from("y"),
            &[0, 1],
            vec![
                Poly::from_coef(-21).to_modp(p),
                Poly::from_coef(-30).to_modp(p),
            ],
        );
        let image_x_is_1 = newton_interp(
            &PolyDomain {
                coef_domain: Int64ModP::new(p),
            },
            Atom::from("y"),
            &[0, 1],
            vec![
                Poly::from_coef(20).to_modp(p),
                Poly::from_coef(17).to_modp(p),
            ],
        );
        let image_x_is_2 = newton_interp(
            &PolyDomain {
                coef_domain: Int64ModP::new(p),
            },
            Atom::from("y"),
            &[0, 1],
            vec![
                Poly::from_coef(-36).to_modp(p),
                Poly::from_coef(-31).to_modp(p),
            ],
        );

        dbg!(&image_x_is_2);

        let res = newton_interp(
            &PolyDomain {
                coef_domain: Int64ModP::new(p),
            },
            Atom::from("x"),
            &[0, 1, 2],
            vec![image_x_is_0, image_x_is_1, image_x_is_2],
        );

        let expected = Poly::sparse(
            ["x", "y"],
            [
                (1, [2, 1]),
                (5, [1, 1]),
                (41, [1, 0]),
                (-9, [0, 1]),
                (-21, [0, 0]),
            ],
        );

        // let expected = "x".poly() * "x".poly() * "y".poly()
        //     + 5.poly() * "x".poly() * "y".poly()
        //     + 41.poly() * "x".poly()
        //     - 9.poly() * "y".poly()
        //     - 21.poly();
        assert!(res == expected.to_modp(p));
    }
}
