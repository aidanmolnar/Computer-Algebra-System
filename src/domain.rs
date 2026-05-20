// From Modern Computer Algebra (Chapter 25)
//  -- Domains (increasing specificity) --
// Ring: Has addition and multiplication
// Integral domain: No nonzero divisors
//     ---> GCDs always exist below here
// UFDs: Every nonzero element can be written as a product of irreducibles, has pseudodivision (page 55 Algorithms for Computer Algebra)
// Euclidean domain: Has degree function and division property
// Field: every nonzero element is a unit (has multiplicative inverse)

// -- Vocab --
// unit: element of integral domain with a multiplicative inverse
// reducible: element of integral domain that can be written as product of non-units (otherwise irreducible)
//   --> units are neither reducible nor irreducible
// associate: element of integral domain where a is associate of b, if a = u*b for some u
// normal form: element of euclidean domain that is representative of all of its associates
// leading unit: a = u*normal(a), where u is leading unit
// lu(0) = 1, normal(0) = 0

pub trait Domain {
    // TODO: Should comparison be a function on the domain?
    type Element: PartialEq;
}

pub trait Monoid: Domain {
    fn identity(&self) -> Self::Element;

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element);

    fn op(&self, mut lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        self.op_assign(&mut lhs, rhs);
        lhs
    }
}

pub trait Group: Monoid {
    fn invert(&self, element: &mut Self::Element);

    fn inverse(&self, mut element: Self::Element) -> Self::Element {
        self.invert(&mut element);
        element
    }
}

// TODO: Assign operations
pub trait Ring: Group {
    fn zero(&self) -> Self::Element {
        self.identity()
    }
    fn add(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        self.op(lhs, rhs)
    }
    fn sub(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        self.op(lhs, self.inverse(rhs))
    }

    fn one(&self) -> Self::Element;
    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element;
}

pub trait EuclideanDomain: Ring {
    // Returns (quotient, remainder)
    fn div(
        &self,
        dividend: Self::Element,
        divisor: Self::Element,
    ) -> (Self::Element, Self::Element);

    fn euclidean_function(&self, element: &Self::Element) -> u32;
}

// TODO: Should this just require that it's a ring?
pub trait NormalForm: Ring {
    // TODO: Maybe make a general implementation
    fn normal(&self, element: Self::Element) -> Self::Element;
    fn inverse_leading_unit(&self, element: &Self::Element) -> Self::Element;
}

pub trait Field: EuclideanDomain {
    fn mul_inverse(&self, element: Self::Element) -> Self::Element;
}

#[derive(Clone, Copy, Debug)]
pub struct ExtendedEuclidResults<D> {
    pub q: D, // Quotient
    pub r: D, // Remainder (gcd)

    // Extended coefficients
    pub s: D,
    pub t: D,
}

// TODO: Change to use normalized version...
pub fn extended_euclidean_algorithm<D>(
    domain: &D,
    f: D::Element,
    g: D::Element,
) -> ExtendedEuclidResults<D::Element>
where
    D: NormalForm + EuclideanDomain,
    D::Element: PartialEq + Clone,
{
    let mut s0 = domain.inverse_leading_unit(&f);
    let mut s1 = domain.zero();
    let mut t0 = domain.zero();
    let mut t1 = domain.inverse_leading_unit(&g);
    let mut q1 = domain.zero();
    let mut r0 = domain.normal(f);
    let mut r1 = domain.normal(g);

    while r1 != domain.zero() {
        let r2;
        (q1, r2) = domain.div(r0, r1.clone());
        let inv_rho = domain.inverse_leading_unit(&r2);

        let r2 = domain.normal(r2);

        let s2 = domain.mul(
            domain.sub(s0, domain.mul(q1.clone(), s1.clone())),
            inv_rho.clone(),
        );
        let t2 = domain.mul(
            domain.sub(t0, domain.mul(q1.clone(), t1.clone())),
            inv_rho, //
        );

        r0 = r1;
        s0 = s1;
        t0 = t1;

        r1 = r2;
        s1 = s2;
        t1 = t2;
    }

    ExtendedEuclidResults {
        q: q1,
        r: r0,
        s: s0,
        t: t0,
    }
}

pub fn euclidean_algorithm<D>(domain: &D, f: D::Element, g: D::Element) -> D::Element
where
    D: NormalForm + EuclideanDomain,
    D::Element: PartialEq + Clone,
{
    extended_euclidean_algorithm(domain, f, g).r
}
