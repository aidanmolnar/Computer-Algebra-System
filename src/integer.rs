use crate::domain::{Domain, EuclideanDomain, Group, Monoid, NormalForm, Ring};

#[derive(Default, Clone, Debug, PartialEq)]
pub struct UInt64;

impl Domain for UInt64 {
    type Element = u64;
}

impl Monoid for UInt64 {
    fn identity(&self) -> Self::Element {
        0
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        *lhs += rhs;
    }

    fn op(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        lhs + rhs
    }
}

#[derive(Default, Clone, Debug, PartialEq)]
pub struct Int64;

impl Domain for Int64 {
    type Element = i64;
}

impl Monoid for Int64 {
    fn identity(&self) -> Self::Element {
        0
    }

    fn op_assign(&self, lhs: &mut Self::Element, rhs: Self::Element) {
        *lhs += rhs;
    }

    fn op(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        lhs + rhs
    }
}

impl Group for Int64 {
    fn invert(&self, element: &mut Self::Element) {
        *element = -*element;
    }

    fn inverse(&self, element: Self::Element) -> Self::Element {
        -element
    }
}

impl Ring for Int64 {
    fn one(&self) -> Self::Element {
        1
    }

    fn mul(&self, lhs: Self::Element, rhs: Self::Element) -> Self::Element {
        lhs * rhs
    }
}

impl EuclideanDomain for Int64 {
    fn div(
        &self,
        dividend: Self::Element,
        divisor: Self::Element,
    ) -> (Self::Element, Self::Element) {
        let q = dividend.div_euclid(divisor);
        let r = dividend - q * divisor;
        (q, r)
    }

    fn euclidean_function(&self, element: &Self::Element) -> u32 {
        element.abs() as u32
    }
}

impl NormalForm for Int64 {
    fn normal(&self, element: Self::Element) -> Self::Element {
        element.abs()
    }

    fn inverse_leading_unit(&self, element: &Self::Element) -> Self::Element {
        if *element == 0 {
            // Leading unit of zero is one by definition
            return 1;
        }
        element.signum()
    }
}
