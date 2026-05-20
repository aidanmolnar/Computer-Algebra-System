use std::collections::BTreeMap;

use crate::domain::{Group, Monoid, Ring};

// Related:
// https://en.wikipedia.org/wiki/Abelian_group#Finitely_generated_abelian_groups

pub fn invert<T, G: Group>(power_group: &G, element: &mut BTreeMap<T, G::Element>) {
    for (_, existing) in element.iter_mut() {
        power_group.invert(existing);
    }
}

pub fn op_assign<T: Ord, G: Monoid>(
    power_group: &G,
    lhs: &mut BTreeMap<T, G::Element>,
    rhs: BTreeMap<T, G::Element>,
) {
    for (element, power) in rhs {
        match lhs.get_mut(&element) {
            Some(x) => {
                power_group.op_assign(x, power);

                // Do not keep elements with zero copies
                if *x == power_group.identity() {
                    lhs.remove(&element);
                }
            }
            None => {
                lhs.insert(element, power);
            }
        }
    }
}

pub fn mul<R: Ring, G: Monoid>(
    coef_ring: &R,
    mono_group: &G,
    lhs: BTreeMap<G::Element, R::Element>,
    rhs: BTreeMap<G::Element, R::Element>,
) -> BTreeMap<G::Element, R::Element>
where
    G::Element: Ord + Clone,
    R::Element: Clone,
{
    let mut result = BTreeMap::new();

    for (lhs_mono, lhs_coef) in lhs.into_iter() {
        for (rhs_mono, rhs_coef) in rhs.clone().into_iter() {
            // Accumulate products into output
            // TODO: sad that this allocates in loop
            let term = BTreeMap::from([(
                mono_group.op(lhs_mono.clone(), rhs_mono),
                coef_ring.mul(lhs_coef.clone(), rhs_coef),
            )]);
            op_assign(coef_ring, &mut result, term);
        }
    }

    result
}
