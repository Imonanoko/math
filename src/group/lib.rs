use num_traits::{FromPrimitive, One, Zero};
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{AddAssign, Div, DivAssign, Rem};
pub fn prime_factors_trial_division<T>(mut n: T) -> HashMap<T, usize>
where
    T: Clone
        + Rem<Output = T>
        + Div<Output = T>
        + FromPrimitive
        + PartialEq
        + PartialOrd
        + Zero
        + One
        + DivAssign
        + AddAssign
        + Eq
        + Hash,
{
    let mut factors = HashMap::new();
    let two = T::from_u64(2).unwrap();
    while n.clone() % two.clone() == T::zero() {
        *factors.entry(two.clone()).or_insert(0) += 1;
        n /= two.clone();
    }
    let mut other_factor = T::from_u64(3).unwrap();
    let mut other_factor_squared = other_factor.clone() * other_factor.clone();
    while other_factor_squared <= n {
        while n.clone() % other_factor.clone() == T::zero() {
            *factors.entry(other_factor.clone()).or_insert(0) += 1;
            n /= other_factor.clone();
        }
        other_factor += two.clone();
        other_factor_squared = other_factor.clone() * other_factor.clone();
    }
    if n > T::one() {
        *factors.entry(n).or_insert(0) += 1;
    }
    factors
}
pub fn prime_factors_pollard_rho<T>(){

}
pub fn find_generators() {}
pub fn is_generator() {}
pub fn order_of() {}
pub fn find_generator_with_order() {}
