use num_traits::{FromPrimitive, One, ToPrimitive, Zero};
use rand::Rng;
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, Rem, Sub};
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
/*
Attempts to find a non-trivial factor of a composite number using Pollard's Rho algorithm
with customizable parameters.

This function is **probabilistic** and may fail to find a factor depending on the choice
of the random constant `c` and the structure of the input `n`. It is often used as a
subroutine in full integer factorization.

# Parameters
- `n`: The composite number to factor.
- `c_range`: The inclusive range of values from which to randomly sample the constant `c`
  used in the polynomial function `f(x) = x² + c (mod n)`.
- `retry_count`: Number of different `c` values to try if a factor is not found initially.
- `max_iterations`: Maximum iterations to try for each `c` before retrying with a new one.

# Returns
- `Some(f)` where `f` is a non-trivial factor of `n`, if successful.
- `None` if no factor was found within the allowed retries and iteration budget.
# Notes
This function assumes `n > 3`. It does not check whether the result is prime.
*/
pub fn pollards_rho<T>(
    n: T,
    c_range: std::ops::RangeInclusive<u64>,
    retry_count: usize,
    max_iterations: usize,
) -> Option<T>
where
    T: Clone
        + PartialEq
        + Eq
        + PartialOrd
        + Rem<Output = T>
        + Sub<Output = T>
        + Mul<Output = T>
        + Add<Output = T>
        + FromPrimitive
        + Zero
        + One,
{
    if n <= T::from_u64(3u64).unwrap() {
        return None;
    }
    let one = T::one();
    let mut rng = rand::thread_rng();

    // Retry multiple times with different constants c
    for _ in 0..retry_count {
        let c = T::from_u64(rng.gen_range(c_range.clone())).unwrap();
        let mut x = T::from_u64(2u64).unwrap();
        let mut y = x.clone();
        let f = |x: &T| (x.clone() * x.clone() + c.clone()) % n.clone();

        let mut d = T::one();

        for _ in 0..max_iterations {
            x = f(&x);
            y = f(&f(&y));
            let diff = if x > y {
                x.clone() - y.clone()
            } else {
                y.clone() - x.clone()
            };
            d = crate::utils::gcd(diff, n.clone());
            if d > one && d < n {
                return Some(d);
            }
        }
    }
    None
}
pub fn find_generators() {}
pub fn is_generator() {}
pub fn order_of() {}
pub fn find_generator_with_order() {}
