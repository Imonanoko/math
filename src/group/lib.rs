use num_traits::{FromPrimitive, One, Zero};
use rand::Rng;
use num_bigint::{BigUint, RandBigInt};
use std::collections::{HashMap, HashSet};
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

/// Attempts to find a non-trivial factor of a composite number using Pollard's Rho algorithm
/// with customizable parameters.
///
/// This function is **probabilistic** and may fail to find a factor depending on the choice
/// of the random constant `c` and the structure of the input `n`. It is often used as a
/// subroutine in full integer factorization.
///
/// # Parameters
/// - `n`: The composite number to factor.
/// - `c_range`: The inclusive range of values from which to randomly sample the constant `c`
///   used in the polynomial function `f(x) = x² + c (mod n)`.
/// - `retry_count`: Number of different `c` values to try if a factor is not found initially.
/// - `max_iterations`: Maximum iterations to try for each `c` before retrying with a new one.
///
/// # Returns
/// - `Some(f)` where `f` is a non-trivial factor of `n`, if successful.
/// - `None` if no factor was found within the allowed retries and iteration budget.
/// # Notes
/// This function assumes `n > 3`. It does not check whether the result is prime.
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

/// Trait that defines a generic group operation for elements of type `T`.
///
/// This trait represents the binary operation of a group (e.g., addition, multiplication).
/// It is expected to be associative and closed within the set.
///
/// # Example
/// ```
/// let result = op.op(&a, &b); // Compute a ⋆ b
/// ```
pub trait GroupOp<T> {
    fn op(&self, a: &T, b: &T) -> T;
}

/// Group operation for modular multiplication: `(a * b) % modulus`.
///
/// This structure represents the multiplicative group modulo `n`.
/// It requires that all elements used are coprime to `modulus`.
///
/// # Example
/// ℤ₇* = {1, 2, 3, 4, 5, 6} under multiplication modulo 7
#[derive(Clone)]
pub struct MulGroup<T> {
    pub modulus: T,
}
impl<T> GroupOp<T> for MulGroup<T>
where
    T: Clone + Rem<Output = T> + Mul<Output = T>,
{
    fn op(&self, a: &T, b: &T) -> T {
        (a.clone() * b.clone()) % self.modulus.clone()
    }
}

/// Group operation for modular addition: `(a + b) % modulus`.
///
/// This structure represents the additive group modulo `n`.
///
/// # Example
/// ℤ₄ = {0, 1, 2, 3} under addition modulo 4
#[derive(Clone)]
pub struct AddGroup<T> {
    pub modulus: T,
}
impl<T> GroupOp<T> for AddGroup<T>
where
    T: Clone + Rem<Output = T> + Add<Output = T>,
{
    fn op(&self, a: &T, b: &T) -> T {
        (a.clone() + b.clone()) % self.modulus.clone()
    }
}

/// Computes the order of a group element under a specified group operation.
///
/// The order is the smallest positive integer `k` such that:
/// `op.op(g^k, g) == identity`
///
/// # Arguments
///
/// * `element` - The group element whose order is to be calculated.
/// * `identity` - The identity element of the group.
/// * `op` - The group operation.
///
/// # Returns
///
/// The order of the element (i.e., the smallest `k` such that `g^k = e`).
///
/// # Panics
///
/// Panics if the number of iterations exceeds 1,000,000,000,000 (to prevent infinite loops).
pub fn order_of<T, Op>(element: T, identity: T, op: Op) -> u64
where
    Op: GroupOp<T>,
    T: Clone + PartialEq,
{
    let mut x = element.clone();
    let mut count = 1u64;
    while x != identity {
        x = op.op(&x, &element);
        count += 1;
        if count > 1_000_000_000_000 {
            panic!("order_of: iteration limit exceeded — possible infinite group or incorrect identity.");
        }
    }
    count
}

/// Checks whether a given element is a generator of the group.
///
/// # Arguments
///
/// * `g` - The candidate element to test.
/// * `identity` - The identity element of the group.
/// * `op` - The group operation.
/// * `group_elements` - A list of all group elements.
///
/// # Returns
///
/// `true` if `g` is a generator, `false` otherwise.
pub fn is_generator<T, Op>(g: T, identity: T, op: Op, group_elements: &[T]) -> bool
where
    Op: GroupOp<T>,
    T: Clone + Eq + Hash,
{
    if g == identity {
        return false;
    }
    let mut generated = HashSet::new();
    let mut x = identity.clone();
    for _ in 0..group_elements.len() {
        generated.insert(x.clone());
        x = op.op(&x, &g);
    }
    let group_elements: HashSet<_> = group_elements.iter().cloned().collect();
    generated == group_elements
}

/// Finds all generators in the given group.
///
/// A generator is an element that can generate the entire group
/// under repeated application of the group operation.
///
/// # Arguments
///
/// * `identity` - The identity element of the group.
/// * `group_elements` - A list of all group elements.
/// * `op` - The group operation (implements `GroupOp<T>`).
///
/// # Returns
///
/// A list of all generator elements in the group.
pub fn find_generators<T, Op>(identity: T, group_elements: &[T], op: Op) -> Vec<T>
where
    T: Clone + Eq + Hash,
    Op: GroupOp<T> + Clone,
{
    group_elements
        .iter()
        .cloned()
        .filter(|g| is_generator(g.clone(), identity.clone(), op.clone(), group_elements))
        .collect()
}
/// Quickly finds all generators in a cyclic group by checking element order.
///
/// A generator is an element whose order equals the group's order.
/// This method is efficient and only works for cyclic groups.
///
/// # Arguments
///
/// * `identity` - The identity element of the group.
/// * `group_order` - The total number of elements in the group.
/// * `group_elements` - All elements of the group.
/// * `op` - The group operation.
///
/// # Returns
///
/// A list of all elements `g` such that `order_of(g) == group_order`.
pub fn find_generators_by_order<T, Op>(
    identity: T,
    group_order: usize,
    group_elements: &[T],
    op: Op,
) -> Vec<T>
where
    T: Clone + PartialEq,
    Op: GroupOp<T> + Clone,
{
    group_elements
        .iter()
        .cloned()
        .filter(|g| order_of(g.clone(), identity.clone(), op.clone()) == group_order as u64)
        .collect()
}

/// Attempts to find a generator `g` of order `q` in the multiplicative group `Z_p*`,
/// using the method `g = h^((p - 1)/q) mod p`, where `h` is chosen randomly.
///
/// ## Parameters
/// - `p`: A prime number defining the multiplicative group `Z_p*`.
/// - `q`: The desired subgroup order. Must be a prime such that `q | (p - 1)`.
///
/// ## Returns
/// - `Some(g)`: A generator of order `q` in `Z_p*` if found within the attempt limit.
/// - `None`: No valid generator was found after the fixed number of attempts.
pub fn find_generator_of_order_q(p: &BigUint, q: &BigUint) -> Option<BigUint> {
    let one = BigUint::one();
    let two = BigUint::from_u32(2u32).unwrap();
    let p_minus_one = p - &one;
    let exponent = &p_minus_one / q;
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let h = rng.gen_biguint_range(&two, &p_minus_one);
        let g = h.modpow(&exponent, p);
        if g != one {
            return Some(g);
        }
    }
    None
}