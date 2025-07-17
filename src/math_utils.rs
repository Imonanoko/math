use num_traits::{One, Signed, ToPrimitive, Zero};
use rand::Rng;
use std::cmp::{PartialEq, PartialOrd};
use std::ops::{BitAnd, Div, Mul, Rem, ShrAssign, Sub};
/*
greatest common divisor
*/
pub fn gcd<T>(mut a: T, mut b: T) -> T
where
    T: Clone + Zero + Rem<Output = T>,
{
    while !b.is_zero() {
        let temp = b.clone();
        b = a % b;
        a = temp;
    }
    a
}
/*
modular exponentiation find base^exp mod modulus.
use square and multiply algorithm
*/
pub fn mod_pow<T>(mut base: T, mut exp: T, modulus: T) -> T
where
    T: Clone
        + Zero
        + One
        + PartialEq
        + PartialOrd
        + Mul<Output = T>
        + Rem<Output = T>
        + ShrAssign<u32>
        + BitAnd<Output = T>,
{
    if modulus == T::one() || base == T::zero() {
        return T::zero();
    }
    // if exp == 0 return 1
    let mut res = T::one();
    base = base % modulus.clone();
    // example 3^13 = (3^8)^1*(3^4)^1*(3^2)^0*(3^1)^1
    while exp > T::zero() {
        if (exp.clone() & T::one()) != T::zero() {
            res = (res * base.clone()) % modulus.clone();
        }
        base = base.clone() * base.clone() % modulus.clone();
        exp >>= 1;
    }
    res
}
/*
Computes the modular inverse of `a` modulo `n`, i.e., finds x such that:
a * x ≡ 1 (mod n)
Returns `None` if inverse does not exist (i.e., when gcd(a, n) ≠ 1)
*/
pub fn mod_inv<T>(a: T, n: T) -> Option<T>
where
    T: Clone
        + PartialEq
        + PartialOrd
        + Zero
        + One
        + Signed
        + Sub<Output = T>
        + Div<Output = T>
        + Rem<Output = T>
        + Mul<Output = T>,
{
    if gcd(a.clone(), n.clone()) != T::one() {
        return None;
    }
    let (mut t, mut new_t) = (T::zero(), T::one());
    let (mut r, mut new_r) = (n.clone(), a.clone());

    while !new_r.is_zero() {
        let quotient = r.clone() / new_r.clone();

        let temp_t = t.clone();
        t = new_t.clone();
        new_t = temp_t - quotient.clone() * new_t;

        let temp_r = r.clone();
        r = new_r.clone();
        new_r = temp_r - quotient * new_r;
    }

    // If t < 0, convert to positive modulo
    if t.is_negative() {
        t = t + n;
    }

    Some(t)
}
/*
Checks whether a number is prime using the Miller-Rabin primality test.

This algorithm is based on Fermat's little theorem, which states that
if n is prime and a is any number such that 1 < a < n,
then a^(n-1) ≡ 1 (mod n).

Miller-Rabin strengthens this test by checking whether the base a
reveals n as a composite through squaring behavior.

Returns `true` if `n` is probably prime, and `false` if definitely composite.
The probability of a false positive is at most 4^-round.

This is a probabilistic algorithm: increasing `round` improves accuracy.
*/
pub fn is_prime_miller_rabin<T>(n: &T, round: usize) -> bool
where
    T: Clone
        + One
        + Zero
        + PartialOrd
        + Eq
        + Sub<Output = T>
        + ShrAssign<u32>
        + BitAnd<Output = T>
        + Rem<Output = T>
        + Mul<Output = T>
        + From<u64>
        + ToPrimitive,
{
    if *n <= T::from(3u64) {
        return *n == T::from(2u64) || *n == T::from(3u64);
    }
    // Return false if n is even
    if (n.clone() & T::one()) == T::zero() {
        return false;
    }
    // Decompose n - 1 as 2^s * d, with d odd
    let one = T::one();
    let two = T::from(2u64);
    let n_minus_one = n.clone() - one.clone();
    let mut d = n_minus_one.clone();
    let mut s = 0u64;
    while (d.clone() & one.clone()) == T::zero() {
        d >>= 1;
        s += 1;
    }
    let mut rng = rand::rng();
    // Perform the test 'round' number of times
    for _ in 0..round {
        let a = loop {
            let rand_val = rng.gen_range(2..n.to_u64().unwrap_or(u64::MAX));
            let a_candidate = T::from(rand_val);
            if a_candidate < *n {
                break a_candidate;
            }
        };
        // Compute x = a^d mod n
        let mut x = mod_pow(a.clone(), d.clone(), n.clone());
        // If x == 1 or x == n - 1, this round passes
        if x == one || x == n_minus_one {
            continue;
        }
        // Repeatedly square x up to s - 1 times
        let mut continue_outer = false;
        for _ in 0..s - 1 {
            x = mod_pow(x.clone(), two.clone(), n.clone());
            if x == n_minus_one {
                continue_outer = true;
                break;
            }
        }
        // If x never becomes n - 1, n is definitely composite
        if continue_outer {
            continue;
        }
        return false;
    }
    true
}
