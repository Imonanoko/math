mod math_utils;
pub mod utils {
    pub use crate::math_utils::{gcd, is_prime_miller_rabin, mod_inv, mod_pow};
}

// use utils::{gcd,mod_pow};
pub mod group;
use group::lib;

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;
    use num_traits::{FromPrimitive,Zero};
    #[test]
    fn test_mod_pow() {
        {
            assert_eq!(utils::mod_pow(2u64, 10u64, 1000u64), 24);
        }
        {
            let base = BigUint::from_u64(5).unwrap();
            let exp = BigUint::from_u64(20).unwrap();
            let modulus = BigUint::from_u64(13).unwrap();

            assert_eq!(
                utils::mod_pow(base, exp, modulus),
                BigUint::from_u64(1).unwrap()
            );
        }
    }
    #[test]
    fn test_mod_inv() {
        assert_eq!(utils::mod_inv(3i64, 11i64), Some(4));
        assert_eq!(utils::mod_inv(10i64, 17i64), Some(12));
        assert_eq!(utils::mod_inv(6i64, 12i64), None);
    }

    #[test]
    fn test_pollards_rho() {
        use num_bigint::BigUint;
        use num_traits::FromPrimitive;

        let n = BigUint::from_u64(8051).unwrap();
        let factor = group::lib::pollards_rho(n.clone());
        assert!(factor.is_some());
        assert!((n.clone() % factor.unwrap()).is_zero());
    }
}
