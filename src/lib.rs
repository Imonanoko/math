mod math_utils;
pub mod utils {
    pub use crate::math_utils::*;
}

// use utils::{gcd,mod_pow};
mod group;
pub mod groups {
    pub use crate::group::lib::*;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::groups::{find_generators, find_generators_by_order, is_generator};
    use crate::utils::crt;
    use num_bigint::{BigUint, BigInt};
    use num_traits::{FromPrimitive, Zero};
    use std::collections::HashSet;
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
        let factor = group::lib::pollards_rho(n.clone(), 1..=10, 50, 1000);
        assert!(factor.is_some());
        assert!((n.clone() % factor.unwrap()).is_zero());
    }
    #[test]
    fn test_order_of() {
        assert_eq!(
            groups::order_of(3u64, 1u64, groups::MulGroup { modulus: 7 }),
            6u64
        );
        assert_eq!(
            groups::order_of(1u64, 0u64, groups::AddGroup { modulus: 4 }),
            4u64
        );
        assert_eq!(
            groups::order_of(2u64, 0u64, groups::AddGroup { modulus: 4 }),
            2u64
        );
    }
    #[test]
    fn test_is_generator() {
        assert!(is_generator(
            3,
            1,
            groups::MulGroup { modulus: 7 },
            &[1, 2, 3, 4, 5, 6]
        ));
        assert!(is_generator(
            5,
            1,
            groups::MulGroup { modulus: 7 },
            &[1, 2, 3, 4, 5, 6]
        ));
    }
    #[test]
    fn test_find_generators() {
        assert_eq!(
            find_generators(1, &[1, 2, 3, 4, 5, 6], groups::MulGroup { modulus: 7 }),
            vec![3, 5]
        );
    }
    #[test]
    fn test_find_generators_by_order() {
        assert_eq!(
            find_generators_by_order(1, 6, &[1, 2, 3, 4, 5, 6], groups::MulGroup { modulus: 7 }),
            vec![3, 5]
        );
        let p = 719;
        let op = groups::MulGroup { modulus: p };
        let id = 1;
        let group_elements: Vec<_> = (1..p).filter(|x| utils::gcd(*x, p) == 1).collect();
        let group_order = group_elements.len();
        let g1 = find_generators(id, &group_elements, op.clone());
        let g2 = find_generators_by_order(id, group_order, &group_elements, op);
        let set1: HashSet<_> = g1.iter().cloned().collect();
        let set2: HashSet<_> = g2.iter().cloned().collect();
        assert_eq!(set1, set2);
    }
    #[test]
    fn test_crt() {
        let residues = vec![BigInt::from(2), BigInt::from(3), BigInt::from(2)];
        let moduli = vec![BigInt::from(3), BigInt::from(4), BigInt::from(5)];

        let x = crt(&residues, &moduli);
        assert_eq!(x, Some(BigInt::from(47)));
    }
}
