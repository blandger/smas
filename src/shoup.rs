/// Structure for precomputations using the Shoup algorithm
#[derive(Debug)]
pub struct ShoupPrecomp {
    pub(crate) constant: u64,
    pub(crate) constant_prime: u64,
    pub(crate) modulus: u64,
}

impl ShoupPrecomp {
    /// Create new precomputed values for Shoup's algorithm
    pub fn new(constant: u64, modulus: u64) -> Self {
        if modulus == 0 {
            assert!(modulus > 0, "Modulus must be greater than 0");
        }

        let constant_prime = ((constant as u128) << 64) / (modulus as u128);

        Self {
            constant,
            constant_prime: constant_prime as u64,
            modulus,
        }
    }

    /// Multiplies the number by the precomputed constant modulo
    pub fn multiply(&self, x: u64) -> u64 {
        // Compute q = ⌊constant × x/2^64⌋
        let q = ((self.constant_prime as u128 * x as u128) >> 64) as u64;

        // Compute r = constant × x - q × modulus
        let r = self.constant.wrapping_mul(x).wrapping_sub(q.wrapping_mul(self.modulus));

        // Final correction
        if r >= self.modulus {
            r - self.modulus
        } else {
            r
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::int::ModInt;
    use super::*;

    #[test]
    fn test_shoup_multiplication() {
        let modulus = 7;
        let constant = 3;
        let pre_comp = ShoupPrecomp::new(constant, modulus);

        let a = ModInt::new(5, modulus);
        assert_eq!(a.mul_shoup(&pre_comp).value(), 1); // (5 * 3) mod 7 = 1

        // Check several different numbers with the same precomputed constant
        let b = ModInt::new(4, modulus);
        assert_eq!(b.mul_shoup(&pre_comp).value(), 5); // (4 * 3) mod 7 = 5

        let c = ModInt::new(6, modulus);
        assert_eq!(c.mul_shoup(&pre_comp).value(), 4); // (6 * 3) mod 7 = 4
    }

    #[test]
    #[should_panic]
    fn test_panic_different_modulus_on_mul_shoup() {
        let a = ModInt::new(5, 7);
        let precomp = ShoupPrecomp::new(3, 11);
        let _ = a.mul_shoup(&precomp);
    }
}
