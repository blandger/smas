use std::fmt::Display;
use std::ops::{Add, Mul};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModInt {
    value: u64,
    modulus: u64,
}

impl Display for ModInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} (mod {})", self.value, self.modulus)
    }
}

/// Structure for precomputations using the Shoup algorithm
#[derive(Debug)]
pub struct ShoupPrecomp {
    constant: u64,
    constant_prime: u64,
    modulus: u64,
}

impl ShoupPrecomp {
    /// Performs precomputations for multiplication using the Shoup algorithm
    pub fn new(constant: u64, modulus: u64) -> Option<Self> {
        if modulus == 0 {
            return None;
        }

        // Compute constant' = ⌊constant × 2^64/modulus⌋
        let constant_prime = ((constant as u128) << 64) / (modulus as u128);

        Some(Self {
            constant,
            constant_prime: constant_prime as u64,
            modulus,
        })
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

impl ModInt {
    /// Creates a new number with the given value and modulus
    pub fn new(value: u64, modulus: u64) -> Option<Self> {
        if modulus == 0 {
            None
        } else {
            Some(Self {
                value: value % modulus,
                modulus,
            })
        }
    }

    /// Get value
    pub fn value(&self) -> u64 {
        self.value
    }

    /// Get modulus
    pub fn modulus(&self) -> u64 {
        self.modulus
    }

    /// Multiplication by the Barrett algorithm
    pub fn barrett_mul(&self, other: &ModInt) -> Option<ModInt> {
        if self.modulus != other.modulus {
            return None;
        }

        // Precomputation for Barrett algorithm
        let k = 64; // Для u64
        let r = (1u128 << k) / self.modulus as u128;

        let x = self.value as u128;
        let y = other.value as u128;

        // Multiplying numbers
        let m = x * y;

        // Estimating the quotient
        let t = (((m as u128) * (r as u128)) >> k) as u64;

        // Calculating the remainder
        let u = m as u64 - t * self.modulus;

        // Final correction
        let result = if u >= self.modulus {
            u - self.modulus
        } else {
            u
        };

        Some(ModInt::new(result, self.modulus).unwrap())
    }

    /// Multiplication by a constant using the Shoup algorithm
    pub fn shoup_mul(&self, pre_comp: &ShoupPrecomp) -> Option<ModInt> {
        if self.modulus != pre_comp.modulus {
            return None;
        }

        let result = pre_comp.multiply(self.value);
        Some(ModInt::new(result, self.modulus).unwrap())
    }
}

// Implementation of addition
impl Add for ModInt {
    type Output = Option<ModInt>;

    fn add(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            return None;
        }

        let sum = self.value + other.value;
        let result = if sum >= self.modulus {
            sum - self.modulus
        } else {
            sum
        };

        Some(ModInt::new(result, self.modulus).unwrap())
    }
}

// Implementation of multiplication
impl Mul for ModInt {
    type Output = Option<ModInt>;

    fn mul(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            return None;
        }

        let result = ((self.value as u128 * other.value as u128) % self.modulus as u128) as u64;
        Some(ModInt::new(result, self.modulus).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        assert_eq!((a + b).unwrap().value(), 1); // (5 + 3) mod 7 = 1
    }

    #[test]
    fn test_multiplication() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        assert_eq!((a * b).unwrap().value(), 1); // (5 * 3) mod 7 = 1
    }

    #[test]
    fn test_barrett_multiplication() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        assert_eq!(a.barrett_mul(&b).unwrap().value(), 1); // (5 * 3) mod 7 = 1
    }

    #[test]
    fn test_shoup_multiplication() {
        let modulus = 7;
        let constant = 3;
        let pre_comp = ShoupPrecomp::new(constant, modulus).unwrap();

        let a = ModInt::new(5, modulus).unwrap();
        assert_eq!(a.shoup_mul(&pre_comp).unwrap().value(), 1); // (5 * 3) mod 7 = 1

        // Check several different numbers with the same precomputed constant
        let b = ModInt::new(4, modulus).unwrap();
        assert_eq!(b.shoup_mul(&pre_comp).unwrap().value(), 5); // (4 * 3) mod 7 = 5

        let c = ModInt::new(6, modulus).unwrap();
        assert_eq!(c.shoup_mul(&pre_comp).unwrap().value(), 4); // (6 * 3) mod 7 = 4
    }

    #[test]
    fn test_different_moduli() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 11).unwrap();
        assert_eq!(a + b, None);
        assert_eq!(a * b, None);
        assert_eq!(a.barrett_mul(&b), None);

        let precomp = ShoupPrecomp::new(3, 11).unwrap();
        assert_eq!(a.shoup_mul(&precomp), None);
    }
}
