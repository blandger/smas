use std::fmt::Display;
use std::ops::{Add, Mul, Div, Sub};
use crate::shoup::ShoupPrecomp;

/// ModInt structure for modular arithmetic operations
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModInt {
    value: u64,
    modulus: u64,
}

/// Implementation for Display, when it's printed by println("{}", &ModInt)
impl Display for ModInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} mod {}", self.value, self.modulus)
    }
}

impl ModInt {
    /// Creates a new number with the given value and modulus
    /// # Parameters:
    /// - `value`: The value of the number.
    /// - `modulus`: The modulus.
    ///
    /// # Panic:
    /// Panics if the modulus is 0.
    pub fn new(value: u64, modulus: u64) -> Self {
        assert!(modulus > 0, "Modulus must be greater than 0");
        Self {
            value: value % modulus,
            modulus,
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

    /// Fast modular exponentiation
    pub fn pow(&self, mut exp: u64) -> Self {
        let mut base = *self;
        let mut result = Self::new(1, self.modulus);

        while exp > 0 {
            if exp & 1 == 1 {
                result = result * base;
            }
            base = base * base;
            exp >>= 1;
        }
        result
    }

    /// Find multiplicative inverse using Extended Euclidean Algorithm
    pub fn inv(&self) -> Self {
        let mut t = 0i64;
        let mut newt = 1i64;
        let mut r = self.modulus as i64;
        let mut newr = self.value as i64;

        while newr != 0 {
            let quotient = r / newr;
            (t, newt) = (newt, t - quotient * newt);
            (r, newr) = (newr, r - quotient * newr);
        }

        if r > 1 {
            assert!(r > 1, "r must be greater than 1");
        }
        if t < 0 {
            t += self.modulus as i64;
        }
        Self::new(t as u64, self.modulus)
    }

    /// Multiplication by the Barrett algorithm.
    /// Barrett's algorithm uses a pre-computed parameter to speed up the operation. Instead of dividing by modulo n, it uses multiplication and shifts.
    pub fn mul_barrett_(&self, other: &ModInt) -> ModInt {
        if self.modulus != other.modulus {
            assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        }

        // Precomputation for Barrett algorithm
        let k = std::mem::size_of::<u64>() * 8; // number of bits length for u64
        let r = (1u128 << k) / self.modulus as u128;

        let x = self.value as u128;
        let y = other.value as u128;

        // Multiplying numbers
        let m = x * y;

        // Estimating the quotient
        let t = ((m * r) >> k) as u64;

        // Calculating the remainder
        let u = m as u64 - t * self.modulus;

        // Final correction
        let result = if u >= self.modulus {
            u - self.modulus
        } else {
            u
        };

        ModInt::new(result, self.modulus)
    }

    /// Fast modular reduction (constant-time).
    fn _mod_fast(value: i64, modulus: u64) -> u64 {
        let modulus = modulus as i64;
        let mut r = value % modulus;

        // Make r positive using constant-time operations.
        r += modulus & ((r >> 63) as u64 as i64); // If r < 0, add modulus.
        r as u64
    }

    /// Barrett multiplication (constant-time).
    pub fn mul_barrett(&self, other: &ModInt) -> ModInt {
        assert_eq!(self.modulus, other.modulus, "Modulus must be equal");

        let k = std::mem::size_of::<u64>() * 8; // Number of bits in u64.
        let r = (1u128 << k) / self.modulus as u128; // Precomputed Barrett reduction factor.

        let x = self.value as u128;
        let y = other.value as u128;

        // Perform multiplication and Barrett reduction.
        let m = x * y;
        let t = ((m * r) >> k) as u64; // Estimate quotient.
        let u = m as u64 - t * self.modulus; // Calculate remainder.

        // Correct the remainder using a mask (constant-time).
        let mask = (u >= self.modulus) as u64; // 0 if false, 1 if true.
        let result = u - (self.modulus & mask);

        ModInt::new(result, self.modulus)
    }

    /// Multiplication by a constant using the Shoup algorithm
    pub fn mul_shoup(&self, pre_comp: &ShoupPrecomp) -> ModInt {
        if self.modulus != pre_comp.modulus {
            assert_eq!(self.modulus, pre_comp.modulus, "Modulus must be equal");
        }

        let q = ((pre_comp.constant_prime as u128 * self.value as u128) >> 64) as u64;
        let r = pre_comp.constant.wrapping_mul(self.value)
            .wrapping_sub(q.wrapping_mul(self.modulus));

        let result = if r >= self.modulus {
            r - self.modulus
        } else {
            r
        };

        ModInt::new(result, self.modulus)
    }

    /// Calculate GCD (Greatest Common Divisor) using Euclidean algorithm
    pub fn gcd(mut a: u64, mut b: u64) -> u64 {
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        a
    }

    /// Checks if the number has multiplicative inverse
    pub fn has_inverse(&self) -> bool {
        if self.value == 0 {
            return false;
        }
        Self::gcd(self.value, self.modulus) == 1
    }
}

// Implementation of addition
impl Add for ModInt {
    type Output = ModInt;

    fn add(self, other: ModInt) -> Self::Output {
        assert_eq!(self.modulus, other.modulus, "Modulus must be equal");

        let sum = self.value + other.value;
        let result = if sum >= self.modulus {
            sum - self.modulus
        } else {
            sum
        };

        ModInt::new(result, self.modulus)
    }
}

// Implement standard arithmetic operations for ModInt
impl Sub for ModInt {
    type Output = ModInt;

    fn sub(self, other: ModInt) -> Self::Output {
        assert_eq!(self.modulus, other.modulus, "Modulus must be equal");

        let diff = self.value.wrapping_sub(other.value); // Wrapping subtraction to avoid underflow.
        let mask = (diff >> 63) as u64; // Check if underflow occurred.
        let result = diff.wrapping_add(mask * self.modulus); // Add modulus if underflow occurred.

        ModInt::new(result, self.modulus)
    }
}

// Implement standard arithmetic operations for ModInt
impl Mul for ModInt {
    type Output = ModInt;

    fn mul(self, other: ModInt) -> Self::Output {
        assert_eq!(self.modulus, other.modulus, "Modulus must be equal");

        let result = ((self.value as u128 * other.value as u128) % self.modulus as u128) as u64;
        ModInt::new(result, self.modulus)
    }
}

// Implement standard arithmetic operations for ModInt
impl Div for ModInt {
    type Output = ModInt;

    fn div(self, other: ModInt) -> Self::Output {
        assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        assert!(other.value != 0, "Division by zero");
        assert!(
            other.has_inverse(),
            "Divisor {} does not have multiplicative inverse modulo {}",
            other.value, other.modulus
        );

        // Division is multiplication by multiplicative inverse
        self * other.inv()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_addition() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 7);
        let result = a + b;
        println!("({} + {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 + 3) mod 7 = 1

        let a = ModInt::new(3, 5);
        let b = ModInt::new(4, 5);
        let result2 = a + b;
        println!("({} + {}) = {}", &a, &b, &result2);
        assert_eq!(2, result2.value()); // (3 + 4) mod 5 = 2
    }

    #[test]
    fn test_addition_commutativity() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 7);
        let result = a + b;
        println!("({} + {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 + 3) mod 7 = 1

        let result2 = b + a;
        println!("({} + {}) = {}", &b, &a, &result2);
        assert_eq!(1, result2.value()); // (3 + 5) mod 7 = 1

        assert_eq!(result, result2);
    }

    #[test]
    fn test_multiplication() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 7);
        let result = a * b;
        println!("({} * {}) = {}", &a, &b, &result);
        assert_eq!(1, (a * b).value()); // (5 * 3) mod 7 = 1

        let a = ModInt::new(4, 6);
        let b = ModInt::new(5, 6);
        let result2 = a * b;
        println!("({} * {}) = {}", &a, &b, &result2);
        assert_eq!(2, result2.value()); // (4 * 5) mod 6 = 2
    }

    #[test]
    fn test_multiplication_commutativity() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 7);
        let result = a * b;
        println!("({} * {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 * 3) mod 7 = 1

        let result2 = b * a;
        println!("({} * {}) = {}", &b, &a, &result);
        assert_eq!(1, result2.value()); // (3 * 5) mod 7 = 1

        assert_eq!(result, result2);
    }

    #[test]
    fn test_add_association() {
        let a = ModInt::new(3, 7);
        let b = ModInt::new(4, 7);
        let c = ModInt::new(5, 7);
        // ((3+4)+5) ≡ (3+(4+5)) mod 7 = 5
        let result1 = (a + b) + c;
        println!("(({} + {}) + {}) = {}", &a, &b, &c, &result1);
        let result2 = a + (b + c);
        println!("({} + ({} + {}) = {}", &a, &b, &c, &result2);
        assert_eq!(result1, result2);
        assert_eq!(5, result1.value());
        assert_eq!(5, result2.value());
    }

    #[test]
    fn test_mul_association() {
        let a = ModInt::new(2, 7);
        let b = ModInt::new(3, 7);
        let c = ModInt::new(4, 7);
        // ((2 * 3) * 4) ≡ (2 * (3 * 4)) mod 7 = 3
        let result1 = (a * b) * c;
        println!("(({} * {}) * {}) = {}", &a, &b, &c, &result1);
        let result2 = a * (b * c);
        println!("({} * ({} * {}) = {}", &a, &b, &c, &result2);
        assert_eq!(result1, result2);
        assert_eq!(3, result1.value());
        assert_eq!(3, result2.value());
    }

    #[test]
    fn test_distributivity() {
        let a = ModInt::new(3, 5);
        let b = ModInt::new(2, 5);
        let c = ModInt::new(4, 5);
        // a⋅(b+c) ≡ (a⋅b)+(a⋅c) mod N
        let result_left = a * (b + c);
        let result_right = (a * b) + (a * c);
        assert_eq!(result_left, result_right);
        assert_eq!(3, result_left.value());
        assert_eq!(3, result_right.value());
    }

    #[test]
    fn test_pow_operation() {
        let base1 = ModInt::new(3, 7);
        let exp1 = 4;
        let result1 = base1.pow(exp1);
        assert_eq!(result1.value(), 4); // 3^4 mod 7 = 81 mod 7 = 4

        let base2 = ModInt::new(2, 5);
        let exp2 = 3;
        let result2 = base2.pow(exp2);
        assert_eq!(result2.value(), 3); // 2^3 mod 5 = 8 mod 5 = 3
    }

    #[test]
    fn test_barrett_multiplication() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 7);
        assert_eq!(1, a.mul_barrett(&b).value()); // (5 * 3) mod 7 = 1
    }

    #[test]
    fn test_barrett_multiplication_2() {
        let a = ModInt::new(17, 23);
        let b = ModInt::new(19, 23);
        assert_eq!(1, a.mul_barrett(&b).value()); // (17 * 19) mod 23 = 1
    }

    #[test]
    #[should_panic]
    fn test_panic_different_modulus_on_add() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 11);
        let _ = a + b;
    }

    #[test]
    #[should_panic]
    fn test_panic_different_modulus_on_mul() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 11);
        let _ = a * b;
    }

    #[test]
    #[should_panic]
    fn test_panic_different_modulus_on_mul_barrett() {
        let a = ModInt::new(5, 7);
        let b = ModInt::new(3, 11);
        let _ = a.mul_barrett(&b);
    }

    #[test]
    fn test_division_prime_modulus() {
        // Using prime modulus 17
        let a = ModInt::new(10, 17);
        let b = ModInt::new(3, 17);
        let result = a / b;
        println!("{a} / {b} = {result}");

        // 10/3 ≡ 10 * 6 ≡ 60 ≡ 9 (mod 17), because 3 * 6 ≡ 1 (mod 17)
        assert_eq!(9, result.value());
    }

    #[test]
    fn test_basic_division() {
        // Test case 1: Basic division
        let a = ModInt::new(8, 13);
        let b = ModInt::new(3, 13);
        // 8/3 mod 13 ≡ 8⋅3(^−1) mod 13 ≡ 8⋅9 mod 13 ≡ 72 mod 13 ≡ 7
        let result = a / b; // 8 / 3 mod 13
        println!("{a} / {b} = {result}");
        assert_eq!(7, result.value()); // 8 * 9 mod 13 = 7
    }

    #[test]
    fn test_division_by_one() {
        // Test case 2: Dividing by 1 should return the same number
        let a = ModInt::new(10, 13);
        let b = ModInt::new(1, 13);
        let result = a / b; // 10 / 1 mod 13
        assert_eq!(result.value(), 10);
    }

    #[test]
    fn test_division_zero_result() {
        // Test case 3: Division where result is 0
        let a = ModInt::new(0, 13);
        let b = ModInt::new(5, 13);
        let result = a / b; // 0 / 5 mod 13
        assert_eq!(result.value(), 0);
    }

    #[test]
    #[should_panic]
    fn test_panic_division_by_zero() {
        // Test case 4: Division by a number without an inverse
        let a = ModInt::new(8, 13);
        let b = ModInt::new(0, 13);
        let _ = a / b; // Division by 0 should panic

    }

    #[test]
    #[should_panic(expected = "Divisor 3 does not have multiplicative inverse modulo 15")]
    fn test_division_no_inverse() {
        let a = ModInt::new(10, 15);
        let b = ModInt::new(3, 15);  // 3 has no inverse mod 15
        let _ = a / b;  // Should panic
    }

    #[test]
    fn test_has_inverse() {
        // Prime modulus - all non-zero numbers have inverse
        let a = ModInt::new(10, 17);
        assert!(a.has_inverse());

        // Composite modulus - only numbers coprime with modulus have inverse
        let b = ModInt::new(4, 15);  // gcd(4,15) = 1
        assert!(b.has_inverse());

        let c = ModInt::new(3, 15);  // gcd(3,15) = 3
        assert!(!c.has_inverse());
    }

    #[test]
    fn tes_gcd() {
        assert_eq!(ModInt::gcd(48, 18), 6); // 6 is gcd
        assert_eq!(ModInt::gcd(17, 5), 1); // the numbers are relatively prime
        assert_eq!(ModInt::gcd(15, 5), 5);
    }
}
