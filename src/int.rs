use std::fmt::Display;
use std::ops::{Add, Mul, Sub};
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

    /// Fast modulus operation that handles negative values.
    fn _mod_fast(value: i64, modulus: u64) -> u64 {
        let modulus = modulus as i64;
        let mut r = value % modulus;
        if r < 0 {
            r += modulus;
        }
        r as u64
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
    pub fn mul_barrett(&self, other: &ModInt) -> ModInt {
        if self.modulus != other.modulus {
            assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        }

        // Precomputation for Barrett algorithm
        let k = 64; // bits length for u64
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

    /// Multiplication by a constant using the Shoup algorithm
    pub fn shoup_mul(&self, pre_comp: &ShoupPrecomp) -> ModInt {
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
}

// Implementation of addition
impl Add for ModInt {
    type Output = ModInt;

    fn add(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        }

        let sum = self.value + other.value;
        let result = if sum >= self.modulus {
            sum - self.modulus
        } else {
            sum
        };

        ModInt::new(result, self.modulus)
    }
}

impl Sub for ModInt {
    type Output = ModInt;

    fn sub(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        }

        let result = if self.value >= other.value {
            self.value - other.value
        } else {
            self.modulus - (other.value - self.value)
        };

        ModInt::new(result, self.modulus)
    }
}

// Implement standard arithmetic operations for ModInt
impl Mul for ModInt {
    type Output = ModInt;

    fn mul(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            assert_eq!(self.modulus, other.modulus, "Modulus must be equal");
        }

        let result = ((self.value as u128 * other.value as u128) % self.modulus as u128) as u64;
        ModInt::new(result, self.modulus)
    }
}


/// Main NTT (Number Theoretic Transform) structure
pub struct NTT {
    modulus: u64,
    root: ModInt,
    root_inv: ModInt,
    root_pw: u64,
}

impl NTT {
    /// Create new NTT instance
    pub fn new(modulus: u64, root: u64, root_pw: u64) -> Self {
        let root = ModInt::new(root, modulus);
        // Fermat's little theorem for inverse
        let root_inv = root.pow(modulus - 2);

        Self {
            modulus,
            root,
            root_inv,
            root_pw,
        }
    }

    fn is_power_of_two(n: usize) -> bool {
        n != 0 && (n & (n - 1)) == 0
    }

    fn bit_reverse(mut num: usize, bits: usize) -> usize {
        let mut result = 0;
        for _ in 0..bits {
            result = (result << 1) | (num & 1);
            num >>= 1;
        }
        result
    }

    /// Forward NTT using Cooley-Tukey butterfly
    pub fn transform_cooley_tukey(&self, mut values: Vec<u64>) -> Vec<u64> {
        let n = values.len();
        if !Self::is_power_of_two(n) || n > (1 << self.root_pw) {
            return vec![];
        }

        // Bit-reverse permutation
        let bits = (n as f64).log2() as usize;
        for i in 0..n {
            let j = Self::bit_reverse(i, bits);
            if i < j {
                values.swap(i, j);
            }
        }

        // Main NTT loop - Cooley-Tukey butterfly
        let mut size = 1;
        while size < n {
            let wlen = self.root.pow((self.modulus - 1) / (2 * size as u64));

            for i in (0..n).step_by(2 * size) {
                let mut w = ModInt::new(1, self.modulus);

                for j in 0..size {
                    let u = ModInt::new(values[i + j], self.modulus);
                    let v = w * ModInt::new(values[i + j + size], self.modulus);

                    values[i + j] = (u + v).value();
                    values[i + j + size] = (u - v).value();

                    w = w * wlen;
                }
            }
            size *= 2;
        }

        values
    }

    /// Inverse NTT using Gentleman-Sande butterfly
    pub fn inverse_transform_gentleman_sande(&self, mut values: Vec<u64>) -> Vec<u64> {
        let n = values.len();
        if !Self::is_power_of_two(n) {
            return vec![];
        }

        let mut size = n / 2;
        while size > 0 {
            let wlen = self.root_inv.pow((self.modulus - 1) / (2 * size as u64));

            for i in (0..n).step_by(2 * size) {
                let mut w = ModInt::new(1, self.modulus);

                for j in 0..size {
                    let u = ModInt::new(values[i + j], self.modulus);
                    let v = ModInt::new(values[i + j + size], self.modulus);

                    values[i + j] = (u + v).value();
                    let temp = (u - v) * w;
                    values[i + j + size] = temp.value();

                    w = w * wlen;
                }
            }
            size /= 2;
        }

        // Apply scaling factor 1/n
        let n_mod = ModInt::new(n as u64, self.modulus);
        let n_inv = n_mod.inv();

        for value in values.iter_mut() {
            *value = (ModInt::new(*value, self.modulus) * n_inv).value();
        }

        values
    }

    /// Multiply polynomials using NTT
    pub fn multiply_polynomials(&self, mut a: Vec<u64>, mut b: Vec<u64>) -> Vec<u64> {
        let n = a.len().max(b.len());
        let n = n.next_power_of_two() * 2;  // Double size to avoid cyclic convolution
        a.resize(n, 0);
        b.resize(n, 0);

        let a_ntt = self.transform_cooley_tukey(a);
        let b_ntt = self.transform_cooley_tukey(b);

        let mut c_ntt = vec![0; n];
        for i in 0..n {
            let prod = ModInt::new(a_ntt[i], self.modulus) *
                    ModInt::new(b_ntt[i], self.modulus);
            c_ntt[i] = prod.value();
        }

        self.inverse_transform_gentleman_sande(c_ntt)
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
    fn test_shoup_multiplication() {
        let modulus = 7;
        let constant = 3;
        let pre_comp = ShoupPrecomp::new(constant, modulus);

        let a = ModInt::new(5, modulus);
        assert_eq!(a.shoup_mul(&pre_comp).value(), 1); // (5 * 3) mod 7 = 1

        // Check several different numbers with the same precomputed constant
        let b = ModInt::new(4, modulus);
        assert_eq!(b.shoup_mul(&pre_comp).value(), 5); // (4 * 3) mod 7 = 5

        let c = ModInt::new(6, modulus);
        assert_eq!(c.shoup_mul(&pre_comp).value(), 4); // (6 * 3) mod 7 = 4
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
    #[should_panic]
    fn test_panic_different_modulus_on_mul_shoup() {
        let a = ModInt::new(5, 7);
        let precomp = ShoupPrecomp::new(3, 11);
        let _ = a.shoup_mul(&precomp);
    }

    /*#[test]
    fn test_ntt_butterfly_patterns() { // TODO - check
        // Using prime modulus 7681 = 1 + 15*2^9
        let modulus = 7681;
        let root = 17;
        let root_pw = 9;

        let ntt = NTT::new(modulus, root, root_pw).unwrap();

        // Test data
        let values = vec![1, 2, 3, 4, 0, 0, 0, 0];

        let transformed = ntt.transform_cooley_tukey(values.clone()).unwrap();
        let restored = ntt.inverse_transform_gentleman_sande(transformed).unwrap();

        for (a, b) in values.iter().zip(restored.iter()) {
            assert_eq!(*a, *b);
        }
    }

    #[test]
    fn test_polynomial_multiplication() { // TODO - check
        let modulus = 7681;
        let root = 17;
        let root_pw = 9;

        let ntt = NTT::new(modulus, root, root_pw).unwrap();

        let poly1 = vec![1, 2];  // 1 + 2x
        let poly2 = vec![3, 4];  // 3 + 4x

        let result = ntt.multiply_polynomials(poly1, poly2).unwrap();

        assert_eq!(result[0], 3);  // Constant term
        assert_eq!(result[1], 10 % modulus);  // x term
        assert_eq!(result[2], 8);  // x^2 term
    }*/
}
