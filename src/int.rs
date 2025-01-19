use std::fmt::Display;
use std::ops::{Add, Mul, Sub};

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

    /// Fast modular exponentiation
    pub fn pow(&self, mut exp: u64) -> Option<Self> {
        let mut base = *self;
        let mut result = Self::new(1, self.modulus)?;

        while exp > 0 {
            if exp & 1 == 1 {
                result = (result * base)?;
            }
            base = (base * base)?;
            exp >>= 1;
        }
        Some(result)
    }

    /// Find multiplicative inverse using Extended Euclidean Algorithm
    pub fn inv(&self) -> Option<Self> {
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
            return None;
        }
        if t < 0 {
            t += self.modulus as i64;
        }
        Self::new(t as u64, self.modulus)
    }

    /// Multiplication by the Barrett algorithm.
    /// Barrett's algorithm uses a pre-computed parameter to speed up the operation. Instead of dividing by modulo n, it uses multiplication and shifts.
    pub fn mul_barrett(&self, other: &ModInt) -> Option<ModInt> {
        if self.modulus != other.modulus {
            return None;
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

        Some(ModInt::new(result, self.modulus).unwrap())
    }

    /// Multiplication by a constant using the Shoup algorithm
    pub fn shoup_mul(&self, pre_comp: &ShoupPrecomp) -> Option<ModInt> {
        if self.modulus != pre_comp.modulus {
            return None;
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

impl Sub for ModInt {
    type Output = Option<ModInt>;

    fn sub(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            return None;
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
    type Output = Option<ModInt>;

    fn mul(self, other: ModInt) -> Self::Output {
        if self.modulus != other.modulus {
            return None;
        }

        let result = ((self.value as u128 * other.value as u128) % self.modulus as u128) as u64;
        Some(ModInt::new(result, self.modulus).unwrap())
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
    /// Create new precomputed values for Shoup's algorithm
    pub fn new(constant: u64, modulus: u64) -> Option<Self> {
        if modulus == 0 {
            return None;
        }

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


/// Main NTT (Number Theoretic Transform) structure
pub struct NTT {
    modulus: u64,
    root: ModInt,
    root_inv: ModInt,
    root_pw: u64,
}

impl NTT {
    /// Create new NTT instance
    pub fn new(modulus: u64, root: u64, root_pw: u64) -> Option<Self> {
        let root = ModInt::new(root, modulus)?;
        // Fermat's little theorem for inverse
        let root_inv = root.pow(modulus - 2)?;

        Some(Self {
            modulus,
            root,
            root_inv,
            root_pw,
        })
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
    pub fn transform_cooley_tukey(&self, mut values: Vec<u64>) -> Option<Vec<u64>> {
        let n = values.len();
        if !Self::is_power_of_two(n) || n > (1 << self.root_pw) {
            return None;
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
            let wlen = self.root.pow((self.modulus - 1) / (2 * size as u64))?;

            for i in (0..n).step_by(2 * size) {
                let mut w = ModInt::new(1, self.modulus)?;

                for j in 0..size {
                    let u = ModInt::new(values[i + j], self.modulus)?;
                    let v = (w * ModInt::new(values[i + j + size], self.modulus)?)?;

                    values[i + j] = (u + v)?.value();
                    values[i + j + size] = (u - v)?.value();

                    w = (w * wlen)?;
                }
            }
            size *= 2;
        }

        Some(values)
    }

    /// Inverse NTT using Gentleman-Sande butterfly
    pub fn inverse_transform_gentleman_sande(&self, mut values: Vec<u64>) -> Option<Vec<u64>> {
        let n = values.len();
        if !Self::is_power_of_two(n) {
            return None;
        }

        let mut size = n / 2;
        while size > 0 {
            let wlen = self.root_inv.pow((self.modulus - 1) / (2 * size as u64))?;

            for i in (0..n).step_by(2 * size) {
                let mut w = ModInt::new(1, self.modulus)?;

                for j in 0..size {
                    let u = ModInt::new(values[i + j], self.modulus)?;
                    let v = ModInt::new(values[i + j + size], self.modulus)?;

                    values[i + j] = (u + v)?.value();
                    let temp = ((u - v)? * w)?;
                    values[i + j + size] = temp.value();

                    w = (w * wlen)?;
                }
            }
            size /= 2;
        }

        // Apply scaling factor 1/n
        let n_mod = ModInt::new(n as u64, self.modulus)?;
        let n_inv = n_mod.inv()?;

        for value in values.iter_mut() {
            *value = (ModInt::new(*value, self.modulus)? * n_inv)?.value();
        }

        Some(values)
    }

    /// Multiply polynomials using NTT
    pub fn multiply_polynomials(&self, mut a: Vec<u64>, mut b: Vec<u64>) -> Option<Vec<u64>> {
        let n = a.len().max(b.len());
        let n = n.next_power_of_two() * 2;  // Double size to avoid cyclic convolution
        a.resize(n, 0);
        b.resize(n, 0);

        let a_ntt = self.transform_cooley_tukey(a)?;
        let b_ntt = self.transform_cooley_tukey(b)?;

        let mut c_ntt = vec![0; n];
        for i in 0..n {
            let prod = (
                ModInt::new(a_ntt[i], self.modulus)? *
                    ModInt::new(b_ntt[i], self.modulus)?
            )?;
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
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        let result = (a + b).unwrap();
        println!("({} + {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 + 3) mod 7 = 1

        let a = ModInt::new(3, 5).unwrap();
        let b = ModInt::new(4, 5).unwrap();
        let result2 = (a + b).unwrap();
        println!("({} + {}) = {}", &a, &b, &result2);
        assert_eq!(2, result2.value()); // (3 + 4) mod 5 = 2
    }

    #[test]
    fn test_addition_commutativity() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        let result = (a + b).unwrap();
        println!("({} + {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 + 3) mod 7 = 1

        let result2 = (b + a).unwrap();
        println!("({} + {}) = {}", &b, &a, &result2);
        assert_eq!(1, result2.value()); // (3 + 5) mod 7 = 1

        assert_eq!(result, result2);
    }

    #[test]
    fn test_multiplication() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        let result = (a * b).unwrap();
        println!("({} * {}) = {}", &a, &b, &result);
        assert_eq!(1, (a * b).unwrap().value()); // (5 * 3) mod 7 = 1

        let a = ModInt::new(4, 6).unwrap();
        let b = ModInt::new(5, 6).unwrap();
        let result2 = (a * b).unwrap();
        println!("({} * {}) = {}", &a, &b, &result2);
        assert_eq!(2, result2.value()); // (4 * 5) mod 6 = 2
    }

    #[test]
    fn test_multiplication_commutativity() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        let result = (a * b).unwrap();
        println!("({} * {}) = {}", &a, &b, &result);
        assert_eq!(1, result.value()); // (5 * 3) mod 7 = 1

        let result2 = (b * a).unwrap();
        println!("({} * {}) = {}", &b, &a, &result);
        assert_eq!(1, result2.value()); // (3 * 5) mod 7 = 1

        assert_eq!(result, result2);
    }

    #[test]
    fn test_add_association() {
        let a = ModInt::new(3, 7).unwrap();
        let b = ModInt::new(4, 7).unwrap();
        let c = ModInt::new(5, 7).unwrap();
        // ((3+4)+5) ≡ (3+(4+5)) mod 7 = 5
        let result1 = ((a + b).unwrap() + c).unwrap();
        println!("(({} + {}) + {}) = {}", &a, &b, &c, &result1);
        let result2 = (a + (b + c).unwrap()).unwrap();
        println!("({} + ({} + {}) = {}", &a, &b, &c, &result2);
        assert_eq!(result1, result2);
        assert_eq!(5, result1.value());
        assert_eq!(5, result2.value());
    }

    #[test]
    fn test_mul_association() {
        let a = ModInt::new(2, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        let c = ModInt::new(4, 7).unwrap();
        // ((2 * 3) * 4) ≡ (2 * (3 * 4)) mod 7 = 3
        let result1 = ((a * b).unwrap() * c).unwrap();
        println!("(({} * {}) * {}) = {}", &a, &b, &c, &result1);
        let result2 = (a * (b * c).unwrap()).unwrap();
        println!("({} * ({} * {}) = {}", &a, &b, &c, &result2);
        assert_eq!(result1, result2);
        assert_eq!(3, result1.value());
        assert_eq!(3, result2.value());
    }

    #[test]
    fn test_distributivity() {
        let a = ModInt::new(3, 5).unwrap();
        let b = ModInt::new(2, 5).unwrap();
        let c = ModInt::new(4, 5).unwrap();
        // a⋅(b+c) ≡ (a⋅b)+(a⋅c) mod N
        let result_left = (a * (b + c).unwrap()).unwrap();
        let result_right = ((a * b).unwrap() + (a * c).unwrap()).unwrap();
        assert_eq!(result_left, result_right);
        assert_eq!(3, result_left.value());
        assert_eq!(3, result_right.value());
    }

    #[test]
    fn test_barrett_multiplication() {
        let a = ModInt::new(5, 7).unwrap();
        let b = ModInt::new(3, 7).unwrap();
        assert_eq!(1, a.mul_barrett(&b).unwrap().value()); // (5 * 3) mod 7 = 1
    }

    #[test]
    fn test_barrett_multiplication_2() {
        let a = ModInt::new(17, 23).unwrap();
        let b = ModInt::new(19, 23).unwrap();
        assert_eq!(1, a.mul_barrett(&b).unwrap().value()); // (17 * 19) mod 23 = 1
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
        assert_eq!(a.mul_barrett(&b), None);

        let precomp = ShoupPrecomp::new(3, 11).unwrap();
        assert_eq!(a.shoup_mul(&precomp), None);
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
