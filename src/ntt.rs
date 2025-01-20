use crate::int::ModInt;

/// Main NTT (Number Theoretic Transform) structure
pub struct NTT {
    pub(crate) modulus: u64,
    pub(crate) root: ModInt,
    pub(crate) root_inv: ModInt,
    pub(crate) root_pw: u64,
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

    // #[test]
    fn test_ntt_butterfly_patterns() { // TODO - check
        // Using prime modulus 7681 = 1 + 15*2^9
        let modulus = 7681;
        let root = 17;
        let root_pw = 9;

        let ntt = NTT::new(modulus, root, root_pw);

        // Test data
        let values = vec![1, 2, 3, 4, 0, 0, 0, 0];

        let transformed = ntt.transform_cooley_tukey(values.clone());
        let restored = ntt.inverse_transform_gentleman_sande(transformed);

        for (a, b) in values.iter().zip(restored.iter()) {
            assert_eq!(*a, *b);
        }
    }

    // #[test]
    fn test_polynomial_multiplication() { // TODO - check
        let modulus = 7681;
        let root = 17;
        let root_pw = 9;

        let ntt = NTT::new(modulus, root, root_pw);

        let poly1 = vec![1, 2];  // 1 + 2x
        let poly2 = vec![3, 4];  // 3 + 4x

        let result = ntt.multiply_polynomials(poly1, poly2);

        assert_eq!(result[0], 3);  // Constant term
        assert_eq!(result[1], 10 % modulus);  // x term
        assert_eq!(result[2], 8);  // x^2 term
    }
}
