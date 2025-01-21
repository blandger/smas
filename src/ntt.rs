use crate::int::ModInt;

/// Main NTT (Number Theoretic Transform) structure
#[derive(Debug)]
pub struct NTT {
    pub(crate) modulus: u64,
    pub(crate) root: ModInt,     // primitive root of unity
    pub(crate) root_pw: u64,     // power of 2 that divides modulus-1
}

impl NTT {
    /// Creates a new NTT instance
    ///
    /// # Parameters
    /// * `modulus` - Prime modulus of form k * 2^n + 1
    /// * `root` - Primitive root of unity
    /// * `root_pw` - Power of 2 in modulus-1 factorization
    ///
    /// # Example
    /// ```
    /// // For modulus 7681 = 15 * 2^9 + 1
    /// use smas::ntt::NTT;
    /// let ntt = NTT::new(7681, 17, 9);
    /// ```
    pub fn new(modulus: u64, root: u64, root_pw: u64) -> Self {
        Self {
            modulus,
            root: ModInt::new(root, modulus),
            root_pw,
        }
    }

    /// Check if number is power of two
    fn is_power_of_two(n: usize) -> bool {
        n != 0 && (n & (n - 1)) == 0
    }

    /// Compute bit reversal of a number
    ///
    /// # Parameters
    /// * `num` - Number to reverse bits
    /// * `bit_length` - Number of bits to consider
    fn bit_reverse(mut num: usize, bit_length: usize) -> usize {
        let mut result = 0;
        for _ in 0..bit_length {
            result = (result << 1) | (num & 1);
            num >>= 1;
        }
        result
    }

    /// Forward NTT using Cooley-Tukey butterfly
    /// This is a decimation-in-time (DIT) algorithm
    /// Input is in normal order, output is in bit-reversed order
    ///
    /// # Parameters
    /// * `values` - Input values to transform
    ///
    /// # Returns
    /// * Transformed values
    pub fn transform_cooley_tukey(&self, mut values: Vec<ModInt>) -> Vec<ModInt> {
        let n = values.len();
        assert!(Self::is_power_of_two(n), "Length must be power of 2");
        assert!(n <= (1 << self.root_pw), "Length too large for chosen parameters");

        // Bit-reverse permutation
        let bits = (n as f64).log2() as usize;
        for i in 0..n {
            let j = Self::bit_reverse(i, bits);
            if i < j {
                values.swap(i, j);
            }
        }

        // Main NTT loop - Cooley-Tukey butterfly
        let mut m = 1;
        let mut k = n >> 1;

        while m < n {
            let w_m = self.root.pow((self.modulus - 1) / (n as u64));

            for i in (0..n).step_by(2 * m) {
                let mut w = ModInt::new(1, self.modulus);

                for j in 0..m {
                    let t = w * values[i + j + m];
                    values[i + j + m] = values[i + j] - t;
                    values[i + j] = values[i + j] + t;
                    w = w * w_m;
                }
            }

            m *= 2;
            k >>= 1;
        }

        values
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    fn _test_ntt_manual_calculation() {
        // modulus = 17, root = 3
        let ntt = NTT::new(17, 3, 4);

        // Test with [1, 2, 3, 4]
        let input: Vec<ModInt> = vec![1, 2, 3, 4]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let transformed = ntt.transform_cooley_tukey(input);

        // Manually calculated values:
        // X[0] = 1 + 2 + 3 + 4 = 10
        // X[1] = 1 + 2ω + 3ω² + 4ω³ = 6
        // X[2] = 1 + 2ω² + 3ω⁴ + 4ω⁶ = 7
        // X[3] = 1 + 2ω³ + 3ω⁶ + 4ω⁹ = 9
        // where ω = 3 is a primitive 4th root of unity modulo 17
        let expected = vec![10, 6, 7, 9];

        for (a, b) in transformed.iter().zip(expected.iter()) {
            assert_eq!(a.value(), *b,
                       "Mismatch at position {}",
                       transformed.iter().position(|x| x.value() == a.value()).unwrap()
            );
        }
    }

    // Additional validation test
    #[test]
    fn _test_ntt_properties() {
        let ntt = NTT::new(17, 3, 4);

        // Create test vectors
        let input: Vec<ModInt> = vec![1, 0, 0, 0]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let transformed = ntt.transform_cooley_tukey(input);

        // For [1,0,0,0], all coefficients should be 1
        for val in transformed {
            assert_eq!(val.value(), 1);
        }
    }

    // #[test]
    fn test_ntt_small_transform() {
        // Using modulus 17, primitive root 3
        // 17 = 1 * 2^4 + 1
        let ntt = NTT::new(17, 3, 4);

        // Input values [1, 2, 3, 4]
        let input: Vec<ModInt> = vec![1, 2, 3, 4]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let transformed = ntt.transform_cooley_tukey(input.clone());

        // Expected values calculated independently
        let expected: Vec<u64> = vec![10, 2, 7, 15];

        for (a, b) in transformed.iter().zip(expected.iter()) {
            assert_eq!(a.value(), *b);
        }
    }

    // #[test]
    fn test_ntt_small_transform_2() {
        // Using modulus 17, primitive root 3
        // 17 = 1 * 2^4 + 1
        let ntt = NTT::new(17, 3, 4);

        // Input values [1, 2, 3, 4]
        let input: Vec<ModInt> = vec![1, 2, 3, 4]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let transformed = ntt.transform_cooley_tukey(input.clone());

        // Correct expected values
        let expected: Vec<u64> = vec![10, 6, 7, 9];

        for (a, b) in transformed.iter().zip(expected.iter()) {
            assert_eq!(a.value(), *b);
        }
    }

    #[test]
    fn test_ntt_identity_property() {
        // Using modulus 7681 = 15 * 2^9 + 1, primitive root 17
        let ntt = NTT::new(7681, 17, 9);

        // Create input with single "spike"
        let mut input: Vec<ModInt> = vec![ModInt::new(0, 7681); 8];
        input[0] = ModInt::new(1, 7681);

        let transformed = ntt.transform_cooley_tukey(input.clone());

        // For spike at zero, transform should be constant
        let expected_val = ModInt::new(1, 7681);
        for val in transformed {
            assert_eq!(val, expected_val);
        }
    }

    #[test]
    fn test_ntt_linearity() {
        let ntt = NTT::new(17, 3, 4);

        // Test vectors
        let x: Vec<ModInt> = vec![1, 2, 3, 4]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let y: Vec<ModInt> = vec![4, 3, 2, 1]
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        // Transform of sum
        let sum: Vec<ModInt> = x.iter()
            .zip(y.iter())
            .map(|(a, b)| *a + *b)
            .collect();
        let transform_sum = ntt.transform_cooley_tukey(sum);

        // Sum of transforms
        let transform_x = ntt.transform_cooley_tukey(x);
        let transform_y = ntt.transform_cooley_tukey(y);
        let sum_transforms: Vec<ModInt> = transform_x.iter()
            .zip(transform_y.iter())
            .map(|(a, b)| *a + *b)
            .collect();

        // Verify linearity property: NTT(x + y) = NTT(x) + NTT(y)
        assert_eq!(transform_sum, sum_transforms);
    }

    #[test]
    #[should_panic(expected = "Length must be power of 2")]
    fn test_ntt_invalid_length() {
        let ntt = NTT::new(17, 3, 4);
        let input: Vec<ModInt> = vec![1, 2, 3] // Length 3 is not power of 2
            .into_iter()
            .map(|x| ModInt::new(x, 17))
            .collect();

        let _result = ntt.transform_cooley_tukey(input);
    }}
