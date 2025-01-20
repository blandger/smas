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
