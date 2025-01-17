use smas::int::{ModInt, ShoupPrecomp};

fn main() {
    // Create precomputed values for a constant
    let modulus = 7;
    let constant = 3;
    let pre_comp = ShoupPrecomp::new(constant, modulus).unwrap();


    // Create a number that we want to multiply by the constant
    let a = ModInt::new(5, modulus).unwrap();
    println!("a = {a}");

    // Perform multiplication
    let result = a.shoup_mul(&pre_comp).unwrap();
    println!("result = {result}");
    assert_eq!(result.value(), 1); // (5 * 3) mod 7 = 1
}
