use smas::int::ModInt;
use smas::shoup::ShoupPrecomp;

fn main() {
    // Create precomputed values for a constant
    let modulus = 7;
    let constant = 3;
    let pre_comp = ShoupPrecomp::new(constant, modulus);


    // Create a number that we want to multiply by the constant
    let a = ModInt::new(5, modulus);
    println!("a = {a}");

    // Perform multiplication
    let result = a.mul_shoup(&pre_comp);
    println!("result = {result}");
    assert_eq!(result.value(), 1); // (5 * 3) mod 7 = 1
}
