use smas::int::ModInt;

fn main() {
    let a = ModInt::new(5, 7).unwrap();
    let b = ModInt::new(3, 7).unwrap();

    // Regular modular addition
    let sum = (a + b).unwrap();
    println!("sum = {sum}");
    assert_eq!(sum.value(), 1);
    assert_eq!(sum.modulus(), 7);

    // Regular modular multiplication
    let product = (a * b).unwrap();
    println!("product = {product}");
    assert_eq!(product.value(), 1);

    // Multiplication by the Barrett algorithm
    let barrett_product = a.barrett_mul(&b).unwrap();
    println!("barrett_product = {barrett_product}");
    assert_eq!(barrett_product.value(), 1);
}
