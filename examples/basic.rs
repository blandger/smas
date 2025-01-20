use smas::int::ModInt;

fn main() {
    let a = ModInt::new(5, 7);
    let b = ModInt::new(3, 7);

    // Regular modular addition
    let sum = a + b;
    println!("sum = {sum}");
    assert_eq!(1, sum.value());
    assert_eq!(7, sum.modulus());

    // Regular modular multiplication
    let product = a * b;
    println!("product = {product}");
    assert_eq!(1, product.value());

    // Multiplication by the Barrett algorithm
    let barrett_product = a.mul_barrett(&b);
    println!("barrett_product = {barrett_product}");
    assert_eq!(1, barrett_product.value());
}
