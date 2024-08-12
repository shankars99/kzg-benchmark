mod kzg;
mod polynomial;

use kzg::KZG;
use polynomial::FpPoly;

fn main() {
    let poly = FpPoly::new(2, 3);
    println!("{:#?}", poly);

    let kzg = KZG::new(poly);
    println!("{:#?}", kzg);
}
