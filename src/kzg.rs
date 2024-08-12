use crate::polynomial::{Fp, FpPoly};

use ff::Field;
use polynomen::Poly;
use rand::{thread_rng, Rng};

#[derive(Debug, Clone)]
pub struct KZG {
    pub f: FpPoly,
    pub q: Option<Poly<f64>>,
    pub gp: Vec<Fp>,
    pub C: Fp,
}

impl KZG {
    pub fn new(poly: FpPoly) -> Self {
        let rng = thread_rng();
        let tau = Fp::random(rng);

        /*
         * @note: Don't write the generation of gp functionally as it messes up the readability
         */

        // [g]
        let mut gp = vec![Fp::from(poly.g)];

        for i in 1..poly.n {
            // [g, g*tau, g*tau^2, ...]
            gp.push(*gp.last().unwrap() * tau);
        }

        // Compute the commitment of the polynomial
        let C = Fp::from(poly.g) * poly.evaluate(tau);
        Self {
            f: poly,
            q: None,
            gp,
            C,
        }
    }

    pub fn opening(&self) -> (Self, Fp) {
        // generate u as a fiat shamir random variable
        let u = Self::fiat_shamir_oracle();

        // f(x) - f(u)
        // As f(u) is a point with no degree. We can subtract is from the degree(0) term in f(x)
        let v = self.f.evaluate(u).into_f64();

        let mut c_q = self.f.c.clone();
        c_q[0] = c_q[0] - v;

        let den_q = Poly::new_from_roots(&[u.into_f64()]);
        let q = c_q / den_q;

        (
            Self {
                q: Some(q),
                ..self.clone()
            },
            u,
        )
    }

    fn fiat_shamir_oracle() -> Fp {
        let val = thread_rng().gen_range(0..2_i64.pow(16)) as u64;
        Fp::from(val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const N: u64 = 2;
    const MAX_COEFFICIENT: u64 = 3;

    fn setup() -> KZG {
        let poly = FpPoly::new(N, MAX_COEFFICIENT);
        KZG::new(poly.clone())
    }

    #[test]
    fn new_kzg() {
        let kzg = setup();
        println!("{:#?}", kzg);

        assert_eq!(kzg.gp.len(), N as usize);
    }

    #[test]
    fn open_kzg() {
        let (kzg, u) = setup().opening();
        // println!("{:#?}", kzg);

        let root_term = Poly::new_from_roots(&[u.into_f64()]);
        // q_prime = q * (x-u)
        let q_prime = root_term * kzg.q.unwrap();
        let v = kzg.f.evaluate(u).into_f64();

        // q_prime_and_v = q_prime + v
        let mut q_prime_and_v = q_prime;
        q_prime_and_v[0] = q_prime_and_v[0] + v;

        assert_eq!(kzg.f.c, q_prime_and_v);
    }
}
