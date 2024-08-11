use ff::{Field, PrimeField};
use polynomen::Poly;
use rand::{thread_rng, Rng};

const G: u64 = 7;

#[derive(PrimeField)]
#[PrimeFieldModulus = "52435875175126190479447740508185965837690552500527637822603658699938581184513"]
#[PrimeFieldGenerator = "7"]
#[PrimeFieldReprEndianness = "little"]
pub struct Fp([u64; 4]);

impl Fp {
    pub fn into_f64(&self) -> f64 {
        let x = self.to_repr();
        let x = x.as_ref();

        let base: u64 = 256;

        let x_as_u64: u64 = (0..8)
            .map(|i| (x[i] as u64) * base.pow(i as u32))
            .fold(0, |acc, term| acc + term);

        x_as_u64 as f64
    }
}

#[derive(Debug, Clone)]
pub struct FpPoly {
    pub c: Poly<f64>,
    pub c_Fp: Vec<Fp>,
    pub g: u64,
    pub n: u64,
}

impl FpPoly {
    pub fn new(n: u64, max_c: u64) -> Self {
        let mut rng = thread_rng();

        // Create the randomly generated coefficient u64 vector then
        // iterate through it and create an equivalent Fp array
        let c: Vec<f64> = (0..n).map(|_| rng.gen_range(0..=max_c) as f64).collect();
        let c_fp: Vec<Fp> = (c.iter()).map(|&r| Fp::from(r as u64)).collect();

        Self {
            c: Poly::new_from_roots(&c),
            c_Fp: c_fp,
            g: G,
            n: c.len() as u64,
        }
    }
    pub fn from_f64(c: Vec<f64>) -> Self {
        let c_fp: Vec<Fp> = (c.iter()).map(|&r| Fp::from(r as u64)).collect();

        Self {
            c: Poly::new_from_coeffs(&c),
            c_Fp: c_fp,
            g: G,
            n: c.len() as u64,
        }
    }

    pub fn evaluate(&self, x: Fp) -> Fp {
        let n = self.c.degree().unwrap() as u64;
        // computes SUM( f_i * x ** i )
        let evaluation: Fp = (0..n)
            .map(|i| self.c_Fp[i as usize] * x.pow([i]))
            .fold(Fp::ZERO, |acc, term| acc + term);

        evaluation
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    const N: u64 = 3;
    const MAX_COEFFICIENT: u64 = 2;

    #[test]
    fn create_poly() {
        let poly = FpPoly::new(N, MAX_COEFFICIENT);
        assert_eq!(poly.c.degree().unwrap(), N as usize);
    }

    #[test]
    fn poly_division() {
        let poly_n = FpPoly::new(N, MAX_COEFFICIENT).c;
        let poly_d = Poly::new_from_roots(&[2.]);
        let poly_q = poly_n.clone() / poly_d.clone();

        let x = 5.0;

        let poly_n_eval = poly_n.eval_by_val(x) as u64;
        let poly_d_eval = poly_d.eval_by_val(x) as u64;
        let poly_q_eval = poly_q.eval_by_val(x) as u64;

        assert_eq!(
            poly_n_eval / poly_d_eval,
            poly_q_eval,
            "fractions smh. re-run and it should be good"
        );
    }

    #[test]
    fn poly_subtraction() {
        let poly_n = FpPoly::new(N, MAX_COEFFICIENT).c;
        let poly_d = Poly::new_from_roots(&[2.]);
        let poly_q = poly_n.clone() - poly_d.clone();

        let x = 5.0;

        let poly_n_eval = poly_n.eval_by_val(x);
        let poly_d_eval = poly_d.eval_by_val(x);
        let poly_q_eval = poly_q.eval_by_val(x);

        assert_eq!(poly_n_eval - poly_d_eval, poly_q_eval);
    }

    #[test]
    fn fp_to_u64() {
        let x = u64::MAX as f64;
        let x_fp = Fp::from(x as u64);

        let y = x_fp.into_f64();
        assert_eq!(x, y, "the two are not equal");
    }
}
