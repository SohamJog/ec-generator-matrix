use ff::Field;
use ff::PrimeField;
use smallvec::SmallVec;

#[derive(Debug)]
pub enum Error {
    SingularMatrix,
}

macro_rules! acc {
    ($m:ident, $r:expr, $c:expr) => {
        $m.data[$r * $m.col_count + $c]
    };
}

fn flatten<T>(m: Vec<Vec<T>>) -> Vec<T> {
    let mut result: Vec<T> = Vec::with_capacity(m.len() * m[0].len());
    for row in m {
        for v in row {
            result.push(v);
        }
    }
    result
}

#[derive(PartialEq, Debug, Clone)]
pub struct Matrix<F: Field> {
    row_count: usize,
    col_count: usize,
    data: SmallVec<[F; 1024]>,
}

fn exp<F: Field + Copy>(mut base: F, mut power: usize) -> F {
    let mut result = F::ONE;
    while power > 0 {
        if power % 2 == 1 {
            result *= base;
        }
        base *= base;
        power /= 2;
    }
    result
}

impl<F: Field + Copy + PrimeField> Matrix<F> {
    pub fn col_count(&self) -> usize {
        self.col_count
    }

    pub fn row_count(&self) -> usize {
        self.row_count
    }
    pub fn new(rows: usize, cols: usize) -> Self {
        let data = SmallVec::from_vec(vec![F::ZERO; rows * cols]);
        Self {
            row_count: rows,
            col_count: cols,
            data,
        }
    }

    pub fn new_with_data(init_data: Vec<Vec<F>>) -> Self {
        let rows = init_data.len();
        let cols = init_data[0].len();
        for r in &init_data {
            assert_eq!(r.len(), cols, "Inconsistent row sizes");
        }
        let data = SmallVec::from_vec(flatten(init_data));
        Self {
            row_count: rows,
            col_count: cols,
            data,
        }
    }

    pub fn identity(size: usize) -> Self {
        let mut result = Self::new(size, size);
        for i in 0..size {
            acc!(result, i, i) = F::ONE;
        }
        result
    }

    pub fn get(&self, r: usize, c: usize) -> F {
        acc!(self, r, c)
    }

    pub fn set(&mut self, r: usize, c: usize, val: F) {
        acc!(self, r, c) = val;
    }

    pub fn multiply(&self, rhs: &Self) -> Self {
        assert_eq!(self.col_count, rhs.row_count);
        let mut result = Self::new(self.row_count, rhs.col_count);
        for r in 0..self.row_count {
            for c in 0..rhs.col_count {
                let mut val = F::ZERO;
                for i in 0..self.col_count {
                    val += acc!(self, r, i) * acc!(rhs, i, c);
                }
                acc!(result, r, c) = val;
            }
        }
        result
    }

    pub fn augment(&self, rhs: &Self) -> Self {
        assert_eq!(self.row_count, rhs.row_count);
        let mut result = Self::new(self.row_count, self.col_count + rhs.col_count);
        for r in 0..self.row_count {
            for c in 0..self.col_count {
                acc!(result, r, c) = acc!(self, r, c);
            }
            for c in 0..rhs.col_count {
                acc!(result, r, self.col_count + c) = acc!(rhs, r, c);
            }
        }
        result
    }

    pub fn sub_matrix(&self, rmin: usize, cmin: usize, rmax: usize, cmax: usize) -> Self {
        let mut result = Self::new(rmax - rmin, cmax - cmin);
        for r in rmin..rmax {
            for c in cmin..cmax {
                acc!(result, r - rmin, c - cmin) = acc!(self, r, c);
            }
        }
        result
    }

    pub fn get_row(&self, row: usize) -> &[F] {
        let start = row * self.col_count;
        &self.data[start..start + self.col_count]
    }

    pub fn swap_rows(&mut self, r1: usize, r2: usize) {
        if r1 == r2 {
            return;
        }
        let start1 = r1 * self.col_count;
        let start2 = r2 * self.col_count;
        for i in 0..self.col_count {
            self.data.swap(start1 + i, start2 + i);
        }
    }

    pub fn is_square(&self) -> bool {
        self.row_count == self.col_count
    }

    pub fn gaussian_elim(&mut self) -> Result<(), Error> {
        for r in 0..self.row_count {
            if acc!(self, r, r).is_zero().into() {
                for r_below in r + 1..self.row_count {
                    if !bool::from(acc!(self, r_below, r).is_zero()) {
                        self.swap_rows(r, r_below);
                        break;
                    }
                }
            }

            if acc!(self, r, r).is_zero().into() {
                return Err(Error::SingularMatrix);
            }

            let inv = acc!(self, r, r).invert().unwrap();
            for c in 0..self.col_count {
                acc!(self, r, c) *= inv;
            }

            for r_below in r + 1..self.row_count {
                let scale = acc!(self, r_below, r);
                if scale.is_zero().into() {
                    continue;
                }
                for c in 0..self.col_count {
                    // acc!(self, r_below, c) -= scale * acc!(self, r, c);
                    let rhs = acc!(self, r, c);
                    acc!(self, r_below, c) -= scale * rhs;
                }
            }
        }

        for d in 0..self.row_count {
            for r_above in 0..d {
                let scale = acc!(self, r_above, d);
                if scale.is_zero().into() {
                    continue;
                }
                for c in 0..self.col_count {
                    // acc!(self, r_above, c) -= scale * acc!(self, d, c);
                    let rhs = acc!(self, d, c);
                    acc!(self, r_above, c) -= scale * rhs;
                }
            }
        }

        Ok(())
    }

    pub fn invert(&self) -> Result<Self, Error> {
        if !self.is_square() {
            panic!("Trying to invert a non-square matrix");
        }
        let mut work = self.augment(&Self::identity(self.row_count));
        work.gaussian_elim()?;
        Ok(work.sub_matrix(0, self.row_count, self.row_count, self.row_count * 2))
    }

    pub fn vandermonde(rows: usize, cols: usize) -> Self {
        let mut result = Self::new(rows, cols);
        for r in 0..rows {
            let r_a = F::from_u128(r as u128);
            for c in 0..cols {
                acc!(result, r, c) = exp(r_a, c);
            }
        }
        result
    }

    pub fn test_generator(rows: usize, cols: usize) -> Self {
        assert!(rows >= cols);
        let mut m = Matrix::new(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                acc!(m, r, c) = F::from_u128((r * cols + c + 1) as u128);
            }
        }
        m
    }
    pub fn cauchy(x: &[F], y: &[F]) -> Self {
        assert!(x.len() > 0 && y.len() > 0);
        let rows = x.len();
        let cols = y.len();

        let mut result = Self::new(rows, cols);
        for r in 0..rows {
            for c in 0..cols {
                let denom = x[r] + y[c];
                assert!(
                    !bool::from(denom.is_zero()),
                    "Cauchy matrix has zero denominator"
                );
                acc!(result, r, c) = denom.invert().unwrap();
            }
        }
        result
    }
}

#[cfg(test)]
mod tests {
    extern crate alloc;

    use super::Matrix;
    use crate::ff::PrimeField;
    use crate::field::F127;
    use alloc::vec;

    macro_rules! matrix {
        (
            $(
                [ $( $x:expr ),+ ]
            ),*
        ) => (
            Matrix::<F127>::new_with_data(vec![ $( vec![$( F127::from_u128($x) ),*] ),* ])
        );
        ($rows:expr, $cols:expr) => (Matrix::<F127>::new($rows, $cols));
    }

    #[test]
    fn test_matrix_col_count() {
        let m1 = matrix!([1, 0, 0]);
        let m2 = matrix!([0, 0, 0], [0, 0, 0]);
        let m3 = Matrix::<F127>::new(1, 4);

        assert_eq!(3, m1.col_count());
        assert_eq!(3, m2.col_count());
        assert_eq!(4, m3.col_count());
    }

    #[test]
    fn test_matrix_row_count() {
        let m1 = matrix!([1, 0, 0]);
        let m2 = matrix!([0, 0, 0], [0, 0, 0]);
        let m3 = Matrix::<F127>::new(1, 4);

        assert_eq!(1, m1.row_count());
        assert_eq!(2, m2.row_count());
        assert_eq!(1, m3.row_count());
    }

    #[test]
    fn test_matrix_swap_rows() {
        {
            let mut m1 = matrix!([1, 2, 3], [4, 5, 6], [7, 8, 9]);
            let expect = matrix!([7, 8, 9], [4, 5, 6], [1, 2, 3]);
            m1.swap_rows(0, 2);
            assert_eq!(expect, m1);
        }
        {
            let mut m1 = matrix!([1, 2, 3], [4, 5, 6], [7, 8, 9]);
            let expect = m1.clone();
            m1.swap_rows(0, 0);
            assert_eq!(expect, m1);
            m1.swap_rows(1, 1);
            assert_eq!(expect, m1);
            m1.swap_rows(2, 2);
            assert_eq!(expect, m1);
        }
    }

    #[test]
    #[should_panic]
    fn test_inconsistent_row_sizes() {
        matrix!([1, 0, 0], [0, 1], [0, 0, 1]);
    }

    #[test]
    #[should_panic]
    fn test_incompatible_multiply() {
        let m1 = matrix!([0, 1], [0, 1], [0, 1]);
        let m2 = matrix!([0, 1, 2]);

        m1.multiply(&m2);
    }

    #[test]
    #[should_panic]
    fn test_incompatible_augment() {
        let m1 = matrix!([0, 1]);
        let m2 = matrix!([0, 1], [2, 3]);

        m1.augment(&m2);
    }

    #[test]
    fn test_matrix_identity() {
        let m1 = Matrix::<F127>::identity(3);
        let m2 = matrix!([1, 0, 0], [0, 1, 0], [0, 0, 1]);
        assert_eq!(m1, m2);
    }

    // #[test]
    // fn test_matrix_multiply() {
    //     let m1 = matrix!([1, 2], [3, 4]);
    //     let m2 = matrix!([5, 6], [7, 8]);
    //     let actual = m1.multiply(&m2);
    //     let expect = matrix!([11, 22], [19, 42]);
    //     assert_eq!(actual, expect);
    // }
    #[test]
    fn test_matrix_multiply() {
        let m1 = matrix!([1, 2], [3, 4]);
        let m2 = matrix!([5, 6], [7, 8]);
        let actual = m1.multiply(&m2);

        // All values below 127, manually computed over GF(127):
        // [1*5 + 2*7 = 5 + 14 = 19, 1*6 + 2*8 = 6 + 16 = 22]
        // [3*5 + 4*7 = 15 + 28 = 43, 3*6 + 4*8 = 18 + 32 = 50]
        let expect = matrix!([19, 22], [43, 50]);
        assert_eq!(actual, expect);
    }

    #[test]
    fn test_matrix_multiply_wraparound() {
        let m1 = matrix!([200, 128], [255, 254]); // will be reduced mod 127: [73, 1], [1, 0]
        let m2 = matrix!([300, 500], [1000, 50]); // reduced: [46, 119], [111, 50]

        let actual = m1.multiply(&m2);

        // Manually reduce all elements modulo 127:
        // m1 = [200 % 127 = 73, 128 % 127 = 1]
        //      [255 % 127 = 1, 254 % 127 = 0]
        //
        // m2 = [300 % 127 = 46, 500 % 127 = 119]
        //      [1000 % 127 = 111, 50 % 127 = 50]
        //
        // result = [
        //     [73*46 + 1*111 = 3358 + 111 = 3469 % 127 = 36,
        //      73*119 + 1*50 = 8687 + 50 = 8737 % 127 = 121],
        //     [1*46 + 0*111 = 46,
        //      1*119 + 0*50 = 119]
        // ]

        let expected = matrix!([40, 101], [46, 119]);
        assert_eq!(actual, expected);
    }

    // #[test]
    // fn test_matrix_inverse_pass_cases() {
    //     {
    //         let m = matrix!([56, 23, 98], [3, 100, 200], [45, 201, 123])
    //             .invert()
    //             .unwrap();
    //         let expect = matrix!([175, 133, 33], [130, 13, 245], [112, 35, 126]);
    //         assert_eq!(m, expect);
    //     }
    //     {
    //         let m = matrix!(
    //             [1, 0, 0, 0, 0],
    //             [0, 1, 0, 0, 0],
    //             [0, 0, 0, 1, 0],
    //             [0, 0, 0, 0, 1],
    //             [7, 7, 6, 6, 1]
    //         )
    //         .invert()
    //         .unwrap();
    //         let expect = matrix!(
    //             [1, 0, 0, 0, 0],
    //             [0, 1, 0, 0, 0],
    //             [123, 123, 1, 122, 122],
    //             [0, 0, 1, 0, 0],
    //             [0, 0, 0, 1, 0]
    //         );
    //         assert_eq!(m, expect);
    //     }
    // }
    #[test]
    fn test_matrix_inverse_pass_cases() {
        let a = matrix!([1, 1, 1], [1, 2, 4], [1, 3, 9]);

        let a_inv = a.invert().unwrap();
        let id = Matrix::<F127>::identity(3);
        let actual = a.multiply(&a_inv);

        assert_eq!(actual, id);
    }

    #[test]
    #[should_panic]
    fn test_matrix_inverse_non_square() {
        matrix!([56, 23], [3, 100], [45, 201]).invert().unwrap();
    }

    #[test]
    #[should_panic]
    fn test_matrix_inverse_singular() {
        matrix!([4, 2], [12, 6]).invert().unwrap();
    }
}
