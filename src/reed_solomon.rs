use crate::field::F127;
use crate::matrix::Matrix;
use ff::PrimeField;

#[derive(Clone, Debug, PartialEq)]
pub struct Share<F: PrimeField> {
    pub index: usize,
    pub value: Option<F>,
}

pub fn encode_rs_with_matrix(message: &[F127], generator: &Matrix<F127>) -> Vec<Share<F127>> {
    let k = message.len();
    assert_eq!(generator.col_count(), k);

    let message_col = Matrix::new_with_data(message.iter().map(|x| vec![*x]).collect());
    let codeword = generator.multiply(&message_col);

    (0..generator.row_count())
        .map(|i| Share {
            index: i,
            value: Some(codeword.get(i, 0)),
        })
        .collect()
}

pub fn decode_rs_with_matrix(
    shares: &[Share<F127>],
    generator: &Matrix<F127>,
    k: usize,
) -> Vec<F127> {
    assert!(shares.len() >= k);

    let known: Vec<&Share<F127>> = shares.iter().filter(|s| s.value.is_some()).collect();
    assert!(known.len() >= k, "Too many erasures");

    let mut rows = vec![];
    let mut vals = vec![];

    for share in known.iter().take(k) {
        let row = generator.get_row(share.index).to_vec();
        rows.push(row);
        vals.push(vec![share.value.unwrap()]);
    }

    let generator_submatrix = Matrix::new_with_data(rows);
    let values = Matrix::new_with_data(vals);

    let inverse = generator_submatrix.invert().expect("Singular matrix");
    let message_matrix = inverse.multiply(&values);

    (0..k).map(|i| message_matrix.get(i, 0)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::F127;
    use crate::matrix::Matrix;

    #[test]
    fn test_rs_encode_decode_with_vandermonde_generator() {
        let message = vec![F127::from_u128(1), F127::from_u128(2), F127::from_u128(3)];
        let n = 5;
        let k = 3;

        let generator = Matrix::<F127>::vandermonde(n, k);

        let mut shares = encode_rs_with_matrix(&message, &generator);
        shares[0].value = None;
        shares[4].value = None;

        let recovered = decode_rs_with_matrix(&shares, &generator, k);
        assert_eq!(recovered, message);
    }

    #[test]
    fn test_rs_encode_decode_with_cauchy_generator() {
        let message = vec![F127::from_u128(9), F127::from_u128(42), F127::from_u128(18)];
        let k = message.len();
        let n = 6;

        // Make sure x and y values are disjoint and don't cause division by zero
        let x: Vec<F127> = (0..n).map(|i| F127::from_u128(i as u128)).collect();
        let y: Vec<F127> = (100..100 + k).map(|i| F127::from_u128(i as u128)).collect();

        let generator = Matrix::<F127>::cauchy(&x, &y);

        let mut shares = encode_rs_with_matrix(&message, &generator);

        // 2 erasures
        shares[1].value = None;
        shares[4].value = None;

        let recovered = decode_rs_with_matrix(&shares, &generator, k);
        assert_eq!(recovered, message);
    }

    // #[test]
    // fn test_rs_encode_decode_with_systematic_identity_generator() {
    //     let message = vec![
    //         F127::from_u128(10),
    //         F127::from_u128(20),
    //         F127::from_u128(30),
    //     ];
    //     let k = message.len();
    //     let n = 6;

    //     let generator = Matrix::<F127>::systematic_identity_generator(n, k);

    //     let mut shares = encode_rs_with_matrix(&message, &generator);

    //     // Simulate erasures
    //     shares[2].value = None;
    //     shares[4].value = None;

    //     let recovered = decode_rs_with_matrix(&shares, &generator, k);
    //     assert_eq!(recovered, message);
    // }

    #[test]
    fn test_rs_with_sequential_vandermonde_systematic() {
        let message = vec![F127::from_u128(1), F127::from_u128(2), F127::from_u128(3)];
        let k = message.len();
        let n = 6;

        let generator = Matrix::<F127>::systematic_with_sequential_vandermonde(n, k);

        let mut shares = encode_rs_with_matrix(&message, &generator);
        shares[1].value = None;
        shares[4].value = None;

        let recovered = decode_rs_with_matrix(&shares, &generator, k);
        assert_eq!(recovered, message);
    }

    #[test]
    fn test_rs_with_random_vandermonde_systematic() {
        let message = vec![
            F127::from_u128(10),
            F127::from_u128(20),
            F127::from_u128(30),
        ];
        let k = message.len();
        let n = 6;

        let generator = Matrix::<F127>::systematic_with_random_vandermonde(n, k);

        let mut shares = encode_rs_with_matrix(&message, &generator);
        shares[0].value = None;
        shares[5].value = None;

        let recovered = decode_rs_with_matrix(&shares, &generator, k);
        assert_eq!(recovered, message);
    }

    #[test]
    #[should_panic]
    fn test_rs_encode_decode_with_test_generator() {
        let message = vec![F127::from_u128(1), F127::from_u128(2), F127::from_u128(3)];
        let n = 5;
        let k = 3;

        let generator = Matrix::<F127>::test_generator(n, k);

        let mut shares = encode_rs_with_matrix(&message, &generator);
        shares[1].value = None;
        shares[3].value = None;

        let recovered = decode_rs_with_matrix(&shares, &generator, k);
        // assert_eq!(recovered, message);
    }
}
