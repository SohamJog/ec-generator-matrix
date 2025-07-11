use criterion::{criterion_group, criterion_main, Criterion};
use ff::PrimeField;
use rand::thread_rng;

extern crate ec_generator_matrix;

use ec_generator_matrix::field::F127;
use ec_generator_matrix::matrix::Matrix;
use ec_generator_matrix::reed_solomon::{decode_rs_with_matrix, encode_rs_with_matrix, Share};
use rand::prelude::IndexedRandom;

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Once;

fn make_shares_with_erasures(
    generator: &Matrix<F127>,
    message: &[F127],
    t: usize,
) -> Vec<Share<F127>> {
    // let mut shares = encode_rs_with_matrix(message, generator);
    // for i in 0..t {
    //     if i < shares.len() {
    //         shares[i].value = None;
    //     }
    // }
    // shares
    // t random erasures
    let mut shares = encode_rs_with_matrix(message, generator);
    let mut rng = thread_rng();
    let indices: Vec<_> = (0..shares.len()).collect();
    for &i in indices.choose_multiple(&mut rng, t) {
        shares[i].value = None;
    }
    shares
}

fn bench_generator_creation(c: &mut Criterion) {
    let k = 30;
    let n = 50;

    c.bench_function("generator::vandermonde", |b| {
        b.iter(|| Matrix::<F127>::vandermonde(n, k));
    });

    c.bench_function("generator::sequential_systematic", |b| {
        b.iter(|| Matrix::<F127>::systematic_with_sequential_vandermonde(n, k));
    });

    c.bench_function("generator::random_systematic", |b| {
        b.iter(|| Matrix::<F127>::systematic_with_random_vandermonde(n, k));
    });

    c.bench_function("generator::random_vandermonde", |b| {
        b.iter(|| Matrix::<F127>::random_vandermonde(n, k));
    });

    c.bench_function("generator::fully_random", |b| {
        b.iter(|| Matrix::<F127>::random_matrix(n, k));
    });

    c.bench_function("generator::systematic_random_rows", |b| {
        b.iter(|| Matrix::<F127>::systematic_with_random_rows(n, k));
    });
}

fn bench_encoding(c: &mut Criterion) {
    let k = 30;
    let n = 50;
    let message: Vec<_> = (0..k).map(|i| F127::from_u128(i as u128)).collect();

    let gen_vandermonde = Matrix::<F127>::vandermonde(n, k);
    let gen_seq = Matrix::<F127>::systematic_with_sequential_vandermonde(n, k);
    let gen_rand = Matrix::<F127>::systematic_with_random_vandermonde(n, k);

    c.bench_function("encode::vandermonde", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_vandermonde));
    });

    c.bench_function("encode::sequential_systematic", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_seq));
    });

    c.bench_function("encode::random_systematic", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_rand));
    });
    let gen_random_vand = Matrix::<F127>::random_vandermonde(n, k);
    let gen_fully_random = Matrix::<F127>::random_matrix(n, k);
    let gen_systematic_random_rows = Matrix::<F127>::systematic_with_random_rows(n, k);

    c.bench_function("encode::random_vandermonde", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_random_vand));
    });

    c.bench_function("encode::fully_random", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_fully_random));
    });

    c.bench_function("encode::systematic_random_rows", |b| {
        b.iter(|| encode_rs_with_matrix(&message, &gen_systematic_random_rows));
    });
}

fn bench_decoding(c: &mut Criterion) {
    let k = 30;
    let n = 50;
    let message: Vec<_> = (0..k).map(|i| F127::from_u128(i as u128)).collect();
    let erasures = n - k;

    let gen_vandermonde = Matrix::<F127>::vandermonde(n, k);
    let gen_seq = Matrix::<F127>::systematic_with_sequential_vandermonde(n, k);
    let gen_rand = Matrix::<F127>::systematic_with_random_vandermonde(n, k);
    let gen_random_vand = Matrix::<F127>::random_vandermonde(n, k);
    let gen_fully_random = Matrix::<F127>::random_matrix(n, k);
    let gen_systematic_random_rows = Matrix::<F127>::systematic_with_random_rows(n, k);

    static PRINT_ONCE: Once = Once::new();

    fn decode_with_count(
        generator: &Matrix<F127>,
        message: &[F127],
        erasures: usize,
        name: &'static str,
        c: &mut Criterion,
    ) {
        let success = AtomicUsize::new(0);
        let panic_count = AtomicUsize::new(0);

        c.bench_function(name, |b| {
            b.iter(|| {
                let shares = make_shares_with_erasures(generator, message, erasures);
                let result = std::panic::catch_unwind(|| {
                    let _ = decode_rs_with_matrix(&shares, generator, message.len());
                });

                match result {
                    Ok(_) => success.fetch_add(1, Ordering::Relaxed),
                    Err(_) => panic_count.fetch_add(1, Ordering::Relaxed),
                };
            });
        });

        PRINT_ONCE.call_once(|| {
            println!("\n=== Decode Benchmark Panic Summary ===");
        });

        println!(
            "{} → Success: {}, Panics: {}",
            name,
            success.load(Ordering::Relaxed),
            panic_count.load(Ordering::Relaxed)
        );
    }

    // decode_with_count(
    //     &gen_vandermonde,
    //     &message,
    //     erasures,
    //     "decode::vandermonde",
    //     c,
    // );
    // decode_with_count(
    //     &gen_seq,
    //     &message,
    //     erasures,
    //     "decode::sequential_systematic",
    //     c,
    // );
    // decode_with_count(
    //     &gen_rand,
    //     &message,
    //     erasures,
    //     "decode::random_systematic",
    //     c,
    // );
    // decode_with_count(
    //     &gen_random_vand,
    //     &message,
    //     erasures,
    //     "decode::random_vandermonde",
    //     c,
    // );

    // decode_with_count(
    //     &gen_fully_random,
    //     &message,
    //     erasures,
    //     "decode::fully_random",
    //     c,
    // );

    decode_with_count(
        &gen_systematic_random_rows,
        &message,
        erasures,
        "decode::systematic_random_rows",
        c,
    );
}

criterion_group!(
    rs_benches,
    // bench_generator_creation,
    // bench_encoding,
    bench_decoding
);
criterion_main!(rs_benches);
