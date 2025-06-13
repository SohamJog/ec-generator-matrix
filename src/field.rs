use ff::{Field, PrimeField};

#[derive(PrimeField)]
#[PrimeFieldModulus = "127"]
#[PrimeFieldGenerator = "3"] // 3 is a generator for GF(127)
#[PrimeFieldReprEndianness = "little"]
pub struct F127([u64; 1]);
