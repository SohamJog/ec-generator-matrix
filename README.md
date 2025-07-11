### Testing first k erasures with k = 30

This benchmark compares the performance of Reed-Solomon encoding and decoding using different types of generator matrices when recovering from exactly `k` erasures (with `k = 30` and `n = 60`).

#### Matrix Types Explained

* **Vandermonde Matrix**: A classical generator matrix constructed using powers of field elements. This is typically dense and not systematic.
* **Sequential Systematic Matrix**: A generator matrix where the first `k` rows form an identity matrix (making it systematic), and the remaining rows are constructed using a Vandermonde structure with sequential inputs.
* **Random Systematic Matrix**: Similar to the sequential systematic matrix, but the non-identity rows are constructed using random field elements instead of sequential ones. This introduces randomness while preserving systematic properties.

| Benchmark                                | Time (µs) |
| ---------------------------------------- | --------- |
| `generate::vandermonde_matrix`           | 15.109 µs |
| `generate::sequential_systematic_matrix` | 6.3452 µs |
| `generate::random_systematic_matrix`     | 8.1763 µs |
| `encode::using_vandermonde`              | 3.6163 µs |
| `encode::using_sequential_systematic`    | 3.5229 µs |
| `encode::using_random_systematic`        | 3.5089 µs |
| `decode::vandermonde_matrix`             | 224.30 µs |
| `decode::sequential_systematic_matrix`   | 125.69 µs |
| `decode::random_systematic_matrix`       | 123.08 µs |
