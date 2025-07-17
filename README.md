# math
## **About This Library**
This library provides mathematical functions designed for large-number computations, such as those involving `BigUint` and other arbitrary-precision integer types.
All functions are implemented using Rust generics to support both built-in types (`u64`, `usize`, etc.) and big number types (`BigUint`, `BigInt`, etc.).
Due to the nature of large-number arithmetic, most generic types are required to implement `Clone` to ensure correctness without sacrificing flexibility. In particular:
- Many operations internally need to preserve input values
- Avoiding excessive ownership moves improves composability
- `Clone` is used in place of `Copy`, as most big number types do not implement `Copy`
This trade-off ensures that the library is both generic and practical for high-precision computations.