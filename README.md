# zinchi

[![Crates.io](https://img.shields.io/crates/v/zinchi.svg)](https://crates.io/crates/zinchi)
[![Documentation](https://docs.rs/zinchi/badge.svg)](https://docs.rs/zinchi)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![CI](https://github.com/OliverBScott/zinchi/workflows/CI/badge.svg)](https://github.com/OliverBScott/zinchi/actions)

A compact binary representation for InChI Keys.

This crate provides a space-efficient binary encoding for International Chemical Identifier (InChI) keys, reducing their size from the standard 27-byte ASCII representation to either 9 or 14 bytes. The implementation is based on the work by John Mayfield (NextMove Software): [Data Compression of InChI Keys and 2D Coordinates](https://www.nextmovesoftware.com/talks/Mayfield_DataCompressionOfInChIKeysAnd2dCoordinates_NIHINCHI_202103.pdf).

> **Note:** This is a personal project created for fun and to explore Rust. While it implements a real compression algorithm, it's primarily a learning exercise rather than a production-critical library.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
zinchi = "0.1"
```

## Usage

### Parsing and Displaying InChI Keys

```rust
use zinchi::InChIKey;

// Parse an InChI key from a string
let key: InChIKey = "ZZJLMZYUGLJBSO-UHFFFAOYSA-N".parse().expect("Failed to parse InChIKey")

// Convert back to string
println!("{}", key);

// Access individual components
println!("Standard: {}", key.is_standard());
println!("Version: {}", key.version());
println!("Protonation: {}", key.get_protonation());
```

### Binary Packing and Unpacking

```rust
use zinchi::InChIKey;

let key: InChIKey = "ZZJLMZYUGLJBSO-UHFFFAOYSA-N".parse()?;

// Pack to binary (9 or 14 bytes)
let packed = key.packed_bytes();
println!("Packed size: {} bytes", packed.len());

// Unpack from binary
let unpacked = InChIKey::unpack_from(&packed)?;
assert_eq!(key, unpacked);
```

### Working with Buffers

```rust
use zinchi::InChIKey;

let key: InChIKey = "ZZJLMZYUGLJBSO-UHFFFAOYSA-N".parse()?;

// Pack into an existing buffer
let mut buffer = [0u8; 14];
let size = key.pack_into(&mut buffer);

// Use only the relevant bytes
let packed_data = &buffer[..size];
```

## InChI Key Format

An InChI key has the format: `AAAAAAAAAAAAAA-BBBBBBBBFV-P`

- **First block** (14 chars): Encodes core molecular constitution (65 bits → 9 bytes)
- **Second block** (8 chars): Encodes stereochemistry and isotopes (37 bits → 5 bytes)
- **Flag** (1 char): `S` for standard, `N` for non-standard
- **Version** (1 char): Currently always `A`
- **Protonation** (1 char): `N` for neutral, or `A`-`M` for protonated states

## Binary Encoding

Standard InChI keys with the common second block `UHFFFAOYSA` (empty stereochemistry hash) are packed into just **9 bytes**. All other InChI keys require **14 bytes**.

This represents a **67-75% reduction** in size compared to the ASCII representation.

### Encoding Details

The first block (14 characters) is decoded into four 14-bit triples and one 9-bit pair, then packed into 9 bytes. The second block (8 characters) is decoded into two 14-bit triples and one 9-bit pair, then packed into 5 bytes. Additional metadata (standard flag, version, protonation) is encoded into spare bits.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- John Mayfield and NextMove Software for the original compression algorithm
- The InChI Trust for the InChI specification

## See Also

- [InChI Trust](https://www.inchi-trust.org/)
- [InChI Technical Manual](https://www.inchi-trust.org/download/104/InChI_TechMan.pdf)
- [Original compression presentation](https://www.nextmovesoftware.com/talks/Mayfield_DataCompressionOfInChIKeysAnd2dCoordinates_NIHINCHI_202103.pdf)
