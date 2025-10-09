//! A compact binary representation for InChI Keys.
//!
//! This crate provides a space-efficient binary encoding for International Chemical Identifier (InChI) keys,
//! reducing their size from the standard 27-byte ASCII representation to either 9 or 14 bytes. The implementation
//! is based on the work by John Mayfield (NextMove Software):
//! [Data Compression of InChI Keys and 2D Coordinates](https://www.nextmovesoftware.com/talks/Mayfield_DataCompressionOfInChIKeysAnd2dCoordinates_NIHINCHI_202103.pdf)
//!
//! # InChI Key Format
//!
//! An InChI key has the format: `AAAAAAAAAAAAAA-BBBBBBBBFV-P`
//! - First block (14 chars): Encoding core molecular constitution
//! - Second block (8 chars): Encoding advanced structural features whichever are applicable (stereochemistry,
//!   isotopic substitution, exact position of mobile hydrogens, metal ligation data)
//! - Flag (1 char): 'S' for standard, 'N' for non-standard
//! - Version (1 char): Currently always 'A'
//! - Protonation (1 char): 'N' for neutral, or 'A'-'M' for protonated states
//!
//! # Binary Encoding
//!
//! Standard InChI keys with the common second block `UHFFFAOYSA` can be packed into just 9 bytes.
//! All other InChI keys require 14 bytes.
//! 
//! # Optional Features
//!
//! - **`serde`**: Enable serialization/deserialization support. When enabled, `InChIKey` serializes
//!   as a string in human-readable formats (JSON, YAML) and uses the compact binary representation
//!   in binary formats (bincode, MessagePack).
//!
//!
//! # Example
//!
//! ```rust
//! use zinchi::InChIKey;
//!
//! let key: InChIKey = "ZZJLMZYUGLJBSO-UHFFFAOYSA-N".parse().unwrap();
//! let packed = key.packed_bytes(); // 9 or 14 bytes
//! let unpacked = InChIKey::unpack_from(&packed).unwrap();
//! assert_eq!(key, unpacked);
//! ```

use std::error::Error;
use std::fmt::{Display, Formatter, Write};
use std::str::FromStr;

/// The packed binary representation of the most common second block "UHFFFAOYSA".
/// This corresponds to structures with no stereochemistry or isotopes (empty hash).
const UHFFFAOYSA: [u8; 5] = [
    0xE3, // 227
    0xB0, // 176
    0xC4, // 196
    0x42, // 66
    0x18, // 24
];

/// Error type for InChI key parsing operations.
#[derive(Debug)]
pub enum InChIKeyParseError {
    /// The key contains invalid characters (non-uppercase letters where expected).
    InvalidCharacter,
    /// The key has incorrect length (expected 27 characters).
    InvalidLength,
    /// The standard/non-standard flag is invalid (expected 'S' or 'N').
    InvalidFlag,
    /// The version character is invalid (expected 'A').
    InvalidVersion,
    /// The protonation flag is invalid (expected 'N' or 'A'-'Z').
    InvalidProtonation,
}

impl Display for InChIKeyParseError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            InChIKeyParseError::InvalidCharacter => write!(f, "invalid character in InChIKey"),
            InChIKeyParseError::InvalidLength => write!(f, "invalid InChIKey length"),
            InChIKeyParseError::InvalidFlag => write!(f, "invalid standard/non-standard flag"),
            InChIKeyParseError::InvalidVersion => write!(f, "invalid InChI version"),
            InChIKeyParseError::InvalidProtonation => write!(f, "invalid protonation flag"),
        }
    }
}

impl Error for InChIKeyParseError {}

/// A compact binary representation of an InChI key.
///
/// This structure encodes an InChI key in 9 or 14 bytes depending on whether
/// it uses the common stereochemistry block "UHFFFAOYSA".
#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct InChIKey {
    /// Packed first block (connectivity hash) - 9 bytes encoding 14 characters.
    fst: [u8; 9],
    /// Packed second block (stereochemistry hash) - 5 bytes encoding 8 characters.
    snd: [u8; 5],
    /// Standard InChI flag (true = 'S', false = 'N').
    std: bool,
    /// InChI version character (currently always 'A' = 0x41).
    ver: u8,
    /// Protonation state ('N' for neutral or 'A'-'Z').
    chg: u8,
}

impl InChIKey {
    /// Returns the size in bytes of the packed binary representation.
    ///
    /// Returns 9 bytes for standard InChI keys with the common second block "UHFFFAOYSA",
    /// and 14 bytes for all other keys.
    pub fn get_packed_size(&self) -> usize {
        9 + 5 * !(self.std && self.snd == UHFFFAOYSA) as usize
    }

    /// Returns `true` if this is a standard InChI key.
    ///
    /// Standard InChI keys are generated with default options and have the flag 'S'.
    /// Non-standard keys (flag 'N') may have custom options applied.
    pub fn is_standard(&self) -> bool {
        self.std
    }

    /// Returns the InChI version character.
    ///
    /// Currently always returns `'A'` for version 1 of the InChI algorithm.
    pub fn version(&self) -> char {
        self.ver as char
    }

    /// Returns the protonation state character.
    ///
    /// Returns `'N'` for neutral molecules, or `'A'` through `'Z'` for various
    /// protonated or deprotonated states.
    pub fn get_protonation(&self) -> char {
        self.chg as char
    }

    /// Packs the InChI key into the provided buffer.
    ///
    /// Writes the binary representation into the buffer and returns the number of bytes written
    /// (either 9 or 14). The buffer must be at least 14 bytes long.
    ///
    /// # Arguments
    ///
    /// * `buffer` - A mutable buffer to write the packed bytes into
    ///
    /// # Returns
    ///
    /// The number of bytes written (9 or 14)
    ///
    /// # Panics
    ///
    /// Panics if the buffer is smaller than 14 bytes.
    ///
    /// # Example
    ///
    /// ```rust
    /// use zinchi::InChIKey;
    ///
    /// let key: InChIKey = "ZZZXOYCAMOTTNS-UHFFFAOYSA-N".parse().unwrap();
    /// let mut buffer = [0u8; 14];
    /// let size = key.pack_into(&mut buffer);
    /// assert_eq!(size, 9);
    /// ```
    pub fn pack_into<B: AsMut<[u8]>>(&self, buffer: &mut B) -> usize {
        let buf = buffer.as_mut();
        assert!(
            buf.len() >= 14,
            "pack_into buffer must be at least 14 bytes"
        );
        buf[..9].copy_from_slice(&self.fst);
        buf[8] |= (self.chg.wrapping_sub(b'A')) << 2;
        if self.get_packed_size() == 9 {
            buf[8] |= 0x02;
            return 9;
        }
        buf[9..14].copy_from_slice(&self.snd);
        if !self.std {
            buf[13] |= 1 << 5;
        }
        14
    }

    /// Returns a vector containing the packed binary representation.
    ///
    /// This is a convenience method that allocates a new vector and calls [`pack_into`](Self::pack_into).
    /// The returned vector will be 9 or 14 bytes depending on the InChI key.
    ///
    /// # Example
    ///
    /// ```rust
    /// use zinchi::InChIKey;
    ///
    /// let key: InChIKey = "ZZZXOYCAMOTTNS-UHFFFAOYSA-N".parse().unwrap();
    /// let packed = key.packed_bytes();
    /// assert_eq!(packed.len(), 9);
    /// ```
    pub fn packed_bytes(&self) -> Vec<u8> {
        let mut buf = vec![0u8; 14];
        let size = self.pack_into(&mut buf);
        buf.truncate(size);
        buf
    }

    /// Unpacks an InChI key from a binary buffer.
    ///
    /// Reconstructs an InChI key from its packed binary representation.
    /// The buffer must be exactly 9 or 14 bytes long.
    ///
    /// # Arguments
    ///
    /// * `buffer` - A buffer containing the packed binary data
    ///
    /// # Returns
    ///
    /// Returns `Ok(InChIKey)` on success, or an error if the buffer has invalid length.
    ///
    /// # Errors
    ///
    /// Returns [`InChIKeyParseError::InvalidLength`] if the buffer is not 9 or 14 bytes long.
    ///
    /// # Example
    ///
    /// ```rust
    /// use zinchi::InChIKey;
    ///
    /// let key: InChIKey = "ZZZXOYCAMOTTNS-UHFFFAOYSA-N".parse().unwrap();
    /// let packed = key.packed_bytes();
    /// let unpacked = InChIKey::unpack_from(&packed).unwrap();
    /// assert_eq!(key, unpacked);
    /// ```
    pub fn unpack_from<B: AsRef<[u8]> + ?Sized>(buffer: &B) -> Result<Self, InChIKeyParseError> {
        let buf = buffer.as_ref();
        if !(buf.len() == 9 || buf.len() == 14) {
            return Err(InChIKeyParseError::InvalidLength);
        }
        let mut fst = [0u8; 9];
        fst.copy_from_slice(&buf[..9]);
        let chg = b'A' + ((buf[8] >> 2) & 0x1F);
        fst[8] &= 0x01; // keep only e-bit
        let mut snd = [0u8; 5];
        let std = if ((buf[8] >> 1) & 0x1) == 1 {
            snd.copy_from_slice(&UHFFFAOYSA);
            true
        } else {
            if buf.len() < 14 {
                return Err(InChIKeyParseError::InvalidLength);
            }
            snd.copy_from_slice(&buf[9..14]);
            snd[4] &= 0x1F; // clear the std flag bit
            ((buf[13] >> 5) & 0x1) != 1
        };
        Ok(InChIKey {
            fst,
            snd,
            std,
            ver: b'A',
            chg,
        })
    }
}

impl FromStr for InChIKey {
    type Err = InChIKeyParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let inchikey = s.strip_prefix("InChIKey=").unwrap_or(s);
        let bytes = inchikey.as_bytes();
        validate_inchikey_bytes(bytes)?;
        let part_one: &[u8; 14] = bytes[..14].try_into().unwrap();
        let part_two: &[u8; 8] = bytes[15..23].try_into().unwrap();
        let mut fst = [0u8; 9];
        let mut snd = [0u8; 5];
        decode_part_one(&mut fst, part_one);
        decode_part_two(&mut snd, part_two);
        let std = bytes[23] == b'S';
        let ver = bytes[24];
        let chg = bytes[26];
        Ok(InChIKey {
            fst,
            snd,
            std,
            ver,
            chg,
        })
    }
}

impl Display for InChIKey {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let mut part_one = [0u8; 14];
        let mut part_two = [0u8; 8];
        encode_part_one(&self.fst, &mut part_one);
        encode_part_two(&self.snd, &mut part_two);
        for &b in &part_one {
            f.write_char(b as char)?;
        }
        f.write_char('-')?;
        for &b in &part_two {
            f.write_char(b as char)?;
        }
        f.write_char(if self.std { 'S' } else { 'N' })?;
        f.write_char(self.ver as char)?;
        f.write_char('-')?;
        f.write_char(self.chg as char)?;
        Ok(())
    }
}

#[cfg(feature = "serde")]
mod serde_impl {
    use super::{InChIKey, InChIKeyParseError};
    use serde::de::{self, Visitor};
    use serde::ser::Serializer;

    impl serde::Serialize for InChIKey {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            if serializer.is_human_readable() {
                serializer.collect_str(self)
            } else {
                let bytes = self.packed_bytes();
                serializer.serialize_bytes(&bytes)
            }
        }
    }

    struct InChIKeyVisitor;

    impl<'de> Visitor<'de> for InChIKeyVisitor {
        type Value = InChIKey;

        fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(formatter, "an InChIKey string or packed InChIKey bytes")
        }

        fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            v.parse()
                .map_err(|e: InChIKeyParseError| de::Error::custom(e.to_string()))
        }

        fn visit_borrowed_str<E>(self, v: &'de str) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            self.visit_str(v)
        }

        fn visit_bytes<E>(self, v: &[u8]) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            InChIKey::unpack_from(v)
                .map_err(|e: InChIKeyParseError| de::Error::custom(e.to_string()))
        }

        fn visit_byte_buf<E>(self, v: Vec<u8>) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            self.visit_bytes(&v)
        }
    }

    impl<'de> serde::Deserialize<'de> for InChIKey {
        fn deserialize<D>(deserializer: D) -> Result<InChIKey, D::Error>
        where
            D: serde::Deserializer<'de>,
        {
            if deserializer.is_human_readable() {
                deserializer.deserialize_str(InChIKeyVisitor)
            } else {
                deserializer.deserialize_bytes(InChIKeyVisitor)
            }
        }
    }

    impl serde::Serialize for InChIKeyParseError {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: Serializer,
        {
            let s = format!("{}", self);
            serializer.serialize_str(&s)
        }
    }
}

/// Decodes a pair of uppercase letters into a 9-bit integer.
///
/// Maps two letters 'AA'-'ZZ' to values 0-675 (26²-1).
#[inline(always)]
fn decode_pair([a, b]: [u8; 2]) -> u16 {
    ((a - b'A') as u16) * 26 + ((b - b'A') as u16)
}

/// Decodes a triple of uppercase letters into a 14-bit integer.
///
/// Maps three letters 'AAA'-'ZZZ' to values 0-16383, skipping invalid ranges
/// used in InChI encoding (values that would decode to strings containing 'AAA').
#[inline(always)]
fn decode_triple([a, b, c]: [u8; 3]) -> u16 {
    let mut i = ((a - b'A') as u16) * 676 + ((b - b'A') as u16) * 26 + (c - b'A') as u16;
    i -= ((i >= 12844) as u16) * 516;
    i -= ((i >= 2704) as u16) * 676;
    i
}

/// Decodes the first block of an InChI key (14 characters) into 9 bytes.
///
/// The first block consists of four triples (12 chars) plus one pair (2 chars).
/// These are decoded and packed into 9 bytes using bit manipulation.
#[inline(always)]
fn decode_part_one(buffer: &mut [u8; 9], s: &[u8; 14]) {
    let [t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, p0, p1] = *s;
    let a = decode_triple([t0, t1, t2]);
    let b = decode_triple([t3, t4, t5]);
    let c = decode_triple([t6, t7, t8]);
    let d = decode_triple([t9, t10, t11]);
    let e = decode_pair([p0, p1]);
    *buffer = [
        a as u8,
        ((a >> 8) & 0x3f) as u8 | ((b << 6) & 0xc0) as u8,
        (b >> 2) as u8,
        ((b >> 10) & 0x0f) as u8 | ((c << 4) & 0xf0) as u8,
        (c >> 4) as u8,
        ((c >> 12) & 0x03) as u8 | ((d << 2) & 0xfc) as u8,
        (d >> 6) as u8,
        e as u8,
        ((e >> 8) & 0x01) as u8,
    ];
}

/// Decodes the second block of an InChI key (8 characters) into 5 bytes.
///
/// The second block consists of two triples (6 chars) plus one pair (2 chars).
/// These are decoded and packed into 5 bytes using bit manipulation.
#[inline(always)]
fn decode_part_two(buffer: &mut [u8; 5], s: &[u8; 8]) {
    let [t0, t1, t2, t3, t4, t5, p0, p1]: [u8; 8] = *s;
    let a = decode_triple([t0, t1, t2]);
    let b = decode_triple([t3, t4, t5]);
    let c = decode_pair([p0, p1]);
    *buffer = [
        a as u8,
        ((a >> 8) & 0x3f) as u8 | ((b << 6) & 0xc0) as u8,
        (b >> 2) as u8,
        ((b >> 10) & 0x0f) as u8 | ((c << 4) & 0xf0) as u8,
        ((c >> 4) & 0x1f) as u8,
    ];
}

/// Encodes a 9-bit integer into a pair of uppercase letters.
///
/// Maps values 0-675 to two letters 'AA'-'ZZ'. The inverse of [`decode_pair`].
#[inline(always)]
fn encode_pair(i: u16, out: &mut [u8]) {
    let i = i & 0x1ff;
    out[1] = b'A' + (i % 26) as u8;
    out[0] = b'A' + ((i / 26) % 26) as u8;
}

/// Encodes a 14-bit integer into a triple of uppercase letters.
///
/// Maps values 0-16383 to three letters 'AAA'-'ZZZ', adding gaps to match
/// InChI encoding. The inverse of [`decode_triple`].
#[inline(always)]
fn encode_triple(i: u16, out: &mut [u8]) {
    let mut i = i & 0x3fff;
    i += ((i >= 2704) as u16) * 676;
    i += ((i >= 12844) as u16) * 516;
    out[2] = b'A' + (i % 26) as u8;
    i /= 26;
    out[1] = b'A' + (i % 26) as u8;
    out[0] = b'A' + ((i / 26) % 26) as u8;
}

/// Encodes 9 bytes into the first block of an InChI key (14 characters).
///
/// Unpacks the 9 bytes into four triples (12 chars) and one pair (2 chars)
/// using bit manipulation, then encodes each as uppercase letters.
#[inline(always)]
fn encode_part_one(buffer: &[u8; 9], out: &mut [u8; 14]) {
    let a = (buffer[0] as u32) | ((buffer[1] as u32 & 0x3f) << 8);
    let b =
        ((buffer[1] as u32 & 0xc0) | ((buffer[2] as u32) << 8) | ((buffer[3] as u32 & 0x0f) << 16))
            >> 6;
    let c =
        ((buffer[3] as u32 & 0xf0) | ((buffer[4] as u32) << 8) | ((buffer[5] as u32 & 0x03) << 16))
            >> 4;
    let d = ((buffer[5] as u32 & 0xfc) | ((buffer[6] as u32) << 8)) >> 2;
    let e = (buffer[7] as u32) | ((buffer[8] as u32 & 0x01) << 8);
    encode_triple(a as u16, &mut out[0..3]);
    encode_triple(b as u16, &mut out[3..6]);
    encode_triple(c as u16, &mut out[6..9]);
    encode_triple(d as u16, &mut out[9..12]);
    encode_pair(e as u16, &mut out[12..14]);
}

/// Encodes 5 bytes into the second block of an InChI key (8 characters).
///
/// Unpacks the 5 bytes into two triples (6 chars) and one pair (2 chars)
/// using bit manipulation, then encodes each as uppercase letters.
#[inline(always)]
fn encode_part_two(buffer: &[u8; 5], out: &mut [u8; 8]) {
    let a = (buffer[0] as u32) | ((buffer[1] as u32 & 0x3f) << 8);
    let b =
        ((buffer[1] as u32 & 0xc0) | ((buffer[2] as u32) << 8) | ((buffer[3] as u32 & 0x0f) << 16))
            >> 6;
    let c = ((buffer[3] as u32 & 0xf0) | ((buffer[4] as u32 & 0x1f) << 8)) >> 4;
    encode_triple(a as u16, &mut out[0..3]);
    encode_triple(b as u16, &mut out[3..6]);
    encode_pair(c as u16, &mut out[6..8]);
}

/// Validates the format of an InChI key byte string.
///
/// Checks that the byte string:
/// - Has exactly 27 characters
/// - Contains hyphens at positions 14 and 25
/// - Contains only uppercase letters A-Z in the hash blocks
/// - Has a valid standard/non-standard flag ('S' or 'N')
/// - Has version 'A'
/// - Has a valid protonation flag ('N' or 'A'-'Z')
///
/// # Errors
///
/// Returns an error if any validation check fails.
fn validate_inchikey_bytes(bytes: &[u8]) -> Result<(), InChIKeyParseError> {
    if bytes.len() != 27 {
        return Err(InChIKeyParseError::InvalidLength);
    }
    let b = bytes;
    if b[14] != b'-' || b[25] != b'-' {
        return Err(InChIKeyParseError::InvalidCharacter);
    }
    if !b[0..14].iter().all(|&c| c.is_ascii_uppercase()) {
        return Err(InChIKeyParseError::InvalidCharacter);
    }
    if !b[15..23].iter().all(|&c| c.is_ascii_uppercase()) {
        return Err(InChIKeyParseError::InvalidCharacter);
    }
    match b[23] {
        b'S' | b'N' => {}
        _ => return Err(InChIKeyParseError::InvalidFlag),
    }
    if b[24] != b'A' {
        return Err(InChIKeyParseError::InvalidVersion);
    }
    let p = b[26];
    if !(p == b'N' || (b'A'..=b'M').contains(&p)) {
        return Err(InChIKeyParseError::InvalidProtonation);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use quickcheck::{quickcheck, Arbitrary, Gen};

    const TEST_KEYS: &[&str] = &[
        "ZZZXOYCAMOTTNS-UHFFFAOYSA-N",
        "ZZQLONWLKXZCDD-UHFFFAOYSA-N",
        "ZZJLMZYUGLJBSO-UHFFFAOYSA-N",
        "ZZJLMZYUGLJBSO-LAEOZQHASA-N",
        "ZTTQDQPASSNLFS-ONEGZZNKSA-N",
    ];

    const INVALID_KEYS: &[&str] = &[
        "SHORT-KEY-N",
        "TOO-LONG-KEY-ABCDEFGHIJKLMNOPQRSTUVWXYZ",
        "INVALID#CHARACTERS-AAAAAA-N",
        "",
    ];

    #[test]
    fn test_from_str_roundtrip() {
        for &key in TEST_KEYS {
            let ik: InChIKey = key.parse().expect("Failed to parse InChIKey");
            let s = ik.to_string();
            assert_eq!(s, key, "Roundtrip string mismatch for {}", key);
        }
    }

    #[test]
    fn test_invalid_parse() {
        for &key in INVALID_KEYS {
            let res: Result<InChIKey, _> = key.parse();
            assert!(res.is_err(), "Invalid key parsed successfully: {}", key);
        }
    }

    #[test]
    fn test_packed_bytes_roundtrip() {
        for &key in TEST_KEYS {
            let ik: InChIKey = key.parse().expect("Failed to parse InChIKey");
            let packed = ik.packed_bytes();
            let unpacked = InChIKey::unpack_from(&packed).unwrap();
            assert_eq!(ik, unpacked, "Packed/unpacked mismatch for {}", key);
        }
    }

    impl Arbitrary for InChIKey {
        fn arbitrary(g: &mut Gen) -> Self {
            let mut fst = [0u8; 9];
            let mut snd = [0u8; 5];
            for byte in &mut fst {
                *byte = *g.choose(&(0..=255).collect::<Vec<_>>()).unwrap();
            }
            fst[8] &= 0x01; // only bit 0 is used (the e-bit)
            for byte in &mut snd {
                *byte = *g.choose(&(0..=255).collect::<Vec<_>>()).unwrap();
            }
            snd[4] &= 0x1F; // only bits 0-4 are used
            let std = *g.choose(&[true, false]).unwrap();
            let chg = if *g.choose(&[true, false]).unwrap() {
                b'N'
            } else {
                *g.choose(&(b'A'..=b'M').collect::<Vec<_>>()).unwrap()
            };
            InChIKey {
                fst,
                snd,
                std,
                ver: b'A',
                chg,
            }
        }
    }

    fn prop_string_roundtrip(key: InChIKey) -> bool {
        let s = key.to_string();
        let parsed: InChIKey = s.parse().unwrap();
        key == parsed
    }

    fn prop_packed_roundtrip(key: InChIKey) -> bool {
        let packed = key.packed_bytes();
        let unpacked = InChIKey::unpack_from(&packed).unwrap();
        key == unpacked
    }

    fn prop_packed_size(key: InChIKey) -> bool {
        let size = key.get_packed_size();
        size == 9 || size == 14
    }

    #[test]
    fn quickcheck_properties() {
        quickcheck(prop_string_roundtrip as fn(InChIKey) -> bool);
        quickcheck(prop_packed_roundtrip as fn(InChIKey) -> bool);
        quickcheck(prop_packed_size as fn(InChIKey) -> bool);
    }

    #[cfg(feature = "serde")]
    mod serde_tests {
        use super::*;

        #[test]
        fn bincode_roundtrip() {
            for &key_str in TEST_KEYS {
                let key: InChIKey = key_str.parse().unwrap();
                let encoded =
                    bincode::serde::encode_to_vec(&key, bincode::config::standard()).unwrap();
                let (decoded, _): (InChIKey, usize) =
                    bincode::serde::decode_from_slice(&encoded, bincode::config::standard())
                        .unwrap();
                assert_eq!(key, decoded, "Bincode roundtrip failed for {}", key_str);
            }
        }

        #[test]
        fn json_roundtrip() {
            for &key_str in TEST_KEYS {
                let key: InChIKey = key_str.parse().unwrap();
                let json = serde_json::to_string(&key).unwrap();
                let decoded: InChIKey = serde_json::from_str(&json).unwrap();
                assert_eq!(key, decoded, "JSON roundtrip failed for {}", key_str);
                assert_eq!(json, format!("\"{}\"", key_str), "JSON format check failed");
            }
        }
    }
}
