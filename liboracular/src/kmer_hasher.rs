use std::cmp::min;
use std::{u64, i64, mem};

// TODO: Make sure this is optional on certain platforms
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

const FIND3: u64 = 7;

lazy_static! {

    static ref CONVERSION_TOCHAR: [char; 256] = {

        let mut conversion: [char; 256] = ['N'; 256];
        conversion[7]  = 'A';
        conversion[0]  = 'T';
        conversion[5]  = 'C';
        conversion[2]  = 'G';

        conversion
    };

    static ref CONVERSION: [u64; 256] = {

        let mut conversion: [u64; 256] = [1; 256];
        conversion[65]  = 7;
        conversion[97]  = 7;
        conversion[84]  = 0;
        conversion[116] = 0;
        conversion[67]  = 5;
        conversion[99]  = 5;
        conversion[71]  = 2;
        conversion[103] = 2;

        conversion
    };

    static ref CONVERSION_I: [i64; 256] = {

        let mut conversion: [i64; 256] = [1; 256];
        conversion[65]  = 7;
        conversion[97]  = 7;
        conversion[84]  = 0;
        conversion[116] = 0;
        conversion[67]  = 5;
        conversion[99]  = 5;
        conversion[71]  = 2;
        conversion[103] = 2;

        conversion
    };

    static ref SHIFT: __m128i = unsafe { _mm_set1_epi64x(3) };
}

// A is 7 is 111
// T is 0 is 000
// C is 5 is 101
// G is 2 is 010
// N is thus: 1 = 001 
// N is also: 3 = 011

#[inline(always)]
pub fn kmerhash(kmer: &[u8]) -> u64 {
    let mut bits: u64 = 0;
    bits = bits.wrapping_add(CONVERSION[usize::from(kmer[0])]);
    &kmer[1..].iter().for_each(|base| {
        bits <<= 3;
        bits = bits.wrapping_add(CONVERSION[usize::from(*base)]);
    });
    bits
}

#[inline(always)]
pub fn kmerhash_smallest(kmer: &[u8]) -> u64 {
    let hash = kmerhash(kmer);
    let rc = calc_rc(kmer.len(), hash);
    min(hash, rc)
}

#[inline(always)]
pub fn calc_rc(k: usize, khash: u64) -> u64 {
    // khash is a kmer already processed with kmerhash
    // k is the k in kmer (thus the seq length)
    let mut rc = !khash.reverse_bits();
    rc = rc >> (64 - (k * 3));
    replace3with1(rc)
    // rc
}

#[inline(always)]
pub fn hash4(k1: &[u8], k2: &[u8], k3: &[u8], k4: &[u8]) -> (u64, u64, u64, u64) {
    let hashes = _hash4(k1, k2, k3, k4);

    (min(hashes.0, calc_rc(k1.len(), hashes.0)),
     min(hashes.1, calc_rc(k2.len(), hashes.1)),
     min(hashes.2, calc_rc(k3.len(), hashes.2)),
     min(hashes.3, calc_rc(k4.len(), hashes.3)))
}

#[inline(always)]
fn replace3with1(mut x: u64) -> u64 {
    let is3: u64 = 3;
    let replace3: u64 = 2;
    let mut distance: u64;
    let mut query: u64;
    for i in 0..=20 {
        distance = i * 3;
        query = FIND3 << distance;
        if (x & query) == (is3 << distance) {
            x = x ^ (replace3 << distance);
        }
    }
    x
}

// AVX can calc 4 at a time
#[inline(always)]
fn _hash4(k1: &[u8], k2: &[u8], k3: &[u8], k4: &[u8]) -> (u64, u64, u64, u64) {
    let results: (u64, u64, u64, u64);

    unsafe {
        let mut hashes = _mm256_setzero_si256();

        let mut add = _mm256_set_epi64x(CONVERSION_I[usize::from(k1[0])],
                                        CONVERSION_I[usize::from(k2[0])],
                                        CONVERSION_I[usize::from(k3[0])],
                                        CONVERSION_I[usize::from(k4[0])]);
        
        hashes = _mm256_add_epi64(hashes, add);

        for i in 1..k1.len() {
            hashes = _mm256_sll_epi64(hashes, *SHIFT);
            add = _mm256_set_epi64x(CONVERSION_I[usize::from(k1[i])],
                                    CONVERSION_I[usize::from(k2[i])],
                                    CONVERSION_I[usize::from(k3[i])],
                                    CONVERSION_I[usize::from(k4[i])]);
            hashes = _mm256_add_epi64(hashes, add);
        }

        results = mem::transmute(hashes)
    }

    (results.3, results.2, results.1, results.0)
}


// Doesn't seem to speed it up, not gives correct result...
// Leaving here to fix in the future (maybe)
#[inline(always)]
fn _rc4(k: usize, k1: u64, k2: u64, k3: u64, k4: u64) -> (u64, u64, u64, u64) {
    let (ik1, ik2, ik3, ik4, x);
    let xs: (u64, u64, u64, u64);

    unsafe {
        ik1 = _bswap64(mem::transmute(k1));
        ik2 = _bswap64(mem::transmute(k2));
        ik3 = _bswap64(mem::transmute(k3));
        ik4 = _bswap64(mem::transmute(k4));
        x = _mm256_set_epi64x(ik1, ik2, ik3, ik4);

        let shift = _mm_set1_epi64x(mem::transmute(64 - (k * 3)));
        xs = mem::transmute(_mm256_sll_epi64(x, shift));
    }

    (replace3with1(xs.3), replace3with1(xs.2), replace3with1(xs.1), replace3with1(xs.0))
    
//    (0, 0, 0, 0)
}

#[inline(always)]
fn _get_pos(pos: usize, hash: u64) -> char {
    let distance = pos * 3;
    return CONVERSION_TOCHAR[((hash & (FIND3 << distance)) >> distance) as usize]
}

#[inline(always)]
pub fn convert_to_kmer(k: usize, hash: u64) -> String {
    let mut kmer: String = String::with_capacity(k);
    for i in (0..k).rev() {
        kmer.push(_get_pos(i, hash));
    }
    // return kmer[(21 - k)..].to_string()
    kmer
    
    /* let mut kmer: String = String::with_capacity(21);
    for i in (0..21).rev() {
        kmer.push(_get_pos(i, hash));
    }
    return kmer[(21 - k)..].to_string() */
}

#[cfg(test)]
mod test {

    use crate::kmer_hasher::{*};

    #[test]
    fn hash_one() {
        let kmer   = b"ACTCACGATCACGATACAAAN";
        let kmerrc = b"NTTTGTATCGTGATCGTGAGT";

        let hash   = kmerhash(kmer);
        let hashrc = kmerhash(kmerrc);
        assert_eq!(8804444413962215417, hash);
        assert_eq!(1153515602050352592, hashrc);
        assert_eq!(8804444413962215417, calc_rc(21, hashrc));
        assert_eq!(1153515602050352592, calc_rc(21, hash));
    }

    #[test]
    fn four_at_a_time() {
        let k = (b"ACTCACGATCACGATACAAAN", 
                 b"TCAGTCACTAGCATACAACTC",
                 b"ACGATCACGATACAAANNNNN",
                 b"TGACTANCATCANTACTTGGT");
        
        let krc = (b"NTTTGTATCGTGATCGTGAGT",
                   b"GAGTTGTATGCTAGTGACTGA",
                   b"NNNNNTTTGTATCGTGATCGT",
                   b"ACCAAGTANTGATGNTAGTCA");

        let mut hashes   = _hash4(k.0, k.0, k.0, k.0);
        let hashesrc     = _hash4(krc.0, krc.1, krc.2, krc.3);

        // 0b0111101000101111101010111000101111101 010 111 000 111 101 111 111 111 001
        // 0b0110101000100110101010111000100110101 010 111 000 110 100 110 110 111 010

        // println!("{:#066b}", hashes.0);
        println!("{:#066b}", hashes.0);
        println!("{:#066b}", calc_rc(21, hashesrc.0));

        assert_eq!(8804444413962215417, hashes.0);
        assert_eq!(hashes.0, calc_rc(21, hashesrc.0));

        hashes = _hash4(k.0, k.1, k.2, k.3);

        assert_eq!(8804444413962215417, hashes.0);
        assert_eq!(851389849605701445, hashes.1);
        assert_eq!(8843027523915715145, hashes.2);
        assert_eq!(425844089580060816, hashes.3);

        assert_eq!(8804444413962215417, calc_rc(21, hashesrc.0));
        assert_eq!(851389849605701445,  calc_rc(21, hashesrc.1));
        assert_eq!(8843027523915715145, calc_rc(21, hashesrc.2));
        assert_eq!(425844089580060816,  calc_rc(21, hashesrc.3));

        assert_eq!(1153515602050352592, hashesrc.0);
        assert_eq!(3350752362468833815, hashesrc.1);
        assert_eq!(1317584511025901904, hashesrc.2);
        assert_eq!(8898905677553234991, hashesrc.3);

        assert_eq!(1153515602050352592, calc_rc(21, hashes.0));
        assert_eq!(3350752362468833815, calc_rc(21, hashes.1));
        assert_eq!(1317584511025901904, calc_rc(21, hashes.2));
        assert_eq!(8898905677553234991, calc_rc(21, hashes.3));
    }

    #[test]
    fn four_at_a_time_vs_one_at_a_time() {
        let k = (b"ACTCACGATCACGATACAAAN", 
                 b"TCAGTCACTAGCATACAACTC",
                 b"ACGATCACGATACAAANNNNN",
                 b"TGACTANCATCANTACTTGGT");
        let results4 = _hash4(k.0, k.1, k.2, k.3);

        assert_eq!(results4.0, kmerhash(k.0));
        assert_eq!(results4.1, kmerhash(k.1));
        assert_eq!(results4.2, kmerhash(k.2));
        assert_eq!(results4.3, kmerhash(k.3));

        let krc = ("NTTTGTATCGTGATCGTGAGT",
                   "GAGTTGTATGCTAGTGACTGA",
                   "NNNNNTTTGTATCGTGATCGT",
                   "ACCAAGTANTGATGNTAGTCA");
        
        let results4rc = _rc4(21, results4.0, results4.1, results4.2, results4.3);

        // assert_eq!(results4rc.0, kmerhash(k.0));
    }

    #[test]
    fn hash4_test() {
        // Hash4 fn should return smallest hash
        let k = (b"ACTCACGATCACGATACAAAN", 
                 b"TCAGTCACTAGCATACAACTC",
                 b"ACGATCACGATACAAANNNNN",
                 b"TGACTANCATCANTACTTGGT");

        let krc = (b"NTTTGTATCGTGATCGTGAGT",
                   b"GAGTTGTATGCTAGTGACTGA",
                   b"NNNNNTTTGTATCGTGATCGT",
                   b"ACCAAGTANTGATGNTAGTCA");

        let hashes   = hash4(k.0, k.1, k.2, k.3);
        let hashesrc = hash4(krc.0, krc.1, krc.2, krc.3);

        let hashes_direct   = _hash4(k.0, k.1, k.2, k.3);
        let hashesrc_direct = _hash4(krc.0, krc.1, krc.2, krc.3);

        println!("{}\n{}", kmerhash(k.0), calc_rc(21, kmerhash(k.0)));
        println!("{}", hashes.0);

        assert_eq!(hashes.0, hashesrc.0, "Hashes 0 do not match");
        assert_eq!(hashes.1, hashesrc.1, "Hashes 1 do not match");
        assert_eq!(hashes.2, hashesrc.2, "Hashes 2 do not match");
        assert_eq!(hashes.3, hashesrc.3, "Hashes 3 do not match");

        assert_eq!(hashes.0, min(hashes_direct.0, hashesrc_direct.0), "Hashes 0 do not match");
        assert_eq!(hashes.1, min(hashes_direct.1, hashesrc_direct.1), "Hashes 1 do not match");
        assert_eq!(hashes.2, min(hashes_direct.2, hashesrc_direct.2), "Hashes 2 do not match");
        assert_eq!(hashes.3, min(hashes_direct.3, hashesrc_direct.3), "Hashes 3 do not match");
    }
}
