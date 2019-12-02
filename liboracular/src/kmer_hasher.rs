use std::cmp::min;
use std::{u64, i64, mem};

// TODO: Make sure this is optional on certain platforms
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

lazy_static! {
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

// A is 7
// T is 0
// C is 5
// G is 2
// N is thus: 1
// N is also: 4... I mean 3

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
    let rc = !khash.reverse_bits();
    rc >> (64 - (k * 3))
}

// TODO: Take advantage of this...
// AVX can calc 4 at a time
fn _hash4(k1: &[u8], k2: &[u8], k3: &[u8], k4: &[u8]) -> (u64, u64, u64, u64) {
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

        mem::transmute(hashes)
    }
}

fn _rc4(k: usize, k1: u64, k2: u64, k3: u64, k4: u64) -> (u64, u64, u64, u64) {
    let (ik1, ik2, ik3, ik4, x);

    unsafe {
        ik1 = _bswap64(mem::transmute(k1));
        ik2 = _bswap64(mem::transmute(k2));
        ik3 = _bswap64(mem::transmute(k3));
        ik4 = _bswap64(mem::transmute(k4));
        x = _mm256_set_epi64x(ik1, ik2, ik3, ik4);

        let shift = _mm_set1_epi64x(mem::transmute(64 - (k * 3)));
        mem::transmute(_mm256_sll_epi64(x, shift))
    }
    
//    (0, 0, 0, 0)
}