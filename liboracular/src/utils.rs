pub fn get_good_sequence_coords (seq: &[u8]) -> Vec<(usize, usize)> {
    let mut start: Option<usize> = None;
    let mut end: usize;
    let mut cur: usize = 0;
    let mut start_coords;
    let mut end_coords;
    let mut coords: Vec<(usize, usize)> = Vec::with_capacity(4096);
    //let results = seq.windows(3).enumerate().filter(|(_y, x)| x != &[78, 78, 78]).map(|(y, _x)| y);

    // Do we need to filter the sequence at all?
    if bytecount::count(&seq, b'N') < 3 {
        coords.push( (0, seq.len()) );
        return coords
    }

    let results = seq.windows(3).enumerate().filter(|(_y, x)| bytecount::count(&x, b'N') < 3).map(|(y, _x)| y);
    for pos in results {
        match start {
            None    => { start = Some(pos); cur = pos; }
            Some(_x) => ()
        };

        if pos - cur > 1 {
            end = cur;
            start_coords = start.unwrap();
            end_coords = end;
            coords.push( (start_coords, end_coords) );
            start = None;
        } else {
            cur = pos;
        }
    }

    // Push final set of coords to the system
    if start != None {
        end = cur;
        start_coords = start.unwrap();
        end_coords = end;
        if end_coords - start_coords > 1 {
            coords.push( (start_coords, end_coords) );
        }
    }

    coords
}

#[inline(always)]
fn _complement_nucl(nucl: u8) -> u8 {
    // Should all be capitalized by now...
    // N -> 78
    // A -> 65
    // C -> 67
    // G -> 71
    // T -> 84
    match &nucl {
        65 => 84, // A -> T
        67 => 71, // C -> G
        84 => 65, // T -> A
        71 => 67, // G -> C
        78 => 78, // Complement of N is N
        _ => 78,  // Everything else -> N
    }
}

// Mutability here because we change everything to uppercase
/// Complement nucleotides -- Reverse is easy enough with Rust internals
#[inline(always)]
pub fn complement_nucleotides(slice: &mut [u8]) {
    for x in slice.iter_mut() {
        *x = _complement_nucl(*x);
    }
}

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::prelude::*;
    use super::*;

    #[test]
    pub fn test_get_good_sequence_coords() {

        let coords = get_good_sequence_coords(b"AAAAAAAAAAAAAAAAAAAANNNAAAAAAAAAAAAAAAAAAAAAAAA");
        println!("{:#?}", coords);
        assert!(coords == [(0, 19), (22, 44)]);

        let coords = get_good_sequence_coords(b"AAAAAAAAAAAAAAAAAAAANNNAAAAAAAAAAAAAAAAAAAAAAAANNN");
        println!("{:#?}", coords);
        assert!(coords == [(0, 19), (22, 46)]);

        let coords = get_good_sequence_coords(b"NNNAAAAAAAAAAAAAAAAAAAANNNAAAAAAAAAAAAAAAAAAAAAAAANNN");
        println!("{:#?}", coords);
        assert!(coords == [(1, 22), (25, 49)]);

        let coords = get_good_sequence_coords(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        println!("{:#?}", coords);
        assert!(coords == [(0, 44)]);

    }

    #[test]
    pub fn test_complement_nucleotides() {
        let mut seq = b"AGTCCCNTNNNNTAAGATTTAGAGACCAAAAA".to_vec();
        complement_nucleotides(&mut seq);
        assert!(seq == b"TCAGGGNANNNNATTCTAAATCTCTGGTTTTT");
        seq.reverse();
        assert!(seq == b"TTTTTGGTCTCTAAATCTTANNNNANGGGACT");

    }

}
