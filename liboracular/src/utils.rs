pub fn get_good_sequence_coords (seq: &[u8]) -> Vec<(usize, usize)> {
    let mut start: Option<usize> = None;
    let mut end: usize;
    let mut cur: usize = 0;
    let mut start_coords;
    let mut end_coords;
    let mut coords: Vec<(usize, usize)> = Vec::with_capacity(8);
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

    coords
}