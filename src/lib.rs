use pyo3::prelude::*;
use pyo3::types::PyAny;
use numpy::ndarray::{Array2, Axis};
use std::collections::{HashMap, HashSet};
use md5;
use log::debug; // Add logging

#[pyfunction]
pub fn build_one_full_read_no_extended_window_mode_rust(
    first_aligned_pos: usize,
    _last_aligned_pos: usize,
    mut full_read: Vec<String>,
    indel_map: &PyAny,  // Accept any Python iterable
    insertion_char: String,
    mut full_qualities: Option<Vec<String>>,
    read_query_name: Option<String>,
    full_read_cigar_hash: Option<String>,
) -> PyResult<(String, Option<Vec<String>>)> {
    if insertion_char != "-" && insertion_char != "X" {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Illegal char representation insertion",
        ));
    }

    let indel_iter = indel_map.iter()?;
    for indel in indel_iter {
        let (name, start, cigar_hash, ref_pos, indel_len, is_del): (
            String,
            usize,
            String,
            usize,
            usize,
            usize,
        ) = indel?.extract()?;

        if read_query_name.as_ref() == Some(&name)
            && start == first_aligned_pos
            && full_read_cigar_hash.as_deref() == Some(cigar_hash.as_str())
        {
            if indel_len > 0 {
                let pos = (ref_pos + 1)
                    .checked_sub(first_aligned_pos)
                    .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyIndexError, _>("Invalid position"))?;

                if pos >= full_read.len() {
                    return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                        "Position out of bounds during deletion",
                    ));
                }

                for _ in 0..indel_len {
                    full_read.remove(pos);
                    if let Some(ref mut quals) = full_qualities {
                        quals.remove(pos);
                    }
                }
            }

            if is_del == 1 {
                let pos = ref_pos
                    .checked_sub(first_aligned_pos)
                    .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyIndexError, _>("Invalid position"))?;

                if pos > full_read.len() {
                    return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                        "Position out of bounds during insertion",
                    ));
                }

                full_read.insert(pos, "-".to_string());
                if let Some(ref mut quals) = full_qualities {
                    quals.insert(pos, "2".to_string());
                }
            }
        }
    }

    Ok((full_read.join(""), full_qualities))
}


#[pyfunction]
pub fn build_one_full_read_with_extended_window_mode_rust(
    first_aligned_pos: usize,
    last_aligned_pos: usize,
    mut full_read: Vec<String>,
    indel_map: &PyAny,
    max_ins_at_pos: &PyAny,
    insertion_char: String,
    mut full_qualities: Option<Vec<String>>,
    read_query_name: Option<String>,
    full_read_cigar_hash: Option<String>,
) -> PyResult<(String, Option<Vec<String>>)> {
    if insertion_char != "-" && insertion_char != "X" {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Illegal char representation insertion",
        ));
    }

    let mut all_inserts = HashMap::new();
    let mut own_inserts = HashSet::new();
    let mut change_in_reference_space_ins = 0;

    // Convert max_ins_at_pos Python dict to Rust HashMap
    let max_ins: HashMap<usize, usize> = max_ins_at_pos.extract()?;

    let indel_iter = indel_map.iter()?;
    for indel in indel_iter {
        let (name, start, cigar_hash, ref_pos, indel_len, is_del): (String, usize, String, usize, usize, usize) = indel?.extract()?;

        if read_query_name.as_ref() == Some(&name)
            && start == first_aligned_pos
            && full_read_cigar_hash.as_deref() == Some(cigar_hash.as_str())
        {
            if indel_len > 0 {
                own_inserts.insert((ref_pos, indel_len));
                change_in_reference_space_ins += indel_len;
                if let Some(&current_max) = all_inserts.get(&ref_pos) {
                    if max_ins[&ref_pos] > current_max {
                        all_inserts.insert(ref_pos, max_ins[&ref_pos]);
                    }
                } else {
                    all_inserts.insert(ref_pos, max_ins[&ref_pos]);
                }
            }

            if is_del == 1 {
                let pos = ref_pos - first_aligned_pos + change_in_reference_space_ins;
                if pos > full_read.len() {
                    return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                        "Position out of bounds during insertion",
                    ));
                }
                full_read.insert(pos, insertion_char.clone());
                if let Some(ref mut quals) = full_qualities {
                    quals.insert(pos, "2".to_string());
                }
            }
        } else if first_aligned_pos <= ref_pos && ref_pos <= last_aligned_pos && indel_len > 0 {
            all_inserts.insert(ref_pos, max_ins[&ref_pos]);
        }
    }

    let mut change_in_reference_space = 0;
    let own_inserts_vec: Vec<(usize, usize)> = own_inserts.into_iter().collect();
    let (own_inserts_pos, own_inserts_len): (Vec<usize>, Vec<usize>) = if !own_inserts_vec.is_empty() {
        let (pos, len) = own_inserts_vec.into_iter().unzip();
        (pos, len)
    } else {
        (Vec::new(), Vec::new())
    };

    let mut sorted_positions: Vec<usize> = all_inserts.keys().cloned().collect();
    sorted_positions.sort_unstable();

    for pos in sorted_positions {
        let n = all_inserts[&pos];
        if own_inserts_pos.contains(&pos) {
            change_in_reference_space += n;
            continue;
        }

        let mut l = max_ins[&pos];
        let mut in_idx = pos + 1 - first_aligned_pos + change_in_reference_space;

        if let Some(index) = own_inserts_pos.iter().position(|&p| p == pos) {
            let k = own_inserts_len[index];
            l -= k;
            in_idx += k;
        }

        for _ in 0..l {
            if in_idx > full_read.len() {
                return Err(PyErr::new::<pyo3::exceptions::PyIndexError, _>(
                    "Position out of bounds during insertion",
                ));
            }
            full_read.insert(in_idx, insertion_char.clone());
            if let Some(ref mut quals) = full_qualities {
                quals.insert(in_idx, "2".to_string());
            }
        }
        change_in_reference_space += max_ins[&pos];
    }

    Ok((full_read.join(""), full_qualities))
}


fn parse_cigar(cigarstring: &str) -> PyResult<Vec<(u32, usize)>> {
    let mut cigartuples = Vec::new();
    let mut num = String::new();

    for c in cigarstring.chars() {
        if c.is_digit(10) {
            num.push(c);
        } else {
            let count = num.parse::<usize>().map_err(|_| {
                PyErr::new::<pyo3::exceptions::PyValueError, _>("Invalid CIGAR string")
            })?;
            let op = match c {
                'M' => 0,
                'I' => 1,
                'D' => 2,
                'S' => 4,
                'H' => 5,
                '=' => 7,
                'X' => 8,
                _ => return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>("Unsupported CIGAR operation")),
            };
            cigartuples.push((op, count));
            num.clear();
        }
    }

    Ok(cigartuples)
}


#[pyfunction]
fn run_one_window_rust(
    reads: Vec<(String, usize, usize, String, Option<Vec<u8>>, String)>, // (qname, ref_start, ref_end, seq, quals, cigarstring)
    reference_name: &str,
    window_start: usize, // 0- based
    window_length: usize,
    control_window_length: usize,
    mut minimum_overlap: usize,
    mut permitted_reads_per_location: HashMap<usize, usize>,
    exact_conformance_fix_0_1_basing_in_reads: bool,
    indel_map: &PyAny,
    max_ins_at_pos: &PyAny,
    extended_window_mode: bool,
    exclude_non_var_pos_threshold: f64,
    counter: usize,
) -> PyResult<(Vec<String>, Vec<Option<Vec<u8>>>, Vec<(String, usize, usize, String)>, Option<Vec<bool>>)> {
    let mut arr = Vec::new();
    let mut arr_read_summary = Vec::new();
    let mut arr_read_qualities_summary = Vec::new();

    let mut base_pair_distr_in_window: Option<Array2<usize>> = None;
    let alphabet = b"ACGT-NX";
    if exclude_non_var_pos_threshold > 0.0 {
        let alphabet = b"ACGT-NX";
        base_pair_distr_in_window = Some(
            Array2::<usize>::zeros((window_length, alphabet.len())),
        );
    }

    let original_window_length = window_length;
    let original_minimum_overlap = minimum_overlap;
    let mut window_length = control_window_length;

    if extended_window_mode {
        minimum_overlap = (minimum_overlap as f64 * window_length as f64 / original_window_length as f64).floor() as usize;
    }

    let mut debug_counter = 0;
    for (qname, first_aligned_pos, ref_end, seq, quals, cigarstring) in reads.iter() {
        let last_aligned_pos = *ref_end - 1; //somehow not needed
        if permitted_reads_per_location.get(&first_aligned_pos).unwrap_or(&0) == &0 {
            continue;
        } else {
            *permitted_reads_per_location.get_mut(&first_aligned_pos).unwrap() -= 1;
        }

        // Adjust last_aligned_pos to match Python's logic (reference_end is exclusive)
        let first_aligned_pos_ref = *first_aligned_pos;

        // Parse CIGAR string into tuples
        let read_cigartuples: Vec<(u32, usize)> = parse_cigar(cigarstring)?;

        // Convert sequence and qualities
        let mut processed_read_cp = String::new();
        let mut full_read: Vec<String> = seq.chars().map(|c| c.to_string()).collect();
        let mut full_qualities: Option<Vec<String>> =
            quals.clone().map(|v| v.into_iter().map(|q| q.to_string()).collect());

        // Handle soft clipping and other CIGAR operations
        for (ct_idx, ct) in read_cigartuples.iter().enumerate() {
            match ct.0 {
                0 | 1 | 2 | 7 | 8 => {}, // Match, insertion, deletion, etc.

                // Soft clipping
                4 => {
                    let num_clipped = ct.1;
                    for _ in 0..num_clipped {
                        let k = if ct_idx == 0 { 0 } else { full_read.len() - 1 };
                        full_read.remove(k);
                        if let Some(ref mut quals) = full_qualities {
                            quals.remove(k);
                        }
                    }

                    if ct_idx != 0 && ct_idx != read_cigartuples.len() - 1 {
                        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
                            "Soft clipping only possible on the edges of a read.",
                        ));
                    }
                },

                // Hard clipping - ignore
                5 => {},

                // Unsupported operations
                _ => {
                    return Err(PyErr::new::<pyo3::exceptions::PyNotImplementedError, _>(
                        format!("CIGAR op code {} not implemented", ct.0),
                    ));
                }
            }
        }

        // Process indels and build the full read
        let full_read_cigar_hash = format!("{:x}", md5::compute(cigarstring.as_bytes()));
        if extended_window_mode {
            let (processed_read, processed_qualities) =
                build_one_full_read_with_extended_window_mode_rust(
                    *first_aligned_pos,
                    last_aligned_pos,
                    full_read.clone(),
                    indel_map,
                    max_ins_at_pos,
                    "-".to_string(),
                    full_qualities.clone(),
                    Some(qname.clone()),
                    Some(full_read_cigar_hash),
                )?;
            processed_read_cp = processed_read.clone();
            full_read = processed_read.chars().map(|c| c.to_string()).collect();
            full_qualities = processed_qualities;
        } else {
            let (processed_read, processed_qualities) =
                build_one_full_read_no_extended_window_mode_rust(
                    *first_aligned_pos,
                    last_aligned_pos,
                    full_read.clone(),
                    indel_map,
                    "-".to_string(),
                    full_qualities.clone(),
                    Some(qname.clone()),
                    Some(full_read_cigar_hash),
                )?;
            processed_read_cp = processed_read.clone();
            full_read = processed_read.chars().map(|c| c.to_string()).collect();
            full_qualities = processed_qualities;
        }

        if first_aligned_pos + minimum_overlap < window_start + 1 + window_length
            && last_aligned_pos  >= window_start + original_minimum_overlap - 2
            && processed_read_cp.len() >= minimum_overlap
        {
            debug_counter+=1;

            let mut num_inserts_right_of_read = 0;
            let mut num_inserts_left_of_read = 0;
            let mut num_inserts_left_of_window = 0;

            if extended_window_mode {
                for (pos, val) in max_ins_at_pos.extract::<HashMap<usize, usize>>()?.iter() {
                    if last_aligned_pos < *pos && *pos < window_start + original_window_length {
                        num_inserts_right_of_read += val;
                    }
                    if window_start <= *pos && *pos < *first_aligned_pos {
                        num_inserts_left_of_read += val;
                    }
                    if *first_aligned_pos <= *pos && *pos < window_start {
                        num_inserts_left_of_window += val;
                    }
                }
            }

            // Calculate slice positions
            let start_cut_out = (window_start as isize + num_inserts_left_of_window as isize - *first_aligned_pos as isize).max(0);
            let end_cut_out = start_cut_out + (window_length as isize - num_inserts_left_of_read as isize);

            // Clamp indices to valid ranges
            let start_cut_out_clamped = start_cut_out.min(full_read.len() as isize).max(0);
            let end_cut_out_clamped = end_cut_out.min(full_read.len() as isize).max(0);

            // Slice the read and qualities
            let mut cut_out_read: Vec<String> =
                full_read[start_cut_out_clamped as usize..end_cut_out_clamped as usize].to_vec();
            let mut cut_out_qualities: Option<Vec<u8>> =
                if let Some(quals) = &full_qualities {
                    Some(
                        quals[start_cut_out_clamped as usize..end_cut_out_clamped as usize]
                            .iter()
                            .map(|q| q.parse::<u8>().unwrap_or(2))
                            .collect(),
                    )
                } else {
                    None
                };

            // Apply LEFT padding if start_cut_out is negative
            let original_start_cut_out = window_start as isize + num_inserts_left_of_window as isize - *first_aligned_pos as isize;
            if original_start_cut_out < 0 {
                let padding_left = (-original_start_cut_out) as usize + num_inserts_left_of_read;
                cut_out_read.splice(0..0, vec!["N".to_string(); padding_left]);
                if let Some(ref mut quals) = cut_out_qualities {
                    quals.splice(0..0, vec![2; padding_left]);
                }
            }

            // Apply RIGHT padding (k)
            let k = (
                (window_start as isize + window_length as isize - 1  // Use current window_length
                - last_aligned_pos as isize
                + num_inserts_right_of_read as isize)
            ).max(0) as usize;

            if k > 0 {
                cut_out_read.extend(vec!["N".to_string(); k]);
                if let Some(ref mut quals) = cut_out_qualities {
                    quals.extend(vec![2; k]);
                }
            }

            // Final length adjustment
            let current_len = cut_out_read.len();
            if current_len < window_length {
                let needed = window_length - current_len;
                cut_out_read.extend(vec!["N".to_string(); needed]);
                if let Some(ref mut quals) = cut_out_qualities {
                    quals.extend(vec![2; needed]);
                }
            } else if current_len > window_length {
                cut_out_read.truncate(window_length);
                if let Some(ref mut quals) = cut_out_qualities {
                    quals.truncate(window_length);
                }
            }

            assert_eq!(cut_out_read.len(), window_length);

            assert_eq!(
                cut_out_read.len(),
                window_length,
                "cut_out_read length {} != window_length {}",
                cut_out_read.len(),
                window_length
            );
            //println!(" AFTER Final adjustment to enforce window_lengthcut_out_read.len() {}", cut_out_read.len());
            // Final assertions
            assert_eq!(cut_out_read.len(), window_length,
                "cut_out_read length {} != window_length {}",
                cut_out_read.len(), window_length
            );
            if let Some(ref quals) = cut_out_qualities {
                assert_eq!(quals.len(), window_length,
                    "cut_out_qualities length {} != window_length {}",
                    quals.len(), window_length
                );
            }

            // Add to output arrays
            arr.push(format!(
                ">{} {}\n{}",
                qname,
                first_aligned_pos + if exact_conformance_fix_0_1_basing_in_reads { 1 } else { 0 },
                cut_out_read.join("")
            ));

            // Update base pair distribution
            if exclude_non_var_pos_threshold > 0.0 {
                if let Some(base_pair_distr_in_window) = &mut base_pair_distr_in_window {
                    for (idx, letter) in cut_out_read.iter().enumerate() {
                        if let Some(pos) = alphabet.iter().position(|&x| x == letter.as_bytes()[0]) {
                            base_pair_distr_in_window[[idx, pos]] += 1;
                        }
                    }
                }
            }

            // Store qualities if present
            if let Some(quals) = cut_out_qualities {
                arr_read_qualities_summary.push(Some(quals));
            } else {
                arr_read_qualities_summary.push(None);
            }

        }
        //Store read summary
        let counter_ref = &counter;
        if &first_aligned_pos_ref >= counter_ref && processed_read_cp.len() >= minimum_overlap{
            if qname == "89.6-3332" {
                println!("in printing -cluase first_aligned_pos_ref: {}",first_aligned_pos_ref);
            }
            arr_read_summary.push((qname.clone(), *first_aligned_pos + 1 , *ref_end, processed_read_cp.clone()));
        }

    }

    // Compute positional filter if applicable
    let mut pos_filter: Option<Vec<bool>> = None;
    if exclude_non_var_pos_threshold > 0.0 && !arr.is_empty() {
        if let Some(base_pair_distr_in_window) = &base_pair_distr_in_window {
            let mut max_freqs = Vec::new();
            let mut total_counts = Vec::new();

            for row in base_pair_distr_in_window.axis_iter(Axis(0)) {
                let max_freq = *row.iter().max().unwrap();
                let total_count = row.sum();
                max_freqs.push(max_freq);
                total_counts.push(total_count);
            }

            pos_filter = Some(
                max_freqs
                    .into_iter()
                    .zip(total_counts)
                    .map(|(max_freq, total_count)| {
                        if total_count == 0 {
                            true // Keep positions with no reads
                        } else {
                            1.0 - (max_freq as f64 / total_count as f64) >= exclude_non_var_pos_threshold
                        }
                    })
                    .collect(),
            );

            // Apply pos_filter to each read
            if let Some(pos_filter) = &pos_filter {
                for quals in arr_read_qualities_summary.iter_mut() {
                    if let Some(quals) = quals {
                        let filtered_quals: Vec<u8> = quals.iter()
                            .enumerate()
                            .filter_map(|(idx, q)| if pos_filter[idx] { Some(*q) } else { None })
                            .collect();
                        *quals = filtered_quals;
                    }
                }
            }


        }
    }

    Ok((arr, arr_read_qualities_summary, arr_read_summary, pos_filter))
}


#[pymodule]
pub fn viloca_rust(_py: Python<'_>, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(build_one_full_read_no_extended_window_mode_rust, m)?)?;
    m.add_function(wrap_pyfunction!(build_one_full_read_with_extended_window_mode_rust, m)?)?;
    m.add_function(wrap_pyfunction!(run_one_window_rust, m)?)?;
    Ok(())
}
