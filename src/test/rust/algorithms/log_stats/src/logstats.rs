use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use super::stats_holder::*;
use regex::Regex;

pub fn simple(file_path : String) {
    let mut thr_stats = build_stats_holder(String::from("thr"));
    let mut d_stats = build_stats_holder(String::from("d"));
    let mut tot_stats = build_stats_holder(String::from("tot"));

    let thr_type = String::from("thr");
    let tot_type = String::from("tot");
    let d_type = String::from("load");

    let re: Regex = Regex::new(r".+?\s(thr|tot|load)\s(\d+)$").unwrap();        

    let mut count = 0;

    if let Ok(lines) = read_lines(file_path) {
        // Consumes the iterator, returns an (Optional) String
        
        for line in lines.flatten() {
            //TODO: is there a more efficient way?
            
            if let Some(caps) = re.captures(&line) {
                let log_type = caps.get(1).unwrap().as_str();
                let t = caps.get(2).unwrap().as_str().parse::<u64>().unwrap();
                if log_type.eq(&thr_type) {
                    thr_stats.accumulate(t);
                    count += 1;
                } else if log_type.eq(&tot_type) {
                    tot_stats.accumulate(t);
                    count += 1;
                } else if log_type.eq(&d_type) {
                    d_stats.accumulate(t);
                    count += 1;
                }
            } 
            //println!("{}", line);
        }
        print!("count={:?}", count);
        if !thr_stats.is_empty() {
            let (mn, stdev) = thr_stats.calc_mean_stdev();
            print!("\nthr mean={:?}, stdev={:?}\n", mn, stdev);
        }
        if !d_stats.is_empty() {
            let (mn, stdev) = d_stats.calc_mean_stdev();
            print!("\nd mean={:?}, stdev={:?}\n", mn, stdev);
        }
        if !tot_stats.is_empty() {
            let (mn, stdev) = tot_stats.calc_mean_stdev();
            print!("\ntot mean={:?}, stdev={:?}\n", mn, stdev);
        }
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
