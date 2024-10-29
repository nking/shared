use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use super::stats_holder::*;
use regex::Regex;

pub fn simple(file_path : String) {
    let mut thr_stats = build_StatsHolder(String::from("thr"));
    let mut d_stats = build_StatsHolder(String::from("d"));
    let mut tot_stats = build_StatsHolder(String::from("tot"));

    let thrType = String::from("thr");
    let totType = String::from("tot");
    let dType = String::from("d");

    let re: Regex = Regex::new(r".+?\s(thr|tot|d)\s(\d+)$").unwrap();        

    if let Ok(lines) = read_lines(file_path) {
        // Consumes the iterator, returns an (Optional) String
        
        for line in lines.flatten() {
            //TODO: is there a more efficient way?
            
            if let Some(caps) = re.captures(&line) {
                let log_type = caps.get(1).unwrap().as_str();
                let t = caps.get(2).unwrap().as_str().parse::<u64>().unwrap();
                if log_type.eq(&thrType) {
                    thr_stats.accumulate(t);
                } else if log_type.eq(&totType) {
                    tot_stats.accumulate(t);
                } else if log_type.eq(&dType) {
                    d_stats.accumulate(t);
                }
            } 
            //println!("{}", line);
        }
        if !thr_stats.isEmpty() {
            let (mn, stdev) = thr_stats.calc_mean_stdev();
            print!("\nthr mean={:?}, stdev={:?}\n", mn, stdev);
        }
        if !d_stats.isEmpty() {
            let (mn, stdev) = d_stats.calc_mean_stdev();
            print!("\nd mean={:?}, stdev={:?}\n", mn, stdev);
        }
        if !tot_stats.isEmpty() {
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
