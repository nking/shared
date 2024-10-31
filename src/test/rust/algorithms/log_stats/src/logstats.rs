use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use super::stats_holder::*;
use regex::Regex;

pub fn simple(file_path : String) {
    let mut thr_stats: StatsHolderMed = build_stats_holder_med(String::from("thr"));
    let mut d_stats: StatsHolderMed = build_stats_holder_med(String::from("d"));
    let mut tot_stats: StatsHolderMed = build_stats_holder_med(String::from("tot"));

    let re: Regex = Regex::new(r".+?\s(thr|tot|load)\s(\d+)$").unwrap();        

    if let Ok(lines) = read_lines(file_path) {
        // Consumes the iterator, returns an (Optional) String
        for line in lines.flatten() {
            
            if let Some(caps) = re.captures(&line) {
    
                let t = caps.get(2).unwrap().as_str().parse::<u64>().unwrap();
                match caps.get(1).unwrap().as_str() {
                    "thr" => {
                        thr_stats.accumulate(t);
                    },
                    "load" => {
                        d_stats.accumulate(t);
                    },
                    "tot" => {
                        tot_stats.accumulate(t);
                    },
                    _ => {}
                };
            } 
            //println!("{}", line);
        }
        if !thr_stats.is_empty() {
            //let (mn, stdev) = thr_stats.calc_mean_stdev();
            //print!("\nthr mean={:?}, stdev={:?}\n", mn, stdev);
            let med = thr_stats.calc_median();
            print!("\nthr median={:?}\n", med);
        }
        if !d_stats.is_empty() {
            let med = d_stats.calc_median();
            print!("\nd median={:?}\n", med);
        }
        if !tot_stats.is_empty() {
            let med = tot_stats.calc_median();
            print!("\ntot median={:?}\n", med);
        }
    }
}

fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
