#[cfg(test)]
pub mod tests {

    use regex::Regex;

    #[test]
    pub fn test_re0() {
    
        let re: Regex = Regex::new(r".+\s(thr|tot|d|load)\s(\d+)$").unwrap(); 
    
        let line = String::from("2024-10-30T00:52:14.807749Z  INFO \
                            run_serial_intrinsics_high_arith_int ThreadId(02) \
                            vec8prod::simd_func: thr 17132000");

        assert!(_test_re(line, &re));

        let line = String::from("2024-10-30T00:52:16.024670Z  INFO \
                                run_serial_intrinsics_high_arith_int ThreadId(02) \
                                vec8prod::simd_func: load 0");

        assert!(_test_re(line, &re));

    }

    fn _test_re(line : String, re: & Regex) -> bool {
        if let Some(caps) = re.captures(&line) {
            let log_type = caps.get(1).unwrap().as_str();
            let t = caps.get(2).unwrap().as_str().parse::<u64>().unwrap();
            print!("log_type={log_type}\ntime={t}\n");
            return true;
        } else {
            print!("ERROR\n");
            return false;
        }
    }
}
