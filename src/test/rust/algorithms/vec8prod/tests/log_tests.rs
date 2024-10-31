#[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
#[cfg(test)]
mod log_tests {
    #[allow(unused_imports)]
    use vec8prod::simd_func::simd_func;
    #[allow(unused_imports)]
    use vec8prod::simd_func::serial_intrinsics_high_arith_int;
    #[allow(unused_imports)]
    use vec8prod::simd_func::serial_per_element_high_arith_int;
    #[allow(unused_imports)]
    use vec8prod::multithread_func::multithread_func;
    
    use rand::Rng;
    use rand::prelude::*;
    use std::{time::SystemTime};
    
    #[test]
    fn log_serial_intrinsics_high_arith_int() {
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_D"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_D", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        
        tracing_subscriber::fmt()
            .with_max_level(tracing::Level::DEBUG)
            // enable thread id to be emitted
            .with_thread_ids(true)
            // enabled thread name to be emitted
            .with_thread_names(true)
            .with_ansi(false)
            .init();
        // .with_ansi(false) avoids escape characters in output
    
        #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
        tracing::info!("serial_intrinsics");
    
        let n_tests = 3;
        for _ in 0..n_tests {
            serial_intrinsics_high_arith_int::<250000>();
        }
    }
    
    // only TIME_TOT feature is logged in this method
    #[test]
    fn log_rand65536_serial() {
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_D"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_D", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        
        tracing_subscriber::fmt()
            .with_max_level(tracing::Level::DEBUG)
            // enable thread id to be emitted
            .with_thread_ids(true)
            // enabled thread name to be emitted
            .with_thread_names(true)
            .with_ansi(false)
            .compact()
            .init();
        // .with_ansi(false) avoids escape characters in output
    
        #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
        tracing::info!("serial_per_element");
    
        let n_tests = 3;
        for _ in 0..n_tests {
            let x = generate_rand_x::<65536>();
    
            #[allow(unused_variables)]
            #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
            let start = std::time::SystemTime::now();
    
            let mut exp_ans = 1.0f32;
            for xi in x.iter() {
                exp_ans = exp_ans * *xi;
            }
    
            #[cfg(feature = "TIME_TOT")]
            match start.elapsed() {Ok(dur) => {tracing::info!("tot {:?}", dur.as_nanos())},Err(_) => {},}
        }   
    }
    
    #[test]
    fn log_rand65536_vec8prod_intrinsics() {
        _log_rand_vec8prod::<65536>("multithreadintrinsics", "intr");
    }
    #[test]
    fn log_rand65536_vec8prod_simd() {
        _log_rand_vec8prod::<65536>("multithreadsimd", "simd");
    }
    #[test]
    fn log_rand65536_vec8prod_no_vectorization() {
        _log_rand_vec8prod::<65536>("multithreadnovec", "novec");
    }
    
    #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
    fn _log_rand_vec8prod<const N: usize>(label: &str, which_test : &str) {
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_D"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_D", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
    
        //RUST_LOG="vec8prod"=trace
        tracing_subscriber::fmt()
            .with_max_level(tracing::Level::DEBUG)
            // enable thread id to be emitted
            .with_thread_ids(true)
            // enabled thread name to be emitted
            .with_thread_names(false)
            .with_ansi(false)
            .compact()
            .init();
        // .with_ansi(false) avoids escape characters in output
    
        #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
        tracing::info!(label);
    
        let n_tests = 3;
        for _ in 0..n_tests {
            let x = generate_rand_x::<N>();
    
            #[allow(unused_variables)]
            #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
            let start = std::time::SystemTime::now();
    
            let mut ans:f32 = 1.;
            if which_test.eq("simd") {
                ans = simd_func::<true>(& N, & x);
            } else  if which_test.eq("intr") {
                ans = simd_func::<false>(& N, & x);
            } else if which_test.eq("novec") {
                ans = multithread_func(&N, &x);
            }
    
            #[cfg(feature = "TIME_TOT")]
            match start.elapsed() {Ok(dur) => {tracing::info!("tot {:?}", dur.as_nanos())},Err(_) => {},}
        
            let mut exp_ans = 1.0f32;
            for xi in x.iter() {
                exp_ans = exp_ans * *xi;
            }
            //println!("main.  exp_ans={}, ans1={}", exp_ans, ans1);
            let r : f32 = ((exp_ans/ans) - 1.0).abs();
            assert!( r < 5E-5);
        }   
    }
    
    #[test]
    fn log_serial_per_element_high_arith_int() {
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_TOT", feature = "TIME_D"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        #[cfg(all(feature = "TIME_D", feature = "TIME_THR"))]
        panic!("Cannot have more than 1 of the logging features set (TIME_TOT, TIME_THR, TIME_D).");
        
        tracing_subscriber::fmt()
            .with_max_level(tracing::Level::DEBUG)
            // enable thread id to be emitted
            .with_thread_ids(true)
            // enabled thread name to be emitted
            .with_thread_names(true)
            .with_ansi(false)
            .init();
        // .with_ansi(false) avoids escape characters in output
    
        #[cfg(any(feature = "TIME_TOT", feature = "TIME_THR", feature = "TIME_D"))]
        tracing::info!("serial_per_element");
    
        let n_tests = 3;
        for _ in 0..n_tests {
            serial_per_element_high_arith_int::<250000>();
        }   
    }
    
    #[allow(dead_code)]
    fn generate_rand_x<const N : usize>() -> [f32; N] {
        let mut x = [0.0f32; N];
        
        let seed: u64 = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH)
        .unwrap().as_secs();
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
        println!("seed={:?}", seed);
    
        let factor = 1.0f32  - f32::MAX.powf(1.0f32/N as f32);
    
        for xi in x.iter_mut() {
            *xi = 1.0f32 + factor * rng.gen::<f32>();
        }
    
        return x;
    }
    
}
