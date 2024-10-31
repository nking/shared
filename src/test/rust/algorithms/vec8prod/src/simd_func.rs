//use std::simd::f32x8;
//use std::simd::num::SimdFloat;
//use std::simd::StdFloat;

//use std::simd::*;
use std::simd::prelude::*;

//use std::simd::u32x8;
//use std::simd::i32x8;
//use std::simd::prelude::Simd;
use std::{panic, thread::scope, time::SystemTime};
use core::arch::x86_64::{__m256, _mm256_mul_ps, _mm256_storeu_ps, __m128, _mm_storeu_ps};
use core::arch::x86_64::{_mm256_loadu_ps, _mm_loadu_ps,  _mm256_permutevar8x32_ps, _mm256_set_epi32, _mm_fmaddsub_ps, _mm_rcp_ps, _mm_sqrt_ps};

use rand::Rng;
use rand::prelude::*;

//#[cfg(target_arch = "x86_64")]
//#[cfg(target_arch = "x86")]

// has ispc too:  https://crates.io/crates/packed_simd
// https://github.com/rust-lang/stdarch

/*pub fn sum_simd_5(a: &[f32], b: &[f32]) -> f32 {
    a.array_chunks::<8>()
        .map(|&a| f32x8::from_array(a))
        .zip(b.array_chunks::<8>().map(|&b| f32x8::from_array(b)))
        .fold(f32x8::splat(0.), |acc, (a, b)| a.mul_add(b, acc))
        .reduce_sum()
}*/

#[allow(non_snake_case)]
pub fn simd_func<const USE_SIMD: bool>(&N : &usize, x : & [f32]) -> f32 {
    //NOTE: we will have errors below in scope spawned thread copies of data
    // if parameter x is typed as x : &mut [f32] even if edited to mut below too

    // see https://doc.rust-lang.org/std/simd/index.html
    //println!("\nsimd module");

    if (N % 8) != 0 {
        panic!("N must be a multiple of 8");
    }

    #[allow(non_snake_case)]
    const N_VEC: usize = 8;

    let n_instances : usize = &N / N_VEC;

    //println!("N={}, N_VEC={}, n_instances={}", N, N_VEC, n_instances);

    let mut prod_results: Vec<f32> = Vec::new();

    scope(|s| {
        let mut split_prev : usize = 0;
        let mut split_next : usize = split_prev + N_VEC;

        for _i in 0..n_instances {
            
            let thr = s.spawn(move || {
                let mut xp:[f32; N_VEC] = [0.0f32; N_VEC];
                xp.copy_from_slice(& x[split_prev..split_next]);

                if USE_SIMD {
                    let res: f32 = simd_partition_thread_8(& mut xp);
                    res 
                } else {
                    let res: f32 = intrinsics_partition_thread(& mut xp);
                    res 
                }
            });
            // unwrap is used to get the result here.
            // unwrap: Extracts the value from an Option or Result type, panicking if the value is None or Err
            let res = thr.join().unwrap();
            //println!("join res={:#?}, and xp[0]={}", res, &xp[N_VEC-1]);
            prod_results.push(res);

            split_prev = split_next;
            split_next = split_prev + N_VEC;
        }
    });

    let mut res : f32 = 1.0f32;
    for r in &prod_results {
        res = res * r;
    }

    return res;
}

#[allow(non_snake_case)]

fn simd_partition_thread_8( x : & mut [f32; 8]) -> f32 {
    // TODO:
    // browse: https://doc.rust-lang.org/std/simd/prelude/trait.SimdFloat.html#tymethod.reduce_product
    
    #[cfg(feature = "TIME_THR")]
    let start = std::time::SystemTime::now();

    const N_VEC: usize = 8;

    #[cfg(feature = "TIME_D")]
    let start = std::time::SystemTime::now();

    let mut a_simd: Simd<f32, N_VEC> = Simd::from_array(*x);

    #[cfg(feature = "TIME_D")]
    match start.elapsed() {Ok(dur) => {tracing::info!("d load {:?}", dur.as_nanos())},Err(_) => {},}

    const SH1 : usize = 1;
    let mut b_simd = a_simd.rotate_elements_left::<SH1>();
    a_simd = a_simd * b_simd;

    const SH2 : usize = 2;
    b_simd = a_simd.rotate_elements_left::<SH2>();
    a_simd = a_simd * b_simd;
       
    const SH4 : usize = 4;
    b_simd = a_simd.rotate_elements_left::<SH4>();
    a_simd = a_simd * b_simd;

    #[cfg(feature = "TIME_D")]
    let start = std::time::SystemTime::now();

    //TODO: is there a faster way to extract only the first element?
    // ? something like _mm_cvtss_f32(_mm256_extractf128_ps(avx_x, 0))
    let r = a_simd.to_array()[0];

    #[cfg(feature = "TIME_D")]
    match start.elapsed() {Ok(dur) => {tracing::info!("d load {:?}", dur.as_nanos())},Err(_) => {},}

    //conditionally present 
    // e.g. cargo test --features TIME_THR -- --nocapture
    #[cfg(feature = "TIME_THR")]
    match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}

    return r;
}

fn intrinsics_partition_thread( x : & mut [f32]) -> f32 {
    
    if !is_x86_feature_detected!("avx2") {
        panic!("this method currently only supports avx2 and arch x86_64");
    }
    
    unsafe {
        
        #[warn(unused_variables)]
        #[cfg(feature = "TIME_THR")]
        let start = std::time::SystemTime::now();

        #[cfg(feature = "TIME_D")]
        let start = std::time::SystemTime::now();

        let mut avx_x: __m256 = _mm256_loadu_ps(x.get_unchecked(0));

        #[cfg(feature = "TIME_D")]
        let start = std::time::SystemTime::now();

        // shift right by 1
        let avx_y:__m256 = _mm256_permutevar8x32_ps(avx_x, 
            _mm256_set_epi32(0,7,6,5,4,3,2,1));
        
        #[cfg(feature = "TIME_D")]
        match start.elapsed() {Ok(dur) => {tracing::info!("shift load {:?}", dur.as_nanos())},Err(_) => {},}

        avx_x = _mm256_mul_ps(avx_x, avx_y);
        
        // shift by 2
        let avx_y:__m256 = _mm256_permutevar8x32_ps(avx_x, 
            _mm256_set_epi32(0,0, 7,6,5,4,3,2));
        
        avx_x = _mm256_mul_ps(avx_x, avx_y);

        // shift by 4
        let avx_y:__m256 = _mm256_permutevar8x32_ps(avx_x, 
            _mm256_set_epi32(0,0, 0, 0, 7,6,5,4));
            
        avx_x = _mm256_mul_ps(avx_x, avx_y);

        //extract first value
        // TODO: is there a faster way than storing register to memory and then
        // getting first element?

        let _ = _mm256_storeu_ps(x.as_mut_ptr(), avx_x);

        #[cfg(feature = "TIME_THR")]
        match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}

        return x[0];
    }
    
}

//============ triplets of methods to compare execution runtimes ======

fn _gen_random_float_array<const N : usize>() -> [f32; N] {
    
    let mut x = [0.0f32; N];

    let seed: u64 = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH)
    .unwrap().as_secs();
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
    println!("seed={:?}", seed);
    
    for xi in x.iter_mut() {
        *xi = 0.0001f32 + rng.gen::<f32>();
    }

    return x;
}

pub fn serial_per_element_high_arith_int<const N : usize>() -> () {
    
    let mut x: [f32; N] = _gen_random_float_array::<N>();

    _serial_per_element_high_arith_int::<N>(& mut x);
}

fn _serial_per_element_high_arith_int<const N : usize>( x : & mut [f32; N]) -> () {
    
    let factors : [f32; 4] = [3.0f32, 3.0f32, 3.0f32, 3.0f32];
    let s : [f32; 4] = [1.5f32, -1.5f32, 1.5f32, -1.5f32];

    #[cfg(feature = "TIME_THR")]
    let start = std::time::SystemTime::now();
    
    // specific to rust: can mutate a value while using the iterator to access it:
    //for xi in x.iter() {
    for j in 0..x.len() {
        for i in 0..4 {
            x[j] = ((x[j]) * factors[0]) + s[i%2];
            x[j] = 1./(x[j]);
            x[j] = (x[j]).sqrt();
        }
    }

    #[cfg(feature = "TIME_THR")]
    match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}
}

pub fn serial_intrinsics_high_arith_int<const N : usize>() -> () {
    let mut x: [f32; N] = _gen_random_float_array::<N>();
    _serial_intrinsics_high_arith_int::<N>(& mut x);
}

fn _serial_intrinsics_high_arith_int<const N : usize>( x : & mut [f32; N]) -> () {
    
    let factors : [f32; 4] = [3.0f32, 3.0f32, 3.0f32, 3.0f32];
    let s : [f32; 4] = [1.5f32, -1.5f32, 1.5f32, -1.5f32];

    unsafe {
        let avx_f: __m128 = _mm_loadu_ps(&factors[0]);
        let avx_s: __m128 = _mm_loadu_ps(&s[0]);

        #[cfg(feature = "TIME_THR")]
        let start = std::time::SystemTime::now();

        for _i in (0..N).step_by(4) {
                        
            #[cfg(feature = "TIME_D")]
            let start = std::time::SystemTime::now();

            //let mut avx_x: __m128 = _mm_loadu_ps(x.get_unchecked(_i));
            let mut avx_x: __m128 = _mm_loadu_ps(& mut x[_i]);
                
            #[cfg(feature = "TIME_D")]
            match start.elapsed() {Ok(dur) => {tracing::info!("load {:?}", dur.as_nanos())},Err(_) => {},}
        
            for _j in 0..4 {
                avx_x = _mm_fmaddsub_ps(avx_x, avx_f, avx_s);
                avx_x = _mm_rcp_ps(avx_x);
                avx_x = _mm_sqrt_ps(avx_x);
            }

            //#[cfg(feature = "TIME_D")]
            //let start = std::time::SystemTime::now();

            //* mut f32
            _mm_storeu_ps(& mut x[_i], avx_x);

            //#[cfg(feature = "TIME_D")]
            //match start.elapsed() {Ok(dur) => {tracing::info!("store {:?}", dur.as_nanos())},Err(_) => {},}
        }

        #[cfg(feature = "TIME_THR")]
        match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}
    }
}

/* 
fn serial_simd_high_arith_int<const N : usize>( x : & [f32]) -> () {
    let mut x: [f32; N] = _gen_random_float_array::<N>();
    _serial_simd_high_arith_int::<N>(&x);
}

fn _serial_simd_high_arith_int<const N : usize>( x : & mut [f32; N]) -> () {
    
    // for SIMD, need to copy x into n_VEC size arrays

    const N_VEC : usize = 4;
    let factors : [f32; 4] = [3.0f32, 3.0f32, 3.0f32, 3.0f32];
    let s : [f32; 4] = [1.5f32, -1.5f32, 1.5f32, -1.5f32];

    let mut split_prev : usize = 0;
    let mut split_next : usize = split_prev + N_VEC;
                            
    for _i in (0..N).step_by(N_VEC) {

        let mut xi:[f32; N_VEC] = [0.0f32; N_VEC];
        xi.copy_from_slice(& x[split_prev..split_next]);

        let mut a_simd: Simd<f32, N_VEC> = Simd::from_array(xi);

        no fmaddsub ...


        split_prev = split_next;
        split_next = split_prev + N_VEC;
    }

}
*/
