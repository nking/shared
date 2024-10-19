//use std::simd::f32x8;
//use std::simd::num::SimdFloat;
//use std::simd::StdFloat;

//use std::simd::*;
use std::simd::prelude::*;

//use std::simd::u32x8;
//use std::simd::i32x8;
//use std::simd::prelude::Simd;
use std::{panic, thread::scope};
use std::arch::x86_64::{_mm256_loadu_ps, _mm256_mul_ps, _mm256_storeu_ps, __m256};

//#[cfg(target_arch = "x86_64")]
//#[cfg(target_arch = "x86")]

// has ispc too:  https://crates.io/crates/packed_simd
// https://github.com/rust-lang/stdarch

//use std::arch::x86_64::*;
//use std::arch::x86_64::{__m256i};
//use std::arch::x86::*;


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
    match start.elapsed() {Ok(dur) => {tracing::info!("D LOAD {:?}", dur.as_nanos())},Err(_) => {},}

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
    let r = a_simd.to_array()[0];

    #[cfg(feature = "TIME_D")]
    match start.elapsed() {Ok(dur) => {tracing::info!("D LOAD {:?}", dur.as_nanos())},Err(_) => {},}

    //conditionally present 
    // e.g. cargo test --features TIME_THR -- --nocapture
    #[cfg(feature = "TIME_THR")]
    match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}

    return r;
}

/*const unsafe fn create_shift_1() -> __m256i /*i32x8*/ {
    unsafe { _mm256_set_epi32(1,2,3,4,5,6,7,0) } 
}*/


#[allow(non_snake_case)]
fn intrinsics_partition_thread( x : & mut [f32]) -> f32 {
    
    unsafe {
        
        #[cfg(feature = "TIME_THR")]
        let start = std::time::SystemTime::now();

        #[cfg(feature = "TIME_D")]
        let start = std::time::SystemTime::now();

        //core::arch::x86_64::
        let mut avx_x: __m256 = _mm256_loadu_ps(x.get_unchecked(0));

        #[cfg(feature = "TIME_D")]
        match start.elapsed() {Ok(dur) => {tracing::info!("load {:?}", dur.as_nanos())},Err(_) => {},}

        let mut n_iter : usize = 0;
        while n_iter < 3 {

            let shift : usize = 1 << n_iter;

            // avx_y needs to be avx_x shifted to the left by shift indices
            // TODO: there should be an efficient way to do this with intrinsics.
            //     need this offset to left by shift indexes:
            //         let mut avx_y : __m256 = avx_x.clone();
            //            then in-place shuffle: shufps
            //            or vpermps
            //            no success with _mm256_permute_ps
            // TODO: consider __m256 _mm256_mask_i32gather_ps
            //https://doc.rust-lang.org/beta/core/arch/x86_64/fn._mm256_mask_i32gather_ps.html
            //
            // for now, will store to x and source avx_y from it with slice
            
            #[cfg(feature = "TIME_D")]
            let start = std::time::SystemTime::now();

            let _ = _mm256_storeu_ps(x.as_mut_ptr(), avx_x);
            
            #[cfg(feature = "TIME_D")]
            match start.elapsed() {Ok(dur) => {tracing::info!("store {:?}", dur.as_nanos())},Err(_) => {},}

            #[cfg(feature = "TIME_D")]
            let start = std::time::SystemTime::now();

            let avx_y: __m256 = _mm256_loadu_ps(x.get_unchecked(shift));
            
            #[cfg(feature = "TIME_D")]
            match start.elapsed() {Ok(dur) => {tracing::info!("load {:?}", dur.as_nanos())},Err(_) => {},}

            avx_x = _mm256_mul_ps(avx_x, avx_y);

            n_iter += 1;
        }

        #[cfg(feature = "TIME_D")]
        let start = std::time::SystemTime::now();

        let _ = _mm256_storeu_ps(x.as_mut_ptr(), avx_x);

        #[cfg(feature = "TIME_D")]
        match start.elapsed() {Ok(dur) => {tracing::info!("load {:?}", dur.as_nanos())},Err(_) => {},}

        //core_arch::simd::f32x8
        //avx_x.as_f32x8
        // under the hood: pub struct __m256(f32, f32, f32, f32, f32, f32, f32, f32);
        //println!("   in simd, res avx_x[0]={:#?}", x[0]);

        #[cfg(feature = "TIME_THR")]
        match start.elapsed() {Ok(dur) => {tracing::info!("thr {:?}", dur.as_nanos())},Err(_) => {},}

        return x[0];
    }
    
}