use std::{panic, thread::scope};

//#![feature(portable_simd)]
// has ispc too:  https://crates.io/crates/packed_simd
// https://github.com/rust-lang/stdarch

//use std::simd;
// to compare to the c code, use the 32 bit implementations:
//use std::{f32, mem};
#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
#[cfg(target_arch = "x86")]
use std::arch::x86::*;

#[allow(non_snake_case)]
pub fn simd_func(&N : &usize, x : &[f32]) -> f32 {
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
                xp.copy_from_slice(&x[split_prev..split_next]);

                let res: f32 = intrinsics_partition_thread(& mut xp);
                res 
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
fn simd_partition_thread( x : &mut [f32]) -> f32 {
    // TODO: implement
    // browse: https://towardsdatascience.com/nine-rules-for-simd-acceleration-of-your-rust-code-part-1-c16fe639ce21
    // browse:  https://monadera.com/blog/faster-rust-with-simd/
    return 0.0f32;
}

#[allow(non_snake_case)]
fn intrinsics_partition_thread( x : &mut [f32]) -> f32 {
    
    unsafe {
    //core::arch::x86_64::
    let mut avx_x: __m256 = _mm256_loadu_ps(x.get_unchecked(0));

    let mut n_iter : usize = 0;
    while n_iter < 3 {

        let shift : usize = 1 << n_iter;

        // avx_y needs to be avx_x shifted to the left by shift indices
        // TODO: there should be an efficient way to do this with intrinsics.
        //     need this offset to left by shift indexes:
        //         let mut avx_y : __m256 = avx_x.clone();
        // TODO: consider __m256 _mm256_mask_i32gather_ps
        //https://doc.rust-lang.org/beta/core/arch/x86_64/fn._mm256_mask_i32gather_ps.html
        //
        // for now, will store to x and source avx_y from it with slice
        let _ = _mm256_storeu_ps(x.as_mut_ptr(), avx_x);
        let avx_y: __m256 = _mm256_loadu_ps(x.get_unchecked(shift));
        
        avx_x = _mm256_mul_ps(avx_x, avx_y);

        n_iter += 1;
    }

    let _ = _mm256_storeu_ps(x.as_mut_ptr(), avx_x);
    //core_arch::simd::f32x8
    //avx_x.as_f32x8
    // under the hood: pub struct __m256(f32, f32, f32, f32, f32, f32, f32, f32);
    //println!("   in simd, res avx_x[0]={:#?}", x[0]);

    return x[0];
}
}