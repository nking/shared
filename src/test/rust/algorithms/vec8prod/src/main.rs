#![feature(portable_simd)]
#![feature(array_chunks)]

pub mod simd_func;
pub mod multithread_func;
mod misc;

#[allow(unused_imports)]
use crate::simd_func::{simd_func};
#[allow(unused_imports)]
use crate::multithread_func::multithread_func;
#[allow(unused_imports)]
use rand::Rng;
#[allow(unused_imports)]
use rand::prelude::*;
#[allow(unused_imports)]
use rand_chacha::ChaCha8Rng;
#[allow(unused_imports)]
use std::time::SystemTime;

//Packages: A Cargo feature that lets you build, test, and share crates
//  Crates: A tree of modules that produces a library or executable
//  Modules and use: Let you control the organization, scope, and privacy of paths
//  Paths: A way of naming an item, such as a struct, function, or module

// crates can be a binary or library type.
//   binary: can compile to executable.  must have a main function.
//           src/bin/
//   library: do not have a main function.   they define functionality to be
//            used by projects.
// crate root is a source file and is the root module of the crate.
//    e.g. main.rs
// a package is a bundle of 1 or more crates that provides a set of functionality.
//   it contains a Cargo.toml file.
//   can contain 1 or more binary crates, but only 1 library crate.
//   package directory has its own src/lib.rs

fn test<const N : usize>( x : &mut [f32]) -> () {
    
    //INIT_TIME();
    //START_TOT_TIME();

    let mut exp_ans = 1.0f32;
    for xi in x.iter_mut() {
        exp_ans = exp_ans * *xi;
    }

    //STOP_TOT_TIME(serial);

    let mut x2 : [f32; N] = [0.0f32; N];
    let mut x3: [f32; N] = [0.0f32; N];
    x2.copy_from_slice(& x[0..N]);
    x3.copy_from_slice(& x[0..N]);

    let ans1 = multithread_func(&N, x);

    //println!("main.  exp_ans={}, ans1={}", exp_ans, ans1);
    
    let mut r : f32 = ((exp_ans/ans1) - 1.0).abs();
    assert!( r < 5E-5);

    let ans2 = simd_func::<false>(& N, &mut x2);
    r  = ((exp_ans/ans2) - 1.0).abs();
    assert!( r < 5E-5);

    let ans3 = simd_func::<true>(& N, &mut x2);
    r = ((exp_ans/ans3) - 1.0).abs();
    assert!( r < 5E-5);

}

#[test]
fn test_16() -> () {

    //INIT_TIME_TITLE(rusttest16);

    const N: usize = 16;
    let mut x = [0.0f32; N];

    for i in 0..N {
        x[i] = (i + 10) as f32;
    }
    
    test::<N>(&mut x);
}

#[test]
fn test_rand_128() -> () {
    
    //INIT_TIME_TITLE(rusttest128);

    // generate 128 random vector of numbers in range [1,1.65] whose product is <= 3.4E38
    const N: usize = 128;
    let mut x = [0.0f32; N];

    let seed: u64 = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH)
    .unwrap().as_secs();
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(seed);
    println!("seed={:?}", seed);
    
    for xi in x.iter_mut() {
        *xi = 0.65f32 + rng.gen::<f32>();
    }

    test::<N>(&mut x);
     
}

fn main() {

}


