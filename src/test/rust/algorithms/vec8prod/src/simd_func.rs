//#![feature(portable_simd)]
// has ispc too:  https://crates.io/crates/packed_simd
// https://github.com/rust-lang/stdarch
//use std::simd;

#[allow(non_snake_case)]
pub fn simd_func(&N : &usize, x : &mut [f32]) -> () {
    // see https://doc.rust-lang.org/std/simd/index.html
    println!("simd module");
}
