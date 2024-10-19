
use vec8prod::simd_func::simd_func;
use vec8prod::multithread_func::multithread_func;
use rand::Rng;
use rand::prelude::*;
use std::time::SystemTime;

#[test]
fn run_all_tests() {

    tracing_subscriber::fmt()
        .with_max_level(tracing::Level::DEBUG)
        // enable thread id to be emitted
        .with_thread_ids(true)
        // enabled thread name to be emitted
        .with_thread_names(true)
        .init();

    test_16();
    test_rand_128();
}

fn test<const N : usize>( x : & [f32]) -> () {
    
    //INIT_TIME();
    //START_TOT_TIME();

    let mut exp_ans = 1.0f32;
    for xi in x.iter() {
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

    let ans2 = simd_func::<false>(& N, & x2);
    r  = ((exp_ans/ans2) - 1.0).abs();
    assert!( r < 5E-5);

    let ans3 = simd_func::<true>(& N, & x2);
    r = ((exp_ans/ans3) - 1.0).abs();
    assert!( r < 5E-5);

}


fn test_16() -> () {

    //INIT_TIME_TITLE(rusttest16);

    const N: usize = 16;
    let mut x = [0.0f32; N];

    for i in 0..N {
        x[i] = (i + 10) as f32;
    }
    
    test::<N>(&x);
}

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

    test::<N>(&x);
    
}


