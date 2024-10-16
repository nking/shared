use std::{panic, thread::scope};

/*
rust is designed to prevent the unguarded concurrent access of the
mutable data (x), so even though I designed the c implementation
of the vec8prod algorithm to operate on 
independent sections of shared array x, rust cannot let the
threads share the array x without guarded synchronization.
In the below loop for _i in 0..n_instances,
there can only be one spawned thread which has a mutable borrow
of x.

I can modify the algorithm for rust in 4 ways:
(1) copy partition the mutable array x into n_instances.
let each thread process its copied partition and return a result.
(2) wrap the entire data array x, unpartitioned into an Arc 
(e.g. Arc<Mutex<x>>. The Arc mutex is cloned (each having same 
reference to the shared mutable data array x), and given to each thread.
A thread accesses the shared mutable data array x by obtaining an exclusive
lock on the data to read and or modify it, then releases the
lock, then another thread can obtain the mutex lock, ...
(2) rayon
(3) crossbeam  
*/

#[allow(non_snake_case)]
pub fn multithread_partition_func(&N : &usize, x : &[f32]) -> f32 {
    //NOTE: we have errors below in scope spawned thread copies of data
    // if parameter x is typed as x : &mut [f32] even if edited to mut below too

    #[allow(non_snake_case)]
    const N_VEC: usize = 8;

    let n_instances : usize = &N / N_VEC;

    let mut prod_results: Vec<f32> = Vec::new();
   
    scope(|s| {
        let mut split_prev : usize = 0;
        let mut split_next : usize = split_prev + N_VEC;

        for _i in 0..n_instances {

            //NOTE, if xp construction and copy from slice is
            // placed here instead, there are 2 copies of N_VEC wide
            // arrays for each _i (one outside of thread spawn and
            // the other inside of thread spawn).
            // Instead, placing xp construction and copy from slice
            // inside of spawned thread results in only 1 copy being made.
            // For this, also need scope to share x
            
            let thr = s.spawn(move || {
                let mut xp:[f32; N_VEC] = [0.0f32; N_VEC];
                xp.copy_from_slice(&x[split_prev..split_next]);

                //println!("i={}, x address={:p}, xp address={:p}", &_i, &x[0], &xp[0]);

                let res: f32 = multithread_partition_thread(& N_VEC, &mut xp); 
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
pub fn multithread_func(&N : &usize, x : &mut [f32]) -> f32 {

    if (N % 8) != 0 {
        panic!("N must be a multiple of 8");
    }
    
    //println!("multithread module");
    //println!("multithread module. N={}, x[0]={}, address={:p}", &N, &x[0], &x[0]);
    // see https://doc.rust-lang.org/std/thread/index.html

    //let n_instances : usize = &N / 8;

    // cannot pass x again as mutable.
    // if pass as not mutable, it needs to return a value
    let r1 = multithread_partition_func(&N, &x);

    return r1;

}

#[allow(non_snake_case)]
fn multithread_partition_thread(&N_VEC : &usize, x : &mut [f32]) -> f32 {
    //println!("   x address={:p}", &x[0]);
    let mut p_id : usize = 2;
    while p_id <= N_VEC {
        let off0: usize = p_id - 1;
        let off1: usize = (p_id/2) - 1;

        let mut j :usize = 0;
        while j < N_VEC {
        
            x[j + off0] = &x[j + off0] * &x[j+off1];

            j += p_id;
        }
        
        p_id = p_id << 1;
    }

    return x[N_VEC - 1];
}
