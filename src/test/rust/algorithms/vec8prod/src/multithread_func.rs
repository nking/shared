use std::{panic, thread};

#[allow(non_snake_case)]
pub fn multithread_partition_func(&N : &usize, x : &[f32]) -> f32 {
    
    #[allow(non_snake_case)]
    const N_VEC: usize = 8;

    let n_instances : usize = &N / N_VEC;

    let mut prod_results: Vec<f32> = Vec::new();
   
    let mut split_prev : usize = 0;
    let mut split_next : usize = split_prev + N_VEC;

    for _i in 0..n_instances {
                    
        let mut xp:[f32; N_VEC] = [0.0f32; N_VEC];
        xp.copy_from_slice(&x[split_prev..split_next]);

        let thr = thread::spawn(move || {
            let res: f32 = multithread_partition_thread(& N_VEC, &mut xp); 
            res 
        });
        // unwrap is needed to get the result
        let res = thr.join().unwrap();
        //println!("join res={:#?}, and xp[0]={}", res, &xp[N_VEC-1]);
        prod_results.push(res);

        split_prev = split_next;
        split_next = split_prev + N_VEC;
    }

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

    // rust is designed to prevent the unguarded concurrent access of the
    // mutable data (x), so even though I designed the algorithm to operate on 
    // independent sections of x, rust cannot let the
    // threads share the array x without guarded synchronization.
    //
    // I can modify the algorithm in 4 ways:
    // (1) copy partition the mutable array x into n_instances.
    // let each thread process its copied partition and return a result.
    // (2) wrap the entire data array x, unpartitioned into an Arc 
    // (e.g. Arc<Mutex<x>> and let rust
    //     perform synchronization (hopefully, only where needed, which in
    //     these independent partitions will be needed nowhere)).
    // (2) rayon
    // (3) crossbeam
    
}

#[allow(non_snake_case)]
fn multithread_partition_thread(&N_VEC : &usize, x : &mut [f32]) -> f32 {
    
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
