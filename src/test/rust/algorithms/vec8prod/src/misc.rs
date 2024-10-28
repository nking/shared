use std::array;

/**
 * given size 16 array, compute the product.
 * this is not efficient, just a look at the pointer math of borrowing
 * and ownership and partial copies of arrays
 */
fn _not_for_use_prod16(a: &[f32; 16]) -> f32 {
    const N_VEC: usize = 16;
    let mut n_iter : usize = 0;
    let &(mut prod) = a;
    while n_iter < 3 {
        let shift : usize = 1 << n_iter;

        let mut b:[f32; N_VEC] = [0.0f32; N_VEC];
        b[0 .. N_VEC - shift].copy_from_slice(&prod[shift..N_VEC]);

        //println!("b={:#?}", b);

        // TODO: replace with simd vector math
        let prod2 = array::from_fn::<f32, N_VEC, _>(|i| prod[i] * b[i]);
        
        prod = prod2;
        //println!("prod={:#?}", prod);

        n_iter += 1;
    }
    
    return prod[0] * prod[8];
}
