The vec8prod algorithm refers to the ISPC algorithm from
Stanford CS 149, Parallel Programming
https://gfxcourses.stanford.edu/cs149
https://gfxcourses.stanford.edu/cs149/fall23/lecture/

The product goal is rewritten for different tools here 
to explore changes in use of memory and flops.
haven't started much of the analysis here...

CLI:
    to build:
        cargo build

    to run tests and get stdout:
        cargo test -- --nocapture

to create timing logs for analysis:
           the time units are nano-seconds, so need to consider the best
           resolution of time for the hardware that the code is running 
           on (a clock cycle) when interpreting results.
           note that the c code logs for vec8prod use units clock cycles instead.

for vec8prod:

you can edit create_logfiles.sh and use it as a shell script 

sh < create_logfiles.sh

it has these statements:

cargo test log_rand65536_serial --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_serial_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_D -- --test-threads 1 --nocapture >& $LOGS_DIR/d_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_D -- --test-threads 1 --nocapture >& $LOGS_DIR/d_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_no_vectorization --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_novec_logs.txt
cargo test log_rand65536_vec8prod_no_vectorization --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_novec_logs.txt

for the higher arithmetic intensity methods:

cargo test log_serial_per_element_high_arith_int --features TIME_THR -- --nocapture >& thr_per_logs.txt
cargo test log_serial_per_element_high_arith_int --features TIME_D -- --nocapture >& d_per_logs.txt
cargo test log_serial_per_element_high_arith_int --features TIME_TOT -- --nocapture >& tot_per_logs.txt
cargo test log_serial_intrinsics_high_arith_int --features TIME_THR -- --nocapture >& thr_intr_logs.txt
cargo test log_serial_intrinsics_high_arith_int --features TIME_D -- --nocapture >& d_intr_logs.txt
cargo test log_serial_intrinsics_high_arith_int --features TIME_TOT -- --nocapture >& tot_intr_logs.txt

