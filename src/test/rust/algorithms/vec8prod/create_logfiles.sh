#!/bin/sh
 
# edit this for an existing directory where you want to put the log files:
export LOGS_DIR='../tmp_log_dir'

cargo test log_rand65536_serial --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_serial_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_intrinsics --features TIME_D -- --test-threads 1 --nocapture >& $LOGS_DIR/d_vecprod_intr_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_simd --features TIME_D -- --test-threads 1 --nocapture >& $LOGS_DIR/d_vecprod_simd_logs.txt
cargo test log_rand65536_vec8prod_no_vectorization --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_vecprod_novec_logs.txt
cargo test log_rand65536_vec8prod_no_vectorization --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_vecprod_novec_logs.txt

