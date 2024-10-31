#!/bin/sh
 
# edit this for an existing directory where you want to put the log files:
export LOGS_DIR='../tmp_log_dir'

cargo test log_serial_intrinsics_high_arith_int --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_moreflops_serialintr_logs.txt
cargo test log_serial_intrinsics_high_arith_int --features TIME_THR -- --test-threads 1 --nocapture >& $LOGS_DIR/thr_moreflops_serialintr_logs.txt
cargo test log_serial_intrinsics_high_arith_int --features TIME_D -- --test-threads 1 --nocapture >& $LOGS_DIR/d_moreflops_serialintr_logs.txt

cargo test log_serial_per_element_high_arith_int --features TIME_TOT -- --test-threads 1 --nocapture >& $LOGS_DIR/tot_moreflops_serial_per_logs.txt
