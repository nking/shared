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
        cargo test run_all_vec8prod_tests --features TIME_TOT -- --nocapture >& tot_logs.txt
        cargo test run_all_vec8prod_tests --features TIME_THR -- --nocapture >& thr_logs.txt
        cargo test run_all_vec8prod_tests --features TIME_D -- --nocapture >& d_logs.txt
           the time units are nano-seconds, so need to consider the best
           resolution of time for the hardware that the code is running 
           on (a clock cycle) when interpreting results.
           note that the c code logs for vec8prod use units clock cycles instead.

for the higher arithmetic intensity methods:
cargo test run_serial_per_element_high_arith_int --features TIME_THR -- --nocapture >& thr_per_logs.txt
cargo test run_serial_per_element_high_arith_int --features TIME_D -- --nocapture >& d_per_logs.txt
cargo test run_serial_intrinsics_high_arith_int --features TIME_THR -- --nocapture >& thr_intr_logs.txt
cargo test run_serial_intrinsics_high_arith_int --features TIME_D -- --nocapture >& d_intr_logs.txt

