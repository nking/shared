The vec8prod algorithm refers to the ISPC algorithm from
Stanford CS 149, Parallel Programming
(look for current year's ISPC lecture on their website:
https://gfxcourses.stanford.edu/cs149)
The product goal is rewritten for different tools here 
to explore changes in use of memory and flops.
haven't started much of the analysis here...

CLI:
    to build:
        cargo build

    to run tests and get stdout:
        cargo test -- --nocapture

    to create timing logs for analysis:
        cargo test --features TIME_TOT -- --nocapture >& tot_logs.txt
        cargo test --features TIME_THR -- --nocapture >& thr_logs.txt
        cargo test --features TIME_D -- --nocapture >& d_logs.txt
           the time units are nano-seconds, so need to consider the best
           resolution of time for the hardware that the code is running 
           on (a clock cycle) when interpreting results.
           note that the c code logs for vec8prod use units clock cycles instead.
