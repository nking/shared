The vec8prod algorithm refers to the ISPC algorithm from
Stanford CS 149, Parallel Programming
https://gfxcourses.stanford.edu/cs149
https://gfxcourses.stanford.edu/cs149/fall23/lecture/

This project's goal is to use different languages and tools to compare
algorithms that calculate the product of an array where the
different algorithms aren't an exhaustive set, but are for comparison
for performance and efficiency statistics.

The project is tucked away into the test branch because it's exploratory
rather than meant to be used as part of the API (...might consider moving
it into its own project for clarity).

The vec8prod algorithm is a divide and conquer design that the Stanford
CS149 authors wrote for ISPC.

see notes.txt here and in sub-directories for analysis.


vec8prod and serial per element product, and serial per intrinsics vector product,
and parallel intrinsics vector product
have compile-time configurable time logging.

   to build and run, cd to the base directory of shared.
     for the ispc divide and conquer code (vec8prod):
        ant cmakeISPCTimeLogs
     for the parallel divide and conquer intrinsics code:
      ant cmakeSIMDTimeLogs
     for the multithread divide and conquer:
      ant cmakeSimSIMDTimeLogs

Then to look at timing logs for a higher arithmetic intensity algorithm
use among some of the tools, *_explore_* methods were made.

The brief analyses of the timing logs are in the file notes.txt.

A similar but smaller analysis is in progress in the rust directory
shared/src/test/rust/algorithms/vec8prod

Regarding analyzing the logs from the c directory:

    to analyze the data logs after cleaning for titles etc,
    can use
    spreadsheet apps (Excel, Numbers, etc.)
    or posix command-line system programs or write code.

    tinkering with awk is fast an easy, and Gemini can tell you how to
    perform something with awk or posix command line etc.
    though it gets messy past simple goals and python or rust etc.
    might be better for not repeating the work later...

    I don't know the comparable to awk for windows operating systems shells
    (and use of Microsoft copilot or an apple AI to explore asking shell questions) currently
    <add info here>...

       This example is for the _d_ logs which contain only the data load/store
       times and so don't need a filter for title sections.

       e.g. to group the data load/store times and count them:
         from a file made from running explore_cache_timing() in the simd tests:

         awk '{print $3}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt | sort -g | uniq -c

       e.g. to get the average:
         awk '{sum += $3} END {avg = sum / NR; print "avg: " avg}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt
           prints avg: 0.963592

       having the avg now:
         awk '{sum += ($3 - 0.964)^2} END {st = sqrt(sum/(NR-1)); print "stdv: " st}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt
          prints stdv: 0.836941

         or combine avg and std into 1 messy script

       to get the median (without correction for even number of numbers):
         awk '{print $3}' bin/test-classes-c/algorithms/vec8prod/ispc8_d_logs.txt | sort -n | awk '{a[NR]=$1} END {print "median: " a[NR/2]}'


       then look at the i values where load times are larger than 3 * stdv above avg = 2.5 + 1 = 3.5
        or avg + 10 * std ~ 11
           awk -F ' '  '($3 > 11) {print $1 , $3}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt

        or to look at the differences in loop index i consecutively for those filtered values,
        pipe a modification of the above into awk again:
           awk 'BEGIN {prev=0}{print $1 - prev; prev = $1}'

           awk -F ' '  '($3 > 11) {print $1}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt | awk 'BEGIN {prev=0}{print $1 - prev; prev = $1}'


       and a messy get average for those 10 sigma events:
              awk -F ' '  '($3 > (0.95+10*0.77)) {print $1}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt | awk 'BEGIN {prev=0}{print $1 - prev; prev = $1}' | awk '{sum += $1} END {avg = sum / NR; print "avg: " avg}'
              prints avg: 1536.03

       then take whatever the average is and replace 1536 in here with it:
          awk -F ' '  '($3 > (0.95+10*0.77)) {print $1}' bin/test-classes-c/algorithms/vec8prod/simd/simd_d_logs.txt | awk 'BEGIN {prev=0}{print $1 - prev; prev = $1}' | awk '{sum += ($1 - 1536)^2} END {st = sqrt(sum/(NR-1)); print "stdv: " st}'
          prints 2045

