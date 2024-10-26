The vec8prod algorithm refers to the ISPC algorithm from
Stanford CS 149, Parallel Programming
(look for current year's ISPC lecture on their website:
https://gfxcourses.stanford.edu/cs149)
The product goal is rewritten for different tools here 
to explore changes in use of memory and flops.

there are ant targets to run the codes, and I haven't documented them
yet.

to analyze the data logs after cleaning for titles etc,
can use 
spreadsheet apps (Excel, Numbers, etc.)
or posix command-line system programs or write code.

tinkering with awk is fast an easy, and Gemini can tell you how to
perform something with awk or posix command line etc.
though it gets messy past simple goals and python or rust etc
might be better for not repeating the work later...

I don't know the comparable to awk for windows operating systems shells
(and use of Microsoft copilot to explore asking questions) currently
<add info here>...

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

