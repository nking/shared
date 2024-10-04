package algorithms.scheduling;

import algorithms.sort.MiscSorter;
import gnu.trove.iterator.TIntDoubleIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;

import java.util.*;
import java.util.stream.IntStream;

/**
 *
 * @author nichole
 */
public class Misc {
        
     /**
     * For a single resource, schedule a set of n tasks where each task is associated with a execution time
     * t_i and a deadline d_i. 
     * The objective is to schedule the tasks, no two overlapping in time, 
     * such that they are all completed before their deadline. 
     * If this is not possible, define the lateness of the ith task to be the amount
     * by which its finish time exceeds its deadline. 
     * The objective is to minimize the maximum lateness over all the tasks.
     * 
     * References:
     * <pre>
     * lecture 7 notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect07-greedy-sched.pdf
     * </pre>
     * 
     * runtime complexity O(N * log_2(N)).
      * The results are optimal in minimizing the maximum lateness.
     * 
     @param duration duration of task. duration a.k.a. burst time a.k.a. execution time
     @param deadline deadline for task
     @param outputStart output array to hold start times for the resulting scheduled index order 
     @param outputLate output array to hold lateness for the resulting scheduled index order.
     * if is on time, element will be 0.
     @return indexes for scheduling order
     */
    public int[] unweightedIntervalMinimizeLateGreedy(double[] duration, double[] deadline,
        double[] outputStart, double[] outputLate) {
        int n = duration.length;
        if (deadline.length != n) {
            throw new IllegalArgumentException("d.length must equal t.length");
        }
        if (outputStart.length != n) {
            throw new IllegalArgumentException("outputStart.length must equal t.length");
        }
        if (outputLate.length != n) {
            throw new IllegalArgumentException("outputLate.length must equal t.length");
        }
        //duration = Arrays.copyOf(duration, duration.length);
        //deadline = Arrays.copyOf(deadline, deadline.length);
        Arrays.fill(outputStart, 0);
        Arrays.fill(outputLate, 0);
        
        //sort tasks by increasing deadline to minimize the lateness
        int[] sortedIndexes = IntStream.range(0, deadline.length).boxed()
                .sorted((i, j) -> {
                    int c = Double.compare(deadline[i], deadline[j]);
                    if (c != 0) return c;
                    return Double.compare(duration[i], duration[j]);
                })
                .mapToInt(ele -> ele)
                .toArray();

        double f_prev = 0; // f is the finish time of previous task
        int i; 
        for (int ii = 0; ii < duration.length; ++ii) {
            //assign task i to start at
            i = sortedIndexes[ii];
            outputStart[i] = f_prev;  // start next task
            f_prev = /*f[i] =*/ outputStart[i] + duration[i];  // its finish time
            //lateness[i] = max(0, f[i] - d[i])     // its lateness
            outputLate[i] = Math.max(0, f_prev - deadline[i]);
        }
        return sortedIndexes;
    }
    
    /**
     * schedule a set of n tasks where each task is associated with a execution time 
     * t_i and a deadline d_i. 
     * The objective is to schedule the tasks, no two overlapping in time, 
     * such that they are all completed before their deadline. 
     * If this is not possible, define the lateness of the ith task to be amount 
     * by which its finish time exceeds its deadline. 
     * The objective is to minimize the maximum lateness over all the tasks.
     * 
     * The algorithm is aka Earliest Finish First (EFF) and Earliest Deadline First (EDF)
     *    (1) sort tasks by finish time
     *    (2) iterate over tasks, scheduling each that starts after the previous ended
     * 
     * References:
     * <pre>
     * lecture 7 notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect07-greedy-sched.pdf
     * </pre>
     * 
     * runtime complexity is O(N * log_2(N)).
     * This is a greedy and optimal solution.
     * 
     @param s start times for tasks
     @param f finish times for tasks
     @return indexes for scheduled non-conflicting tasks
     */
    public int[] unweightedIntervalNoConflicts(double[] s, double[] f) {
        int n = s.length;
        
        //sort tasks by increasing finish times .  O(N * log_2(N))
        int[] sortedIndexes = IntStream.range(0, n).boxed()
                .sorted((i, j) -> {
                    int c = Double.compare(f[i], f[j]);
                    if (c != 0) return c;
                    return Double.compare(s[i], s[j]);
                })
                .mapToInt(ele -> ele)
                .toArray();

        double f_prev = -1; // f is the finish time of previous task
        int i; 
        int[] scheduled = new int[n];
        int count = 0;
        for (i = 0; i < n; ++i) {
            if (s[i] > f_prev) {
                scheduled[count++] = sortedIndexes[i];
                f_prev = f[i];
            }
        }
        scheduled = Arrays.copyOfRange(scheduled, 0, count);
        return scheduled;
    }

    /**
     * The objective is to compute any maximum sized subset of non-overlapping intervals that produce the highest
     * sum of values (profits).

     * The algorithm uses dynamic programming and has runtime complexity O(n^2).
     * 
     * The problem is adapted from the lecture notes of David Mount for CMSC 451
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect10-dp-intv-sched.pdf
     *
     * The algorithm solution here is a simpler dynamic program.
     * 
       <pre>
         Roughly:
         
           Through these 2 short examples, one can see that a dynamic solution
           avoiding exponential comparisons of every permutation through re-use
           the answers from sub-problems, should be possible.
           
          First, sort tasks by finish time.

          example:
            0 ---------|
            1    ------------|*
            2             -----|
            3          ------------|
            0+2 is possible. store total weight.
            0+3 is possible and has larger weight than 0+2. store total weight.

          Mount's example:
              0  1  2  3  4  5  6  7  8  9
           0  ---------|
           1     ------------|*
           2           --------|
           3        ---------------|
           4                   ------|*
           5                       -----|
             indexes that can be appended after 0: 2,4,5
             indexes that can be appended after 1: 4,5
             indexes that can be appended after 2: 4,5
             indexes that can be appended after 3: 5
             indexes that can be appended after 4:

             start from i=5.  best combination = [5], weight=w[5]
                        i=4.  best combination = [4], weight=w[4]
                        i=3.  best combination = [3,5], weight=w[3]+memo[5]
                        i=2.  combinations max([2,4], [2,5]) = max(w[2]+memo[4], w[2]+memo[5])
                        i=1.  combinations max([1,4], [1,5]) = max(w[1]+memo[4], w[1]+memo[5])
                        i=0.  combinations max([0,2], [0,4], [0,5]) = max(w[0]+memo[2], w[0]+memo[4], w[0]+memo[5])

           Then use dynamic programming with tabulation to store the best schedule starting at n-1
              as scheduline including that task and its profits.
           Then proceed to n-2 and find the max profit schedule already calcuated for higher indexes that can
              fit after this task's finish time.
           repeating until reach first task.

           Then the highest profit schedule in tabulation is kept.
     * </pre>
     * 
     @param s interval start times
     @param f interval finish times
     @param v interval weights
     @return indexes of scheduled intervals.
     */
    public int[] weightedIntervalBottomUp(double[] s, double[] f, double[] v) {
        //interval [si, fi] of start and finish times

        int n = f.length;

        // ascending order sort by f
        // runtime complexity is O(n*log_2(n))
        int[] sortedIdxs = IntStream.range(0,n).boxed()
                .sorted((i, j)-> Double.compare(f[i], f[j]))
                .mapToInt(ele -> ele)
                .toArray();

        //System.out.printf("sorted indexes=%s\n", Arrays.toString(sortedIdxs));

        // use dynamic programming with tabulation
        double[] tabProfits = new double[n];
        Map<Integer, Set<Integer>> tabIndices = new HashMap<>();

        int i;
        int j;
        double max;
        int jjMax;

        // runtime complexity is O(n^2)
        for (int ii = n - 1; ii >= 0; --ii) {
            i = sortedIdxs[ii];
            tabIndices.putIfAbsent(ii, new HashSet<>());
            tabIndices.get(ii).add(i);
            tabProfits[ii] = v[i];

            max = Double.NEGATIVE_INFINITY;
            jjMax = -1;
            for (int jj = n-1; jj > ii; --jj) { // tabProfits[jj] will already exist and hold best max sum for its part of the schedule to end
                j = sortedIdxs[jj];
                if (s[j] >= f[i]) { // task j can be appended after task i
                    if (tabProfits[jj] > max) {
                        jjMax = jj;
                        max = tabProfits[jj];
                    }
                }
            }// end loop over jj

            if (jjMax > -1) {
                tabProfits[ii] += tabProfits[jjMax];
                tabIndices.get(ii).addAll( tabIndices.get(jjMax) );
            }
        }//end loop over i

        // choose sched w/ max profit:
        max = Double.NEGATIVE_INFINITY;
        jjMax = -1;
        for (i = 0; i < n; ++i) {
            if (tabProfits[i] > max) {
                max = tabProfits[i];
                jjMax = i;
            }
        }

        int[] sched = new int[tabIndices.get(jjMax).size()];
        i = 0;
        for (int idx : tabIndices.get(jjMax)) {
            sched[i++] = idx;
        }

        Arrays.sort(sched);

        return sched;
    }

    /**
     * given one machine and n tasks (where n = duration.length and the task properties
     * are duration, deadline and profit v) find a schedule which maximizes the summed
     * profits v.  A profit v_i is only received for a task finished before its deadline,
     * else there is not a penalty for lateness but no sum is added to the total profit,
     * so all tasks should be scheduled if possible.  The machine can only process one task at a time
     * and without interruption (no preemption).

     The algorithm is optimal.
     This algorithm uses dynamic programming.
     The runtime complexity is (n^2).

     The algorithm arises from question 34-4 of Cormen, Leiserson, Rivest, and Stein,
     "Introduction to Algorithms", fourth edition

     * @param duration non-negative amount of times to complete each task
     * @param deadline non-negative deadlines for each task
     * @param v non-negative profits for each task completed before its deadline.
     * @param outputSchedule output array of length n to be populated by this algorithm with the order for scheduling tasks.
     * @param outLastOnTimeIdx output array of length 1 holding the index of outputSchedule which
     *                      is the last task in the schedule that completes before its deadline.
     * @return the summed profits for the tasks scheduled which will complete on time.
     */
    public static double weightedDynamicSingleResource(int[] duration, double[] deadline, double[] v, int[] outputSchedule, int[] outLastOnTimeIdx) {
        int n = duration.length;
        if (deadline.length != n) {
            throw new IllegalArgumentException("deadline.length must equal duration.length");
        }
        if (v.length != n) {
            throw new IllegalArgumentException("v.length must equal duration.length");
        }
        if (outputSchedule.length != n) {
            throw new IllegalArgumentException("outputSchedule.length must equal duration.length");
        }
        if (outLastOnTimeIdx.length != 1) {
            throw new IllegalArgumentException("outLastOnTimeIdx.length must equal 1");
        }

        // ascending order sort by deadline
        // runtime complexity is O(n*log_2(n))
        int[] sortedIndexes = IntStream.range(0, n).boxed()
                .sorted((i, j) -> Double.compare(deadline[i], deadline[j]))
                .mapToInt(ele -> ele).toArray();

        // see notes on this algorithm at bottom of method

        // interval [si, fi] of start and finish times
        // where finish time f is calculated as time + duration of task.

        // tabProfits key = tab index, value = map w/ key=f, value=summed profits (==summed v).
        //     f is a summed property too.
        TIntObjectMap<TIntDoubleMap> tabProfits = new TIntObjectHashMap<TIntDoubleMap>();

        // tabProfitsF key = tab index, value = set of fs added to this schedule. the f values are used
        // to find summed profits in the tabProfits map.
        TIntObjectMap<TIntSet> tabProfitsF = new TIntObjectHashMap<TIntSet>();

        // populating tabProfits and tabProfitsF is at worst 2*(2^n))

        // dynamic programming using maps instead of dense matrix

        int f;
        double sumP;

        int count = 0;

        TIntIterator iter;
        int fPrev;
        double pPrev;
        for (int ii = 0; ii < n; ++ii) {
            final int i = sortedIndexes[ii];

            // update tabProfits and tabProfitsF
            tabProfits.put(ii, new TIntDoubleHashMap()); // map has key=f, value=v
            tabProfitsF.put(ii, new TIntHashSet());

            if (ii == 0) {
                if (duration[i] <= deadline[i]) {
                    // include if possible
                    f = duration[i];
                    sumP = v[i];
                    tabProfitsF.get(ii).add(f);
                    tabProfits.get(ii).put(f, sumP);
                }
                // exclude
                f = 0;
                sumP = 0;
                tabProfitsF.get(ii).add(f);
                tabProfits.get(ii).put(f, sumP);

                count += 2;
                continue;
            }

            ++count;

            iter = tabProfitsF.get(ii-1).iterator();
            while (iter.hasNext()) {
                fPrev = iter.next();
                pPrev = tabProfits.get(ii - 1).get(fPrev);
                // tentative f
                f = fPrev + duration[i];

                count++;

                if (f <= deadline[i]) {
                    // === include task i ====
                    sumP = pPrev + v[i];
                    if (tabProfits.get(ii).containsKey(f)) {
                        // if entry already exists, take max of this and that
                        if (tabProfits.get(ii).get(f) < sumP) {
                            tabProfits.get(ii).put(f, sumP);
                        }
                    } else {
                        tabProfits.get(ii).put(f, sumP);
                        tabProfitsF.get(ii).add(f);
                    }
                }

                // === exclude task i by bringing in previous values and not adding this i to it ====
                // by storing fPrev and pPrev after a check for existing entry in current tabProfits for i
                if (tabProfits.get(ii).containsKey(fPrev)) {
                    if (tabProfits.get(ii).get(fPrev) < pPrev) {
                        tabProfits.get(ii).put(fPrev, pPrev);
                    }
                } else {
                    tabProfits.get(ii).put(fPrev, pPrev);
                    tabProfitsF.get(ii).add(fPrev);
                }
            } // end loop over f
        } // end loop over i

        System.out.printf("count=%d\n", count);

        // get max of tabProfits[n-1]
        int maxF = -1;
        double p;
        double maxP = Double.NEGATIVE_INFINITY;
        TIntDoubleIterator iter2 = tabProfits.get(n - 1).iterator();
        while (iter2.hasNext()) {
            iter2.advance();
            f = iter2.key();
            p = iter2.value();
            if (p > maxP) {
                maxP = p;
                maxF = f;
            }
        }

        // recover indexes for schedule with backtracking
        f = maxF;
        p = maxP;
        TIntList sched = new TIntArrayList();
        for (int ii = n - 1; ii >= 0; --ii) {
            int i = sortedIndexes[ii];
            if (p <= 0) {
                break;
            }
            if ((ii-1 >= 0) && p == tabProfits.get(ii-1).get(f)) {
                continue;
            } else {
                sched.add(i);
                p -= v[i];
                f -= duration[i];
            }
        }
        assert(f == 0);

        sched.reverse();

        outLastOnTimeIdx[0] = sched.size() - 1;

        if (sched.size() < n) {
            TIntSet schedSet = new TIntHashSet(sched);
            // add the remaining tasks.  since they are ordered by deadline already,
            // this will minimize the lateness
            for (int i = 0; i < n; ++i) {
                if (!schedSet.contains(i)) {
                    sched.add(i);
                }
            }
        }

        for (int i = 0; i < n; ++i) {
            outputSchedule[i] = sched.get(i);
        }

        return maxP;
    }

    // s and f are sorted by ascending order of f before passed to this method
    // runtime complexity is O(n^2)
    private int[] calcP(double[] s, double[] f, int[] sortedIndexes) {
        // iterating from highest index to lowest,
        // find for each s, highest previous index in which f[i-...] < s_i
        int i, j;
        int[] p = new int[f.length+1];
        for (int ii = s.length - 1; ii > -1; ii--) {
            i = sortedIndexes[ii];
            for (int jj = ii - 1; jj > -1; jj--) {
                j = sortedIndexes[jj];
                //System.out.printf("calcP: %d,%d) f[%d]=%.2f s[%d]=%.2f\n", i,j, j, f[j], i, s[i]);
                if (f[j] <= s[i]) {
                    p[i] = j+1;
                    //System.out.printf("   p[%d]=%d\n", i, p[i]);
                    break;
                }
            }
        }
        return p;
    }
    
    /**
     * sect 16.5 of Cormen, Leiserson, Rivest, and Stein
     * 
     * there are numTasks number of tasks, each of which takes 1 unit of time
     * and has its own deadline and penalty for missing the deadline.
     * 
     * minimize the total penalty incurred for missed deadlines.
     * 
     * early tasks: finishes before or at deadline.
     * 
     * late tasks: finish after their deadlines.
     * 
     * early-first form:  early tasks precede late tasks.
     * 
     * canonical form: 
     *     early tasks precede late tasks
     *     and early tasks are in monotonically increasing order of deadlines.
     *     (1) put schedule in early first form
     *     (2) swap sequential pairs in the early list when d_{k} .gt. d_{k+1}
     *     (3) list the early tasks
     *     (4) list the late tasks in any order
     * 
     * if no tasks are late, the set is independent.
     * 
     * the early set by themselves is an independent set.
     * 
     * let L = set of all sets on independent tasks.
     * 
     * N_t(A) = number of tasks t=0,1,2...n in set A whose deadline is t or earlier.
     * 
     * if N_t(A) .gt. t then there is no way to schedule all tasks within deadline.
     * 
     * the problem of maximizing the sum of penalties for the early tasks
     *  is the same as minimizing the sum of penalties for the late tasks.
     * 
     * algorithm with r.t. O(N * log_2(N)):
     * 
     * (1) Use the Greedy algorithm to find a maximum weight independent set of
     *     tasks A.
     * (2) create an optimal schedule having the tasks in A as its early tasks.
    
     @param deadlines values must be between 1 and numTasks, inclusive
     @param penalties penalties for missing a deadline
     @return order of scheduled tasks
     */
    public int[] weightedGreedySingleResource(int[] deadlines, int[] penalties) {
        
        //O(N * log_2(N))
        // schedule the largest penalties first if they fit before deadline (duration of each task is a time unit of 1).
        // the result is a list that may be < n in length
        int[] indexes2 = greedy(deadlines, penalties);

        //System.out.println("greedy algorithm resulting indexes=" + Arrays.toString(indexes2));

        //O(N * log_2(N))
        // sort by ascending deadlines, break ties by descending penalities
        int[] sortedIndexes2 = IntStream.range(0, indexes2.length).boxed()
                .sorted((ii, jj) -> {
                    int c = Double.compare(deadlines[indexes2[ii]], deadlines[indexes2[jj]]);
                    if (c != 0) return c;
                    return Double.compare(penalties[indexes2[jj]],  penalties[indexes2[ii]]);
                })
                .mapToInt(ele -> ele).toArray();

        int[] scheduled2 = new int[deadlines.length];

        TIntSet allI2 = new TIntHashSet();
        int i;
        for (i = 0; i < deadlines.length; ++i) {
            allI2.add(i);
        }

        for (i = 0; i < sortedIndexes2.length; ++i) {
            // transform indexes back to original array indexes
            scheduled2[i] = sortedIndexes2[i];
            allI2.remove(scheduled2[i]);
        }
        //System.out.println("appending late tasks in any order:");
        TIntIterator iter = allI2.iterator();
        while (iter.hasNext()) {
            scheduled2[i] = iter.next();
            i++;
        }

        return scheduled2;
    }
  
    /**
     * schedule the largest penalties first if they fit before deadline (duration of
     * each task is a time unit of 1).
     * runtime complexity is O(N * log_2(N))
     @param deadlines deadlines of each task
     @param penalties penalties of each task
     @return task schedule order
     */
    private int[] greedy(int[] deadlines, int[] penalties) {

        int n = deadlines.length;

        //O(N * log_2(N))
        // sort w, m into monotonically decreasing order by penalties
        int[] sortedIndexes = IntStream.range(0, n).boxed()
                .sorted((i, j) -> {return penalties[j] - penalties[i];})
                .mapToInt(ele -> ele)
                .toArray();

        /*
        System.out.println("sorted by decr penalty:");
        for(i = 0; i < penalties.length; ++i) {
            oIdx = origIndexes[i];
            System.out.printf("a%d deadline=%d penalty=%d\n", oIdx+1, deadlines[i], penalties[i]);
        }
        */

        //TODO: consider overloading method for durations.  then need an array b to complement a, holding integrated
        // durations.
       
        TIntList a = new TIntArrayList();
        int prevA = -1;
        int i;
        for (int ii = 0; ii < n; ++ii) {
            i = sortedIndexes[ii];
            //System.out.printf("i=%d, oIdx=%d, a%d f_i=%d, (%d,%d): ", i, oIdx, oIdx+1, (a.size()+1), deadlines[i], penalties[i]);
            if (deadlines[i] >= (a.size()+1)) {
                //done early
                a.add(i);
                prevA = i;
                //System.out.println("  early accept");
            } else if (!a.isEmpty() 
                && (deadlines[prevA] > deadlines[i]) && (deadlines[prevA] >= (a.size() + 1))
                ) {
                a.add(i);
                prevA = i;
                //System.out.println("  accept");
            } /*else {
                System.out.println("  reject");
            }*/
        }
        int[] ao = new int[a.size()];
        for (i = 0; i < ao.length; ++i) {
            ao[i] = a.get(i);
        }
        return ao;
    }

    /**
    Interval Partitioning:
    Given an infinite number of possible exclusive resources to use, 
    schedule all the activities using the smallest number of resources.
    The activity requests each have a start and finish time.
    Let the resources be a collection R, partitioned into d disjoint subsets R_0,...R_{d-1}
    such that events of R_j are mutually non-conflicting, for each j: 0 ≤ j ≤ (d-1).

    References:
    <pre>
    lecture 7 notes of David Mount for CMSC 451       
    Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
    https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect07-greedy-sched.pdf
    </pre>

     runtime complexity is O(n^2).
     The approach is greedy but produces an optimal solution.

     <pre>
     https://en.wikipedia.org/wiki/Interval_scheduling
     Interval scheduling is a class of problems in computer science, particularly in the area of
     algorithm design. The problems consider a set of tasks. Each task is represented by an
     interval describing the time in which it needs to be processed by some machine
     (or, equivalently, scheduled on some resource).
     </pre>

     @param s start times
     @param f finish times
     @return indexes of resources to schedule the requests on.
     */
    public int[] intervalPartitionGreedy2(double[] s, double[] f) {

        /*
        (1) sort the requests by increasing order of start times. 
        (2) assign to each request the smallest color (possibly a new color) 
            such that it conflicts with no other requests of this color class. 
        */

        int n = s.length;
        
        // runtime complexity O(n*log_2(n))
        //sort requests by increasing start times
        int[] sortedIndexes = IntStream.range(0, n).boxed()
                .sorted((i, j) -> Double.compare(s[i], s[j]))
                .mapToInt(ele -> ele).toArray();

        // a color for each request
        int[] c = new int[s.length];
        
        TIntSet excl = new TIntHashSet();
        int color;
        int i, j;
        for (int ii = 0; ii < n; ++ii) {
            excl.clear();
            i = sortedIndexes[ii];
            for (int jj = 0; jj < ii; ++jj) {
                //j is always smaller than i so s[j] <= s[i].
                //  then order is (sj,fj)  (si,fi)
                j = sortedIndexes[jj];
                //if ([s[j],f[j]] overlaps [s[i],f[i]]) 
                if (s[i] < f[j]) {
                    excl.add(c[j]);
                }
            }
            //Let c be the smallest color NOT in E
            for (color = 0; color < n; ++color) {
                if (!excl.contains(color)) {
                    break;
                }
            }
            c[i] = color;
        }

        return c;
    }

    /**
     Interval Partitioning:
     Given an infinite number of possible exclusive resources to use,
     schedule all the activities using the smallest number of resources.
     The activity requests each have a start and finish time.
     Let the resources be a collection R, partitioned into d disjoint subsets R_0,...R_{d-1}
     such that events of R_j are mutually non-conflicting, for each j: 0 ≤ j ≤ (d-1).

     runtime complexity is O(n*log_2(n)).

     The algorithm used is sometimes called EFF (early finish first) for unlimited resources.

     <pre>
     https://en.wikipedia.org/wiki/Interval_scheduling
     Interval scheduling is a class of problems in computer science, particularly in the area of
     algorithm design. The problems consider a set of tasks. Each task is represented by an
     interval describing the time in which it needs to be processed by some machine
     (or, equivalently, scheduled on some resource).
     </pre>
     *
     @param s start times
     @param f finish times
     @return indexes of resources to schedule the requests on.
     */
    public int[] intervalPartitionGreedy(double[] s, double[] f) {
        int n = s.length;
        if (f.length != n) {
            throw new IllegalArgumentException("s and f must be same length");
        }

        // runtime complexity O(n*log(n))
        //sort requests by increasing finish times
        int[] sortedIndexes = IntStream.range(0, n).boxed()
                .sorted((i, j) -> Double.compare(f[i], f[j]))
                .mapToInt(ele -> ele).toArray();

        List<List<Integer>> resources = new ArrayList<List<Integer>>();
        List<Integer> current = new ArrayList<>();
        resources.add(current);
        double fPrev = -1;
        int i;
        //O(n)
        for (int ii = 0; ii < n; ++ii) {
            i = sortedIndexes[ii];
            if (s[i] < fPrev) {
                current = new ArrayList<>();
                resources.add(current);
            }
            current.add(i);
            fPrev = f[i];
        }
        //O(n)
        // merge resources, bottom up
        double start, finish;
        int jj, j;
        for (int ii = resources.size() - 1; ii > 0; --ii) {
            current = resources.get(ii);
            start = s[current.get(0)];
            // search lower indexes for feasible concatenation
            for (jj = ii - 1; jj >= 0; --jj) {
                List<Integer> lower = resources.get(jj);
                finish = f[lower.get(lower.size() - 1)];
                if (start >= finish) {
                    // append current to lower
                    lower.addAll(current);
                    resources.remove(ii);
                    break;
                }
            }
        }

        //O(n)
        //for each index, write the scheduled resource number
        int[] schedResource = new int[n];
        for (i = 0; i < resources.size(); ++i) {
            current = resources.get(i);
            for (j = 0; j < current.size(); ++j) {
                schedResource[current.get(j)] = i;
            }
        }
        return schedResource;
    }

    /**
     Given exclusive use of 1 resource, find the maximum size of mutually compatible activities
     S = {a1, a2, ... an} where each activity has a start time si and finish time fi.
     0 <= si < fi < inf.

     References:
     <pre>
     Chap 15 of "Introduction to Algorithms", fourth edition
     Cormen, Leiserson, Rivest, and Stein,
     </pre>

     runtime complexity is O(n) if already sorted, else O(n*log_2(n)).

     <pre>
     https://en.wikipedia.org/wiki/Interval_scheduling
     Interval scheduling is a class of problems in computer science, particularly in the area of
     algorithm design. The problems consider a set of tasks. Each task is represented by an
     interval describing the time in which it needs to be processed by some machine
     (or, equivalently, scheduled on some resource).
     </pre>
     *
     @param s start times
     @param f finish times
     @param isSortedByF
     @return indexes of resources to schedule the requests on.
     */
    public static int[] intervalPartitionGreedySingleResource(double[] s, double[] f, boolean isSortedByF) {
        int n = f.length;
        if (s.length != n) {
            throw new IllegalArgumentException("s.length must equal f.length");
        }

        int[] sortedIndexes = null;
        if (isSortedByF) {
            sortedIndexes = new int[n];
            for (int i = 0; i < n; ++i) {
                sortedIndexes[i] = i;
            }
        } else {
            sortedIndexes = IntStream.range(0, n).boxed()
                    .sorted((i, j) -> Double.compare(f[i], f[j]))
                    .mapToInt(ele -> ele)
                    .toArray();
        }

        TIntList a = new TIntArrayList();
        a.add(sortedIndexes[0]);
        int i, j = sortedIndexes[0];
        for (int ii = 1; ii < n ; ++ii) {
            i = sortedIndexes[ii];
            if (s[i] >= f[j]) {
                a.add(i);
                j = i;
            }
        }

        return a.toArray();
    }

}
