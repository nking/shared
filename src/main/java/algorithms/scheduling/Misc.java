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
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class Misc {
        
     /**
     * schedule a set of n tasks where each task is associated with a execution time 
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
        duration = Arrays.copyOf(duration, duration.length);
        deadline = Arrays.copyOf(deadline, deadline.length);
        Arrays.fill(outputStart, 0);
        Arrays.fill(outputLate, 0);
        
        //sort tasks by increasing deadline to minimize the lateness
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(deadline, duration);

        double f_prev = 0; // f is the finish time of previous task
        int i; 
        for (i = 0; i < duration.length; ++i) {
            //assign task i to start at 
            outputStart[i] = f_prev;  // start next task
            f_prev = /*f[i] =*/ outputStart[i] + duration[i];  // its finish time
            //lateness[i] = max(0, f[i] - d[i])     // its lateness
            outputLate[i] = Math.max(0, f_prev - deadline[i]);
        }
        return indexes;
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
     * 
     @param s start times for tasks
     @param f finish times for tasks
     @return indexes for scheduled non-conflicting tasks
     */
    public int[] unweightedIntervalNoConflicts(double[] s, double[] f) {
        int n = s.length;
        s = Arrays.copyOf(s, s.length);
        f = Arrays.copyOf(f, f.length);
        
        //sort tasks by increasing finish times .  O(N * log_2(N))
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(f, s);
        double f_prev = -1; // f is the finish time of previous task
        int i; 
        int[] scheduled = new int[n];
        int count = 0;
        for (i = 0; i < n; ++i) {
            if (s[i] > f_prev) {
                scheduled[count++] = indexes[i];
                f_prev = f[i];
            }
        }
        scheduled = Arrays.copyOfRange(scheduled, 0, count);
        return scheduled;
    }
    
    /**
     * The objective is to compute any maximum sized subset of non-overlapping intervals.
     * Weighted Interval Scheduling: 
     * given a set S = {1, . . . , n} of n activity requests, 
     * where each activity is expressed as an interval [s_i, f_i] from a given 
     * start time si to a given finish time f_i
     * and each request is associated with a numeric weight or value v_i.
     * 
     * The objective is to find a set of non-overlapping requests such that sum 
     * of values of the scheduled requests is maximum.
     * 
     * This code uses dynamic programming and has runtime complexity O(n^2).
     * 
     * The code follows the lecture notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect10-dp-intv-sched.pdf
     
     @param s interval start times
     @param f interval finish times
     @param v interval weights
     @return indexes of scheduled intervals.
     */
    public int[] weightedIntervalBottomUp(double[] s, double[] f, double[] v) {
        //interval [si, fi] of start and finish times
        s = Arrays.copyOf(s, s.length);
        f = Arrays.copyOf(f, f.length);
        v = Arrays.copyOf(v, v.length);
        
        int n = f.length;
        
        // ascending order sort by f
        // runtime complexity is O(log_2(n))
        int[] origIndexes = sort2(f, s, v);
        //System.out.printf("sorted indexes=%s\n", Arrays.toString(origIndexes));
        
        // p[i] is the largest index such that f[p(i)] < s[i]
        //     p[i] is < i
        // runtime complexity is less than O(n^2)
        int[] p = calcP(s, f);
        
        int[] pred = new int[n+1];
        
        double[] memo = new double[n+1];
        Arrays.fill(memo, -1);
        memo[0] = 0;
        int j;
        double leaveWeight, takeWeight;
        // runtime complexity is O(n)
        for (j = 0; j < n; ++j) {
            leaveWeight = memo[j];                // total weight if we leave j
            takeWeight = v[j] + memo[p[j]];         // total weight if we take j
            //System.out.printf("j=%d lw=M[j]=%.2f tw=v[j]+M[p[j]]=%.2f+%.2f=%.2f (where p[j]=%d) ", 
            //    j, leaveWeight, v[j], memo[p[j]], takeWeight, p[j]);
            if (leaveWeight > takeWeight) {
                memo[j + 1] = leaveWeight;              // better to leave j
                pred[j+1] = j;                   // previous is j-1
            } else {
                memo[j + 1] = takeWeight;               // better to take j
                pred[j+1] = p[j];                  // previous is p[j]
            }
           // System.out.printf("  M[j+1]=%.2f\n", memo[j+1]);
        }
         
        //System.out.printf("memo=%s\n", FormatArray.toString(memo, "%.3f"));
        //System.out.printf("   p=%s\n", Arrays.toString(p));
        //System.out.printf("pred=%s\n", Arrays.toString(pred));
        
        int[] sched = new int[j];
        int count = 0;
        j = pred.length-1;
        while (j > 0) {
            if (pred[j] == p[j-1]) {
                //System.out.printf("  sched[%d]=%d\n", j-1, origIndexes[j-1]);
                sched[count++] = origIndexes[j - 1];
            }
            j = pred[j];
        }
        sched = Arrays.copyOfRange(sched, 0, count);
        return sched;
    }
    
    /**
     * The objective is to compute any maximum sized subset of non-overlapping intervals.
     * Weighted Interval Scheduling: 
     * given a set S = {1, . . . , n} of n activity requests, 
     * where each activity is expressed as an interval [s_i, f_i] from a given 
     * start time si to a given finish time f_i
     * and each request is associated with a numeric weight or value v_i.
     * 
     * The objective is to find a set of non-overlapping requests such that sum 
     * of values of the scheduled requests is maximum.
     * 
     * This code uses dynamic programming and has runtime complexity O(n^2).
     * 
     * The problem is from the lecture notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect10-dp-intv-sched.pdf
     * 
     * His pseudocode is present in the version of this method called 
     * weightedIntervalBottomUp().
     * 
     * runtime complexity is O(n^2)
     * 
     * The version here is a simpler dynamic programming solution:
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

             so memo can be a 1-dimensional array
             can also store the indexes in a map with key=integer and value=integer hashset

        sort tasks by finish time.

        for i=[n-1,0]
          max = int.min
          jmax = -1
          for j=[i+1, n) { // memo[j] will already exist and hold best max sum for its part of the schedule to end
            if (s[j] .geq. f[i]) { // task j can be appended after task i
              if (memo[j] .gt. max) {
                jmax = j;
                max = memo[j];
              }
            }
          }
          set = new hashset int();
          map.put(i, set);
          set.add(i);
          if (jmax==-1) {
            memo[i] = w[i];
          } else {
            memo[i] = w[i] + max;
            set.add(jmax);
          }
        }
     * </pre>
     * 
     @param s interval start times
     @param f interval finish times
     @param v interval weights
     @return indexes of scheduled intervals.
     */
    public int[] weightedIntervalBottomUp2(double[] s, double[] f, double[] v) {
        //interval [si, fi] of start and finish times
        s = Arrays.copyOf(s, s.length);
        f = Arrays.copyOf(f, f.length);
        v = Arrays.copyOf(v, v.length);
        
        int n = f.length;
        
        // ascending order sort by f
        // runtime complexity is O(n*log_2(n))
        int[] origIndexes = sort2(f, s, v);
        //System.out.printf("sorted indexes=%s\n", Arrays.toString(origIndexes));

        // holds max values
        double[] memo = new double[n];
        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();
        
        int i;
        int j;
        double max;
        int jMax;
        TIntSet set;

        // runtime complexity is O(n^2)
        for (i = n - 1; i >= 0; --i) {
          max = Double.NEGATIVE_INFINITY;
          jMax = -1;
          for (j = i + 1; j < n; ++j) { // memo[j] will already exist and hold best max sum for its part of the schedule to end
            if (s[j] >= f[i]) { // task j can be appended after task i
              if (memo[j] > max) {
                jMax = j;
                max = memo[j];
              }
            }
          }// end loop over j
          set = new TIntHashSet();
          map.put(i, set);
          set.add(i);
          if (jMax == -1) {
            memo[i] = v[i];
          } else {
            memo[i] = v[i] + memo[jMax];
            set.add(jMax);
          }
        }//end loop over i
        
        max = Double.NEGATIVE_INFINITY;
        jMax = -1;
        for (i = 0; i < n; ++i) {
            if (memo[i] > max) {
                max = memo[i];
                jMax = i;
            }
        } 
                 
        int[] sched = map.get(jMax).toArray();
        for (i = 0; i < sched.length; ++i) {
            sched[i] = origIndexes[sched[i]];
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

     The runtime complexity for worse case is O(2^n), though the n may be less than the number of tasks
     when tasks cannot be scheduled before their deadlines.

     The algorithm is in NP-complete.
     A proof can be extrapolated from SUBSET-SUM (or more difficultly with VERTEX-COVER) decision problem which
     is already known to be NP-complete.

     The VERTEX-COVER proof follows from the schedule being a maximal independent set and that
     a vertex cover is the complement of the maximal independent set.

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
    public static double weightedOptimal(int[] duration, double[] deadline, double[] v, int[] outputSchedule, int[] outLastOnTimeIdx) {
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

        duration = Arrays.copyOf(duration, n);
        deadline = Arrays.copyOf(deadline, n);
        v = Arrays.copyOf(v, n);

        // ascending order sort by deadline
        // runtime complexity is O(n*log_2(n))

        int[] origIndexes = sort2(deadline, duration, v);
        //System.out.printf("sorted indexes=%s\n", Arrays.toString(origIndexes));

        // see notes on this algorithm at bottom of method

        // interval [si, fi] of start and finish times
        // where finish time f is calculated as time + duration of task.

        // memo needs a composite key of index and f
        // memo value is the summed profit
        TIntObjectMap<TIntDoubleMap> memo = new TIntObjectHashMap<TIntDoubleMap>();

        // fMap key is index, and value is set of fi values stored from evaluation at index.  the fi values are used
        // to find summed profits in the memo map.
        TIntObjectMap<TIntSet> fMap = new TIntObjectHashMap<TIntSet>();

        // populating memo and fMap is at worst 2*(2^n))

        // dynamic programming using maps instead of dense matrix

        int i;
        int f;
        double sumP;
        TIntSet prevFMapSet, fMapSet;
        TIntDoubleMap memoFPMap;

        // init:
        // include i=0 task
        i = 0;
        fMapSet = new TIntHashSet();
        memoFPMap = new TIntDoubleHashMap();
        if (duration[i] <= deadline[i]) {
            f = duration[i];
            sumP = v[i];
            fMapSet.add(f);
            memoFPMap.put(f, sumP);
        }
        // exclude i=0 task
        f = 0;
        sumP = 0;
        fMapSet.add(f);
        memoFPMap.put(f, sumP);

        // update memo and fMap
        memo.put(i, memoFPMap);
        fMap.put(i, fMapSet);

        TIntIterator iter;
        int fPrev;
        double pPrev;
        for (i = 1; i < n; ++i) {
            // get f's from prev index
            prevFMapSet = fMapSet;
            fMapSet = new TIntHashSet();
            memoFPMap = new TIntDoubleHashMap();//int f, double v
            memo.put(i, memoFPMap);
            fMap.put(i, fMapSet);

            iter = prevFMapSet.iterator();
            while (iter.hasNext()) {
                fPrev = iter.next();
                pPrev = memo.get(i - 1).get(fPrev);
                // tentative f
                f = fPrev + duration[i];
                if (f <= deadline[i]) {
                    // === include task i ====
                    sumP = pPrev + v[i];
                    if (memoFPMap.containsKey(f)) {
                        // if entry already exists, take max of this and that
                        if (memoFPMap.get(f) < sumP) {
                            memoFPMap.put(f, sumP);
                        }
                    } else {
                        memoFPMap.put(f, sumP);
                        fMapSet.add(f);
                    }
                }

                // === exclude task i ====
                // by storing fPrev and pPrev after a check for existing entry in current memo for i
                if (memoFPMap.containsKey(fPrev)) {
                    if (memoFPMap.get(fPrev) < pPrev) {
                        memoFPMap.put(fPrev, pPrev);
                    }
                } else {
                    memoFPMap.put(fPrev, pPrev);
                    fMapSet.add(fPrev);
                }
            } // end loop over f
        } // end loop over i

        // get max of memo[n-1]
        memoFPMap = memo.get(n - 1);
        int maxF = -1;
        double p;
        double maxP = Double.NEGATIVE_INFINITY;
        TIntDoubleIterator iter2 = memoFPMap.iterator();
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
        for (i = n - 1; i >= 0; --i) {
            if (p <= 0) {
                break;
            }
            if ((i > 0) && p == memo.get(i-1).get(f)) {
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
            for (i = 0; i < n; ++i) {
                if (!schedSet.contains(i)) {
                    sched.add(i);
                }
            }
        }

        // transform to original indexes
        for (i = 0; i < n; ++i) {
            outputSchedule[i] = origIndexes[sched.get(i)];
        }

        return maxP;
    }

    // s and f are sorted by ascending order of f before passed to this method
    // runtime complexity is O(n^2)
    private int[] calcP(double[] s, double[] f) {
        // iterating from highest index to lowest,
        // find for each s, highest previous index in which f[i-...] < s_i
        int i, j;
        int[] p = new int[f.length+1];
        for (i = s.length - 1; i > -1; i--) {
            for (j = i - 1; j > -1; j--) {
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
    public int[] weightedGreedy(int[] deadlines, int[] penalties) {
        
        //O(N * log_2(N))
        //System.out.println("starting greedy algorithm to find indep sets");
        int[] indexes2 = greedy(Arrays.copyOf(deadlines, deadlines.length), 
            Arrays.copyOf(penalties, penalties.length));
        //System.out.println("greedy algorithm resulting indexes=" + Arrays.toString(indexes2));
        
        // rewrite d2 and p2 to be only the greedy results.   indexes2.length is <= deadline.
        int[] d2 = new int[indexes2.length];
        int[] p2 = new int[indexes2.length];
        int[] i2 = new int[indexes2.length];
        int i;
        for (i = 0; i < d2.length; ++i) {
            d2[i] = deadlines[indexes2[i]];
            p2[i] = penalties[indexes2[i]];
            i2[i] = indexes2[i]; // storing indexes2
        }
        
        // sort by increasing deadlines. O(N2 * log_2(N2))
        //these indexes are w.r.t. the truncated d2 and p2, that is, i2
        indexes2 = mergesortIncreasingADecreasingB(d2, p2);
        
        int[] scheduled = new int[deadlines.length];
        TIntSet allI = new TIntHashSet();
        for (i = 0; i < deadlines.length; ++i) {
            allI.add(i);
        }
        for (i = 0; i < indexes2.length; ++i) {
            // transform indexes back to original array indexes
            scheduled[i] = i2[indexes2[i]];
            //System.out.printf("  a%d (%d, %d)\n", 
            //    scheduled[i]+1, deadlines[scheduled[i]], penalties[scheduled[i]]);
            allI.remove(scheduled[i]);
        }
        //System.out.println("appending late tasks in any order:");
        TIntIterator iter = allI.iterator();
        while (iter.hasNext()) {
            scheduled[i] = iter.next();
            //System.out.printf("  a%d (%d, %d)\n", 
            //    scheduled[i]+1, deadlines[scheduled[i]], penalties[scheduled[i]]);
            i++;
        }
        
        return scheduled;
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

        deadlines = Arrays.copyOf(deadlines, deadlines.length);
        penalties = Arrays.copyOf(penalties, penalties.length);

        //O(N * log_2(N))
        // sort w, m into monotonically decreasing order by penalties w
        int[] origIndexes = sortDecr(penalties, deadlines);
        int i, oIdx;
        
        /*
        System.out.println("sorted by decr penalty:");
        for(i = 0; i < penalties.length; ++i) {
            oIdx = origIndexes[i];
            System.out.printf("a%d deadline=%d penalty=%d\n", oIdx+1, deadlines[i], penalties[i]);
        }
        */
       
        TIntList a = new TIntArrayList();
        int prevA = -1;
        for(i = 0; i < deadlines.length; ++i) {
            //oIdx = origIndexes[i];
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
        // rewrite indexes in context of original array indexes:
        int[] ao = new int[a.size()];
        for (i = 0; i < ao.length; ++i) {
            ao[i] = origIndexes[a.get(i)];
        }
        return ao;
    }

    private static int[] sort2(double[] a, double[] b, double[] c) {
        int[] oIdxs = new int[a.length];
        int i;
        for (i = 0; i < a.length; ++i) {
            oIdxs[i] = i;
        }
        mergesort(a, oIdxs, 0, a.length - 1);
        double[] t = new double[b.length];
        for (i = 0; i < a.length; ++i) {
            t[i] = b[oIdxs[i]];
        }
        System.arraycopy(t, 0, b, 0, b.length);
        for (i = 0; i < a.length; ++i) {
            t[i] = c[oIdxs[i]];
        }
        System.arraycopy(t, 0, c, 0, c.length);
        return oIdxs;
    }

    private static int[] sort2(double[] a, int[] b, double[] c) {

        int[] oIdxs = MiscSorter.mergeSortIncreasing(a);
        int i;
        int[] tb = new int[b.length];
        for (i = 0; i < a.length; ++i) {
            tb[i] = b[oIdxs[i]];
        }
        System.arraycopy(tb, 0, b, 0, b.length);

        double[] tc = new double[b.length];
        for (i = 0; i < a.length; ++i) {
            tc[i] = c[oIdxs[i]];
        }
        System.arraycopy(tc, 0, c, 0, c.length);

        return oIdxs;
    }

    private static void mergesort(double[] a, int[] b, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = (idxLo + idxHi) >> 1;
            mergesort(a, b, idxLo, idxMid);           
            mergesort(a, b, idxMid + 1, idxHi);       
            merge(a, b, idxLo, idxMid, idxHi);
        }
    }

    private static void merge(double[] a, int[] b, int idxLo, int idxMid, int idxHi) {
        double[] aL = Arrays.copyOfRange(a, idxLo, idxMid + 2);
        double[] aR = Arrays.copyOfRange(a, idxMid + 1, idxHi + 2); 
        aL[aL.length - 1] = Double.POSITIVE_INFINITY;
        aR[aR.length - 1] = Double.POSITIVE_INFINITY;

        int[] bL = Arrays.copyOfRange(b, idxLo, idxMid + 2);
        int[] bR = Arrays.copyOfRange(b, idxMid + 1, idxHi + 2); 
        bL[bL.length - 1] = Integer.MAX_VALUE;
        bR[bR.length - 1] = Integer.MAX_VALUE;
        
        int posL = 0;
        int posR = 0;
        for (int k = idxLo; k <= idxHi; k++) {
            if (aL[posL] <= aR[posR]) {
                a[k] = aL[posL];
                b[k] = bL[posL];
                posL++;
            } else {
                a[k] = aR[posR];
                b[k] = bR[posR];
                posR++;
            }
        }
    }

    private int[] mergesortIncreasingADecreasingB(int[] a, int[] b) {
        int[] indexes = new int[a.length];
        int i;
        for (i = 0; i < a.length; ++i) {
            indexes[i] = i;
        }
        mergesortIncreasingADecreasingB(a, b, indexes, 0, a.length-1);
        return indexes;
    }

    private void mergesortIncreasingADecreasingB(int[] a, int[] b, int[] c, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = (idxHi + idxLo)/2;
            mergesortIncreasingADecreasingB(a, b, c, 0, idxMid);
            mergesortIncreasingADecreasingB(a, b, c, idxMid + 1, idxHi);
            mergeIncreasingADecreasingB(a, b, c, idxLo, idxMid, idxHi);
        }
    }

    private void mergeIncreasingADecreasingB(int[] a, int[] b, int[] c, 
        int idxLo, int idxMid, int idxHi) {
        
        int[] aL = Arrays.copyOfRange(a, idxLo, idxMid + 2);
        int[] aR = Arrays.copyOfRange(a, idxMid + 1, idxHi + 2); 
        aL[aL.length - 1] = Integer.MAX_VALUE;
        aR[aR.length - 1] = Integer.MAX_VALUE;

        int[] bL = Arrays.copyOfRange(b, idxLo, idxMid + 2);
        int[] bR = Arrays.copyOfRange(b, idxMid + 1, idxHi + 2); 
        bL[bL.length - 1] = Integer.MAX_VALUE;
        bR[bR.length - 1] = Integer.MAX_VALUE;
        
        int[] cL = Arrays.copyOfRange(c, idxLo, idxMid + 2);
        int[] cR = Arrays.copyOfRange(c, idxMid + 1, idxHi + 2); 
        cL[cL.length - 1] = Integer.MAX_VALUE;
        cR[cR.length - 1] = Integer.MAX_VALUE;
        
        int posL = 0;
        int posR = 0;
        for (int k = idxLo; k <= idxHi; ++k) {
            if (aL[posL] < aR[posR]) {
                a[k] = aL[posL];
                b[k] = bL[posL];
                c[k] = cL[posL];
                posL++;
            } else if (aL[posL] > aR[posR]) {
                a[k] = aR[posR];
                b[k] = bR[posR];
                c[k] = cR[posR];
                posR++;
            } else {
                // they're equal, so break ties by values of b
                if (bL[posL] >= bR[posR]) {
                    a[k] = aL[posL];
                    b[k] = bL[posL];
                    c[k] = cL[posL];
                    posL++;
                } else {
                    a[k] = aR[posR];
                    b[k] = bR[posR];
                    c[k] = cR[posR];
                    posR++;
                }
            }
        }
    }
    
    
    private int[] sortDecr(int[] a, int[] b) {
        int[] oIdxs = new int[a.length];
        for (int i = 0; i < a.length; ++i) {
            oIdxs[i] = i;
        }
        quicksortDecr(a, oIdxs, 0, a.length - 1);
        int[] t = new int[b.length];
        int i;
        for (i = 0; i < a.length; ++i) {
            t[i] = b[oIdxs[i]];
        }
        System.arraycopy(t, 0, b, 0, b.length);
        return oIdxs;
    }
    
    private void quicksortDecr(int[] a, int[] b, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = partitionDecr(a, b, idxLo, idxHi);
            quicksortDecr(a, b, idxLo, idxMid-1);
            quicksortDecr(a, b, idxMid+1, idxHi);
        }
    }
    private int partitionDecr(int[] a, int[] b, int idxLo, int idxHi) {
        int xa = a[idxHi];
        int i = idxLo - 1;  
        int swap;
        for (int j = idxLo; j < idxHi ; j++ ) {
            if (a[j] >= xa) { 
                i++;
                swap = a[i];
                a[i] = a[j];
                a[j] = swap;
                swap = b[i];
                b[i] = b[j];
                b[j] = swap;
            }
        }
        swap = a[i + 1];
        a[i + 1] = a[idxHi];
        a[idxHi] = swap;
        swap = b[i + 1];
        b[i + 1] = b[idxHi];
        b[idxHi] = swap;
        return i + 1;
    }

    /*
    Note: can compare greedy interval partitioning algorithm to the left-edge algorithm which
     * sorts by increasing finish times, then loops over each request
     * to add all sequential non-conflicting requests to a resource, then
     * start a new resource for a conflict.
     * Then one attempts to merge resources, by visiting them in reverse order.
     * After all requests have been placed in a resource, one attempts to
     * merge the non-conflicting resources by visiting them in reverse order
     * and comparing to the previous resource (that is, starting with the last resource
     * created, and then the one before it, etc until the first is visited).
     * The left-edge algorithm runtime is similar to this interval partitioning greedy algorithm.
     */

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

     runtime complexity is O(n^2)
     *
     * 
     @param s start times
     @param f finish times
     @return indexes of resources to schedule the requests on.
     */
    public int[] intervalPartitionGreedy(double[] s, double[] f) {     
        
        double[] s2 = Arrays.copyOf(s, s.length);
        
        /*
        (1) sort the requests by increasing order of start times. 
        (2) assign to each request the smallest color (possibly a new color) 
            such that it conflicts with no other requests of this color class. 
        */
        
        // runtime complexity O(n)
        //sort requests by increasing start times
        int[] indexes = MiscSorter.mergeSortIncreasing(s2);
        double[] f2 = new double[s2.length];
        int i, j;
        for (i = 0; i < f.length; ++i) {
            f2[i] = f[indexes[i]];
        }
        
        //System.out.println("indexes sorted by start times = " + Arrays.toString(indexes));
        
        // a color for each request
        int[] c = new int[s.length];
        
        TIntSet excl;        
        int color;
        // runtime complexity O(n*n/2)
        for (i = 0; i < f2.length; ++i) {
            excl = new TIntHashSet();
            for (j = 0; j < i; ++j) {
                //j is always smaller than i so s[j] <= s[i].  
                //  then order is (sj,fj)  (si,fi)
                
                //if ([s[j],f[j]] overlaps [s[i],f[i]]) 
                if (s2[i] < f2[j]) {
                    excl.add(c[j]);
                    /*System.out.printf("conflict for i2=%d, j2=%d (s2[%d]<f2[%d])=(%.1f, %.2f)\n",
                        i, j, i, j, s2[i], f2[i]);
                    System.out.printf("==> i=%d, j=%d (s2[%d]<f2[%d])=(%.1f, %.2f)\n",
                        indexes[i], indexes[j], indexes[i], indexes[j], 
                        s[indexes[i]], f[indexes[i]]);
                    */
                }
            }
            //Let c be the smallest color NOT in E
            for (color = 0; color < f2.length; ++color) {
                if (!excl.contains(color)) {
                    break;
                }
            }
            c[i] = color;
        }
        //rewrite c in terms of original indexes of method argument's unsorted (s,f)
        int[] c2 = new int[f2.length];
        for (i = 0; i < f2.length; ++i) {
            c2[indexes[i]] = c[i];
        }
        return c2;
    }

}
