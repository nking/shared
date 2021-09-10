package algorithms.scheduling;

import algorithms.sort.MiscSorter;
import algorithms.util.FormatArray;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import java.util.SortedSet;
import java.util.TreeSet;

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
     * If this is not possible, define the lateness of the ith task to be amount 
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
     * @param t duration of task
     * @param d deadline for task
     * @param outputStart output array to hold start times for the resulting scheduled index order 
     * @param outputLate output array to hold lateness for the resulting scheduled index order
     * @return indexes for scheduling order
     */
    public int[] unweightedIntervalMinimizeLateGreedy(double[] t, double[] d,
        double[] outputStart, double[] outputLate) {
        int n = t.length;
        if (d.length != n) {
            throw new IllegalArgumentException("d.length must equal t.length");
        }
        if (outputStart.length != n) {
            throw new IllegalArgumentException("outputStart.length must equal t.length");
        }
        if (outputLate.length != n) {
            throw new IllegalArgumentException("outputLate.length must equal t.length");
        }
        t = Arrays.copyOf(t, t.length);
        d = Arrays.copyOf(d, d.length);
        Arrays.fill(outputStart, 0);
        Arrays.fill(outputLate, 0);
        
        //sort tasks by increasing deadline 
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(d, t);
        double f_prev = 0; // f is the finish time of previous task
        int i; 
        for (i = 0; i < t.length; ++i) {
            //assign task i to start at 
            outputStart[i] = f_prev;  // start next task
            f_prev = /*f[i] =*/ outputStart[i] + t[i];        // its finish time
            //lateness[i] = max(0, f[i] - d[i])     // its lateness
            outputLate[i] = Math.max(0, f_prev - d[i]);
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
     * 
     * References:
     * <pre>
     * lecture 7 notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect07-greedy-sched.pdf
     * </pre>
     * 
     * @param s start times for tasks
     * @param f finish times for tasks
     * @return indexes for scheduled non-conflicting tasks
     */
    public int[] unweightedIntervalNoConflicts(double[] s, double[] f) {
        int n = s.length;
        s = Arrays.copyOf(s, s.length);
        f = Arrays.copyOf(f, f.length);
        
        //sort tasks by increasing finish times 
        int[] indexes = MiscSorter.mergeBy1stArgThen2nd(f, s);
        double f_prev = Double.NEGATIVE_INFINITY; // f is the finish time of previous task
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
     * This code uses dynamic programming and has runtime complexity less than O(n^2).
     * 
     * The code follows the lecture notes of David Mount for CMSC 451 
     * Design and Analysis of Computer Algorithms (with some corrections for pseudocode indexes).
     * https://www.cs.umd.edu/class/fall2017/cmsc451-0101/Lects/lect10-dp-intv-sched.pdf
     * 
     * @param s interval start times
     * @param f interval finish times
     * @param v interval weights
     * @return indexes of scheduled intervals.
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
        
        // p is largest index p such that f_p < s_n
        // p(j) is the largest index such that f_p(j) < s_j
        // runtime complexity is less than O(n^2)
        int[] p = calcP(s, f);
        
        int[] pred = new int[n+1];
        
        double[] M = new double[n+1];
        Arrays.fill(M, -1);
        M[0] = 0;
        int j;
        double leaveWeight, takeWeight;
        // runtime complexity is O(n)
        for (j = 0; j < n; ++j) {
            leaveWeight = M[j];                // total weight if we leave j
            takeWeight = v[j] + M[p[j]];         // total weight if we take j
        //    System.out.printf("j=%d lw=M[j]=%.2f tw=v[j]+M[p[j]]=%.2f+%.2f=%.2f (where p[j]=%d) ", 
        //        j, leaveWeight, v[j], M[p[j]], takeWeight, p[j]);
            if (leaveWeight > takeWeight) {
                M[j + 1] = leaveWeight;              // better to leave j
                pred[j+1] = j;                   // previous is j-1
            } else {
                M[j + 1] = takeWeight;               // better to take j
                pred[j+1] = p[j];                  // previous is p[j]
            }
        //    System.out.printf("  M[j+1]=%.2f\n", M[j+1]);
        }
         
        //System.out.printf("M=%s\n", FormatArray.toString(M, "%.3f"));
        //System.out.printf("p=%s\n", Arrays.toString(p));
        //System.out.printf("pred=%s\n", Arrays.toString(pred));
        
        j = pred.length-1;
        int[] sched = new int[j];
        int count = 0;
        while (j > 0) {
            if (pred[j] == p[j-1]) {
                sched[count++] = origIndexes[j];
            }
            j = pred[j];
        }
        sched = Arrays.copyOfRange(sched, 0, count);
        return sched;
    }

    
    private int[] calcP(double[] s, double[] f) {
        // iterating from highest index to lowest,
        // find for each s, highest previous index in which f[i-...] < s_i
        int i, j;
        int[] p = new int[f.length+1];
        //Arrays.fill(p, -1);
        for (i = s.length - 1; i > -1; i--) {
            for (j = i - 1; j > -1; j--) {
                //System.out.printf("%d,%d) f[%d]=%.2f s[%d]=%.2f\n", i,j, j, f[j], i, s[i]);
                if (f[j] <= s[i]) {
                    p[i] = j+1;
                    break;
                }
            }
        }
        return p;
    }
    
    /**
     * sect 16.5 of Cormen et al.
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
     *     (2) swap sequential pairs in the early list when d_{k} > d_{k+1}
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
     * if N_t(A) > t then there is no way to schedule all tasks within deadline.
     * 
     * the problem of maximizing the sum of penalties for the early tasks
     *  is the same as minimizing the sum of penalties for the late tasks.
     * 
     * algorithm with r.t. O(n^2):
     * (1) Use the Greedy algorithm to find a maximum weight independent set of
     *     tasks A.
     * (2) create an optimal schedule having the tasks in A as its early tasks.
    
     * @param deadlines values must be between 1 and numTasks, inclusive
     * @param penalties
     * @return 
     */
    public int[] weightedGreedy(int[] deadlines, int[] penalties) {
        
        //deadlines = Arrays.copyOf(deadlines, deadlines.length);
        //penalties = Arrays.copyOf(penalties, penalties.length);
        
        int[] d2 = Arrays.copyOf(deadlines, deadlines.length);
        int[] p2 = Arrays.copyOf(penalties, penalties.length);
        
        //System.out.println("starting greedy algorithm to find indep sets");
        int[] indexes2 = greedy(d2, p2);
        //System.out.println("greedy algorithm resulting indexes=" + Arrays.toString(indexes2));
        
        // rewrite d2 and p2 to be only the greedy results
        d2 = new int[indexes2.length];
        p2 = new int[indexes2.length];
        int[] i2 = new int[indexes2.length];
        int i;
        for (i = 0; i < d2.length; ++i) {
            d2[i] = deadlines[indexes2[i]];
            p2[i] = penalties[indexes2[i]];
            i2[i] = indexes2[i];
        }
        
        // sort by increasing deadlines
        //these indexes are w.r.t. the truncated d2 and p2, that is, i2
        indexes2 = mergesortIncreasingADecreasingB(d2, p2); 
        //expecting: (2,60), (3,40), (4,70), (4,50), (6,10), (1,30), (4,20)
                   
        //System.out.println("sorted by increasing deadline:");
        TIntSet allI = new TIntHashSet();
        for (i = 0; i < deadlines.length; ++i) {
            allI.add(i);
        }
        int[] scheduled = new int[deadlines.length];
        for (i = 0; i < indexes2.length; ++i) {
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
  
    private int[] greedy(int[] deadlines, int[] penalties) {
               
        deadlines = Arrays.copyOf(deadlines, deadlines.length);
        penalties = Arrays.copyOf(penalties, penalties.length);
        
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
       
        // the schedule of indexes, in a datastructure that sorts upon insert
        SortedSet<Integer> a = new TreeSet<>();
        for(i = 0; i < deadlines.length; ++i) {
            oIdx = origIndexes[i];
            //System.out.printf("a%d f_i=%d, (%d,%d): ", oIdx+1, (a.size()+1), deadlines[i], penalties[i]);
            if (deadlines[i] >= (a.size()+1)) {
                //done early
                a.add(i);
            //    System.out.println("  accept");
            } else if ((a.size() > 0) && 
                (deadlines[i] < deadlines[a.last()])) {
                
                // check whether this increase of start time by 1 unit would push
                //    out the last item (which is the same last item that we just
                //    compared in this conditional clause.
                //    if that were true, do not add this item as it conflicts.
                if (deadlines[a.last()] >= (a.size() + 1)) {
                    a.add(i);
            //        System.out.println("  accept");
            //    } else {
            //        System.out.println("  reject");
                }                
            //} else {
            //    System.out.println("  reject");
            }
        }
        // rewrite indexes in context of original arrays:
        int[] ao = new int[a.size()];
        int count = 0;
        for (Integer ai : a) {
            ao[count++] = origIndexes[ai];
        }
        return ao;
    }

    private int[] sort2(double[] a, double[] b, double[] c) {
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
    
    private void mergesort(double[] a, int[] b, int idxLo, int idxHi) {
        if (idxLo < idxHi) {
            int idxMid = (idxLo + idxHi) >> 1;
            mergesort(a, b, idxLo, idxMid);           
            mergesort(a, b, idxMid + 1, idxHi);       
            merge(a, b, idxLo, idxMid, idxHi);
        }
    }

    private void merge(double[] a, int[] b, int idxLo, int idxMid, int idxHi) {
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

     runtime complexity is between O(n*log_2(n)) and O(n^2), with a worse
     case of O(n^2).
     * 
     * Note: can compare this algorithm to the left-edge algorithms which
     * sorts by increasing finish times, then loops over each request
     * to add all sequential non-conflicting requests to a resource, then
     * start a new resource for a conflict.
     * Then one attempts to merge resources, by visiting them in reverse order.
     * After all requests have been placed in a resource, one attempts to
     * merge the non-conflicting resources by visiting them in reverse order
     * and comparing to the previous resource (that is, starting with the last resource
     * created, and then the one before it, etc until the first is visited).
     * The left-edge algorithm runtime is similar to this interval partitioning greedy algorithm.
     * 
     * @param s start times
     * @param f finish times
     * @return indexes of resources to schedule the requests on.
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
        // runtime complexity < O(n^2), but worse case is O(n^2)
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
