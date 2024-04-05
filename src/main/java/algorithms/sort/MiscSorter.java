package algorithms.sort;

import algorithms.misc.MiscMath0;
import algorithms.util.AngleUtil;
import algorithms.util.FormatArray;

import java.util.Arrays;
import java.util.stream.IntStream;

/*
TODO: revise most of this to use lambdas (streaming api) to sort indexes then
use on the multiple arrays.
e.g.
get indexes for a descending sort of arrayA:

int[] sortedIdxs
                = IntStream.range(0, n).boxed()
                .sorted((i, j)-> Integer.compare(arrayA[j], arrayA[i]))
                .mapToInt(ele->ele).toArray();

 */
/**
 * 
   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
   * 
   
  NOTE: for a more efficient version of merge-sort that uses an auxiliary array
  to hold intermediate partial sorting, and uses insertion sort for the
  shortest segments in the divide and conquer pattern, see
  algorithms.correlation.UnivariateDistance._sortCheck().

 * @author nichole
 */
public class MiscSorter {

    public static int[] returnSortedIndexes(int[] a, boolean ascending) {

        // range:
        //     Returns a sequential ordered IntStream from startInclusive (inclusive)
        //     to endExclusive (exclusive) by an incremental step of 1.
        // boxed:
        //     Returns a Stream consisting of the elements of this stream, each boxed to an Integer
        // sorted:
        //     Returns a stream consisting of the elements of this stream, sorted according to the provided Comparator
        // mapToInt:
        //     Returns an IntStream consisting of the results of applying the given function to the elements of this stream
        // toArray:
        //     Returns an array containing the elements of this stream.

        if (ascending) {
            return IntStream.range(0, a.length).boxed()
                    .sorted((i, j)-> Integer.compare(a[i], a[j]))
                    .mapToInt(ele->ele).toArray();
        }
        return IntStream.range(0, a.length).boxed()
                .sorted((i, j)-> Integer.compare(a[j], a[i]))
                .mapToInt(ele->ele).toArray();
    }

    public static int[] returnSortedIndexes(double[] a, boolean ascending) {
        if (ascending) {
            return IntStream.range(0, a.length).boxed()
                    .sorted((i, j)-> Double.compare(a[i], a[j]))
                    .mapToInt(ele->ele).toArray();
        }
        return IntStream.range(0, a.length).boxed()
                .sorted((i, j)-> Double.compare(a[j], a[i]))
                .mapToInt(ele->ele).toArray();
    }

    /**
     * copies the input array into an array of size one more column dimension and in that added dimension,
     * stores the original indexes.
     * @param a array whose rows will be sorted by column index. a cannot be jagged, that is, all rows must be same length;
     * @param sortByCol sort the rows of the array by this column index
     * @return a sorted array of size [a.length][a[0].length + 1] where the last column has been added to hold the
     * original row indexes.
     */
    public static int[][] ascSort(int[][] a, int sortByCol) {
        int n = a.length;
        int m = a[0].length;
        int[][] out = new int[n][m + 1];
        for(int i = 0; i < n; ++i) {
            System.arraycopy(a[i], 0, out[i], 0, m);
            out[i][m] = i;
        }
        Arrays.sort(out, (a1, a2)->(a1[sortByCol] - a2[sortByCol]));
        return out;
    }

    /**
     * copies the input array into an array of size one more column dimension and in that added dimension,
     * stores the original indexes.
     * @param a array whose rows will be sorted by column index. a cannot be jagged, that is, all rows must be same length;
     * @param sortByCol sort the rows of the array by this column index
     * @return a sorted array of size [a.length][a[0].length + 1] where the last column has been added to hold the
     * original row indexes.
     */
    public static double[][] ascSort(double[][] a, int sortByCol) {
        int n = a.length;
        int m = a[0].length;
        double[][] out = new double[n][m + 1];
        for(int i = 0; i < n; ++i) {
            System.arraycopy(a[i], 0, out[i], 0, m);
            out[i][m] = i;
        }
        Arrays.sort(out, (a1, a2)->(Double.compare(a1[sortByCol], a2[sortByCol])));
        return out;
    }

    /**
     * copies the input array into an array of size one more column dimension and in that added dimension,
     * stores the original indexes.
     * @param a array whose rows will be sorted by column index. a cannot be jagged, that is, all rows must be same length;
     * @param sortByCol sort the rows of the array by this column index
     * @return a sorted array of size [a.length][a[0].length + 1] where the last column has been added to hold the
     * original row indexes.
     */
    public static float[][] ascSort(float[][] a, int sortByCol) {
        int n = a.length;
        int m = a[0].length;
        float[][] out = new float[n][m + 1];
        for(int i = 0; i < n; ++i) {
            System.arraycopy(a[i], 0, out[i], 0, m);
            out[i][m] = i;
        }
        Arrays.sort(out, (a1, a2)->(Float.compare(a1[sortByCol], a2[sortByCol])));
        return out;
    }

    /**
     * use merge sort to sort a1 in decreasing value order and return the 
     * indexes of the original a1 indexes in the sorted order.
     * 
     @param a1
     @return indexes 
     */
    public static int[] mergeSortDecreasing(double[] a1) {
        int[] indexes = new int[a1.length];
        for (int i = 0; i < indexes.length; ++i) {
            indexes[i] = i;
        }
     
        sortBy1stArg(a1, indexes, 0, a1.length - 1, false);
        
        return indexes;
    }

    public static class SortMetaData {
        int[] indexes = null;
        int numberOfInversions = -1;
    }

    public static SortMetaData mergeSortCountInversions(double[] a1) {
        int[] indexes = new int[a1.length];
        for (int i = 0; i < indexes.length; ++i) {
            indexes[i] = i;
        }

        SortMetaData d = new SortMetaData();
        d.indexes = indexes;
        d.numberOfInversions = sortBy1stArg(a1, indexes, 0, a1.length - 1, true);

        return d;
    }
    
    /**
     * use merge sort to sort a1 in increasing value order and return the 
     * indexes of the original a1 indexes in the sorted order.
     * 
     @param a1
     @return indexes 
     */
    public static int[] mergeSortIncreasing(double[] a1) {
        int[] indexes = new int[a1.length];
        for (int i = 0; i < indexes.length; ++i) {
            indexes[i] = i;
        }
     
        sortBy1stArg(a1, indexes, 0, a1.length - 1, true);
        
        return indexes;
    }

    /**
     * ascending sort array a by non-recursive merge-sort.
     * It has O(n*lg(n)) runtime complexity and O(n) space complexity where
     * n is a.length;
     * @param y array to be ascending sorted
     */
    public static void nonRecursiveMergeSort(double[] y) {
        int n = y.length;

        // adapted from https://www.baeldung.com/cs/non-recursive-merge-sort
        // and Chaudhuri & Wu 2018

        double[] t = new double[n];
        int idxT;
        int gap;
        int st1, e1, st2, e2;
        for (int len = 1; len < n; len*=2) {
            gap = 2*len;
            for (int j = 0; j < n; j+=gap) {
                st1 = j;
                e1 = Math.min(j + len - 1, n);
                st2 = j + len;
                e2 = Math.min(j + 2*len - 1, n-1);
                if (st2 >= n) {
                    break;
                }
                //if (e2 >= n) {
                //    e2 = n - 1;
                //}
                // temp length minimum size is 1, maximum size is n
                // -- begin merge.
                idxT = 0;
                while (st1 <= e1 && st2 <= e2) {
                    if (y[st1] <= y[st2]) {
                        t[idxT] = y[st1];
                        ++st1;
                    } else {
                        // inversion: a pair is a[st1], a[st2]
                        System.out.printf("inversion (%.3e,%.3e)\n", y[st1], y[st2]);
                        t[idxT] = y[st2];
                        ++st2;
                    }
                    ++idxT;
                }
                while (st1 <= e1) {
                    t[idxT] = y[st1];
                    ++st1;
                    ++idxT;
                }
                while (st2 <= e2) {
                    t[idxT] = y[st2];
                    ++st2;
                    ++idxT;
                }
                // -- end merge.
                for (int j2 = 0; j2 <= (e2 - st1 + 1); ++j2) {
                    y[j+j2] = t[j2];
                }
            }
        }
    }

    /**
     * calculate the order of indexes of y for ascending sort.  The method does not sort y, and returns
     * the indexes of y which put it in an acending sorted order.
     * The method has O(n*lg(n)) runtime complexity and O(n) space complexity where
     * n is y.length;
     * @param y array
     * @return the sorted index numbers of y
     */
    public static int[] nonRecursiveMergeSortIndexes(double[] y) {
        int n = y.length;

        // adapted from https://www.baeldung.com/cs/non-recursive-merge-sort
        // and Chaudhuri & Wu 2018

        //% The columns of idx are buffers to store sort indices and output buffer
        //%of merge-sort
        //%we alternate between them and avoid unnecessary copying
        int col = 0;
        int[][] idxN2 = new int[n][2];
        for (int i = 0; i < n; ++i) {
            idxN2[i] = new int[2];
            idxN2[i][col] = i;
        }
        int[] idxN = new int[n];

        int gap;
        int st1, e1, st2, e2;
        int idx1, idx2;
        int colS = 1;
        int k, kf;
        for (int len = 1; len < n; len*=2) {
            gap = 2*len;
            for (int z = 0; z < idxN2.length; ++z) {
                idxN[z] = idxN2[z][col];
            }
            k = 0;

            for (int j = 0; j < n; j+=gap) {
                st1 = j;
                e1 = Math.min(j + len - 1, n);
                st2 = j + len;
                e2 = Math.min(j + 2*len - 1, n-1);
                if (st2 >= n) {
                    break;
                }

                // -- begin merge.
                while (st1 <= e1 && st2 <= e2) {
                    k++;
                    idx1 = idxN[st1];
                    idx2 = idxN[st2];
                    if (y[idx1] <= y[idx2]) {
                        idxN2[k-1][colS] = idx1;
                        ++st1;
                    } else {
                        // an inversion pair is y[idx1], y[idx2]
                        //System.out.printf("inversion idxL1=%d idxL2=%d (%.3e,%.3e)\n", idx1, idx2, y[idx1], y[idx2]);
                        idxN2[k-1][colS] = idx2;
                        ++st2;
                    }
                }
                if (st1 <= e1) {
                    kf = k + e1 - st1 + 1;
                    int c = st1;
                    for (int z = k; z < kf; ++z) {
                        idxN2[z][colS] = idxN[c];
                        c++;
                    }
                    k = kf;
                } else if (st2 <= e2) {
                    kf = k + e2 - st2 + 1;
                    //idx( ( k+1):kf, s ) = idx_r( st2 : e2, : );
                    int c = st2;
                    for (int z = k; z < kf; ++z) {
                        idxN2[z][colS] = idxN[c];
                        c++;
                    }
                    k = kf;
                }
            }// end for j=

            col ^= 1;
            colS ^= 1;
        } // end for len=

        int lastCol = colS^1;

        System.out.printf("idxN2=\n%s\n", FormatArray.toString(idxN2, "%d"));
        System.out.printf("last col used = %d\n", lastCol);

        // file idxN and return it
        for (int i = 0; i < n; ++i) {
            idxN[i] = idxN2[i][lastCol];
        }
        return idxN;
    }
    /**
     * sort array a by the order given in the values of indexes.
     * @param a the array to be sorted
     * @param indexes array holding the indexes for the sort order
     */
    public static void sortByIndexes(double[] a, int[] indexes) {
        int n = a.length;
        if (indexes.length != a.length) {
            throw new IllegalArgumentException("a.length must equal indexes.length");
        }
        double[] s = new double[n];
        for (int i = 0; i < n; ++i) {
            s[i] = a[indexes[i]];
        }
        System.arraycopy(s, 0, a, 0, n);
    }
    /**
     * sort array a by the order given in the values of indexes.
     * @param a the array to be sorted
     * @param indexes array holding the indexes for the sort order
     */
    public static void sortByIndexes(float[] a, int[] indexes) {
        int n = a.length;
        if (indexes.length != a.length) {
            throw new IllegalArgumentException("a.length must equal indexes.length");
        }
        float[] s = new float[n];
        for (int i = 0; i < n; ++i) {
            s[i] = a[indexes[i]];
        }
        System.arraycopy(s, 0, a, 0, n);
    }

    /**
     * sort array a by the order given in the values of indexes.
     * @param a the array to be sorted
     * @param indexes array holding the indexes for the sort order
     * @return the sorted values of a
     */
    public static void sortByIndexes(int[] a, int[] indexes) {
        int n = a.length;
        if (indexes.length != a.length) {
            throw new IllegalArgumentException("a.length must equal indexes.length");
        }
        int[] s = new int[n];
        for (int i = 0; i < n; ++i) {
            s[i] = a[indexes[i]];
        }
        System.arraycopy(s, 0, a, 0, n);
    }

    /**
     * use merge sort to sort a1 and a2 by a1 in ascending order, breaking ties 
     * by a2 and return the indexes of the original indexes 
     * (i.e. 0 through a1.length-1) in the sorted order.
     * 
     @param a1
     @param a2
     @return indexes 
     */
    public static int[] mergeBy1stArgThen2nd(double[] a1, double[] a2) {
        
        if (a1.length != a2.length) {
            throw new IllegalArgumentException("a1 and a2 must be same length");
        }
        
        int[] indexes = new int[a1.length];
        for (int i = 0; i < indexes.length; ++i) {
            indexes[i] = i;
        }
        
        sortBy1stArgThen2nd(a1, a2, indexes, 0, a1.length - 1);
        
        return indexes;
    }
    
    /**
     * use quicksort to sort a by ascending values and
     * perform the same operations on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     @param a
     @param b
     */
    public static void sortBy1stArg(int[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
    
    
    /**
     * sort by increasing value a1 and apply same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArgThen2nd(float[] a1, float[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, 0, a1.length - 1);
    }
    
     
    /**
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param idxLo starting index of sorting of a1, inclusive
     @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArgThen2nd(float[] a1, float[] a2, int idxLo, 
        int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeBy1stArgThen2nd(float[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l == r) {
                float lx = a2Left[leftPos];
                float rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }

    /**
     *
     @param a
     @param b
     */
    public static void sortBy1stArg(int[] a, float[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg(a, b, 0, a.length - 1);
    }
       
    /**
     * use mergesort to sort by increasing values of a and apply same swaps to
     * b.
     @param a
     @param b 
     */
    public static void sortBy1stArg2(int[] a, int[] b) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException("a and b must be the same length");
        }
        
        sortBy1stArg2(a, b, 0, a.length - 1);
    }
    
    /**
     * use quicksort to 
       sort a from index idxLo to idxHi, inclusive, by ascending values and
     * perform the same operations on b.  Uses the optimized
     * qsort3 from the book "Programming in Pearls" by Jon Bentley.
     @param a
     @param b
     @param idxLo
     @param idxHi 
     */
    public static void sortBy1stArg(int[] a, int[] b, int idxLo, int idxHi) {
        
        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return;
        }
        if (idxLo < idxHi) {

            int x = a[idxLo];
            int store = idxLo;
            int idxMid = idxHi + 1;

            while (true) {
                do {
                    store++;     
                } while ((store <= idxHi) && (a[store] < x));
                do {
                    idxMid--;
                } while (a[idxMid] > x);
                if (store > idxMid) {
                    break;
                }
                int swap = a[store];
                a[store] = a[idxMid];
                a[idxMid] = swap;
                int swap2 = b[store];
                b[store] = b[idxMid];
                b[idxMid] = swap2;
            }
            int swap = a[idxLo];
            a[idxLo] = a[idxMid];
            a[idxMid] = swap;
            int swap2 = b[idxLo];
            b[idxLo] = b[idxMid];
            b[idxMid] = swap2;
         
            sortBy1stArg(a, b, idxLo, idxMid - 1);

            sortBy1stArg(a, b, idxMid + 1, idxHi);
        }
    }

    /**
     * use mergesort to sort by decreasing value a1 and apply 
       same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     */
    public static void sortByDecr(int[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
              
    }
    
    /**
     *
     @param a1
     @param a2
     */
    public static void sortBy1stArgDecrThen2ndIncr(int[] a1, int[] a2) {
        
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortBy1stArgDecrThen2ndIncr(a1, a2, 0, a1.length - 1);
    }

    /**
     * use mergesort to sort by decreasing value a1 and apply 
       same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     */
    public static void sortByDecr(float[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        
        sortByDecr(a1, a2, 0, a1.length - 1);
              
    }
    
/**
     * sort by increasing value a1 and apply same changes to a2.
     * 
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArg(float[] a1, float[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArg(a1, a2, 0, a1.length - 1);
    }
    
    /**
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param idxLo starting index of sorting of a1, inclusive
     @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArg(float[] a1, float[] a2, int idxLo, 
        int idxHi) {

        int indexMid = -1;

        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArg(a1, a2, idxLo, indexMid);
            sortBy1stArg(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeBy1stArg(float[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        a2Left[nLeft] = Float.MAX_VALUE;
        a1Left[nLeft] = Float.MAX_VALUE;
        a2Right[nRight] = Float.MAX_VALUE;
        a1Right[nRight] = Float.MAX_VALUE;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];

            if (l <= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    
    /**
     * sort by increasing value a1 and apply same changes to a2 using merge-sort.
     * Ties are further sorted by increasing values of a2.
     * The original indexes, sorted by the same swaps as a1 and a2 are returned.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @return the original indexes, sorted by the same swaps that a1 and a2 received.
     */
    public static int[] sortBy1stArg2(long[] a1, long[] a2) {
        if (a1.length != a2.length) {
            throw new IllegalArgumentException("a1 and a2 must be the same length");
        }
        int[] idxs = new int[a1.length];
        for (int i = 0; i < a1.length; ++i) {
            idxs[i] = i;
        }
        sortBy1stArg2(a1, a2, idxs, 0, a1.length - 1);
        return idxs;
    }
    
    /**
     * sort by increasing value a1 and apply same changes to a2 using merge-sort.
     * Ties are further sorted by increasing values of a2.
     * The original indexes, sorted by the same swaps as a1 and a2 are returned.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param idxs the original indexes, sorted by the same swaps that a1 and a2 received.
     @param iLo
     @param iHi
     */
    protected static void sortBy1stArg2(long[] a1, long[] a2, int[] idxs, int iLo, int iHi) {
        if (iLo < iHi) {
            int iMid = iLo + (int)((iHi - iLo)/2.);
            sortBy1stArg2(a1, a2, idxs, iLo, iMid);
            sortBy1stArg2(a1, a2, idxs, iMid + 1, iHi);
            mergeSortBy1stArg2(a1, a2, idxs, iLo, iMid, iHi);
        }
    }
        
    private static void mergeSortBy1stArg2(long[] a, long[] b, int[] c, 
        int iLo, int iMid, int iHi) {
        
        // 0-based: [iLo, iMid] inclusive
        //          [iMid+1, iHi] inclusive
                
        /*
        int nLeft = iMid - iLo + 1; 
        int nRight = iHi - iMid; 
        
        long[] aLeft = new long[nLeft + 1];
        long[] bLeft = new long[nLeft + 1];
        int[] cLeft = new int[nLeft + 1];
        long[] aRight = new long[nRight + 1];
        long[] bRight = new long[nRight + 1];
        int[] cRight = new int[nRight + 1];
        System.arraycopy(a, iLo, aLeft, 0, nLeft); // indexes inclusive [iLo, iLo + iMid - iLo + 1 - 1] = [iLo, iMid]
        System.arraycopy(b, iLo, bLeft, 0, nLeft);
        System.arraycopy(c, iLo, cLeft, 0, nLeft);
        
        System.arraycopy(a, iMid + 1, aRight, 0, nRight); // indexes inclusive [iMid + 1, iMid + 1 + iHi - iMid - 1] = [iMid + 1, iHi]
        System.arraycopy(b, iMid + 1, bRight, 0, nRight);
        System.arraycopy(c, iMid + 1, cRight, 0, nRight);
        */
        
        long[] aLeft = Arrays.copyOfRange(a, iLo, iMid + 2);
        long[] bLeft = Arrays.copyOfRange(b, iLo, iMid + 2);
        int[] cLeft = Arrays.copyOfRange(c, iLo, iMid + 2);
        
        long[] aRight = Arrays.copyOfRange(a, iMid + 1, iHi + 2);
        long[] bRight = Arrays.copyOfRange(b, iMid + 1, iHi + 2);
        int[] cRight = Arrays.copyOfRange(c, iMid + 1, iHi + 2);
        
        aLeft[aLeft.length - 1] = Long.MAX_VALUE;
        bLeft[bLeft.length - 1] = Long.MAX_VALUE;
        cLeft[cLeft.length - 1] = Integer.MAX_VALUE;
        aRight[aRight.length - 1] = Long.MAX_VALUE;
        bRight[bRight.length - 1] = Long.MAX_VALUE;
        cRight[cRight.length - 1] = Integer.MAX_VALUE;
        
        int iL = 0;
        int iR = 0;
        
        long l;
        long r;
        
        for (int k = iLo; k <= iHi; ++k) {
            l = aLeft[iL];
            r = aRight[iR];
            if (l < r) {
                a[k] = l;
                b[k] = bLeft[iL];
                c[k] = cLeft[iL];
                iL++;
            } else if (l > r) {
                a[k] = r;
                b[k] = bRight[iR];
                c[k] = cRight[iR];
                iR++;
            } else {
                // a1[iL] == a1[iR] so break ties using a2
                l = bLeft[iL];
                r = bRight[iR];
                if (l <= r) {
                    a[k] = aLeft[iL];
                    b[k] = bLeft[iL];
                    c[k] = cLeft[iL];
                    iL++;
                } else {
                    a[k] = aRight[iR];
                    b[k] = bRight[iR];
                    c[k] = cRight[iR];
                    iR++;
                }
            }
        }
    }
        
    /**
     * use mergesort to sort a1 by increasing values and apply same changes to
     * a2.
     @param a1
     @param a2
     @param idxLo
     @param idxHi 
     */
    public static void sortBy1stArg2(int[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortBy1stArg2(a1, a2, idxLo, indexMid);
            
            sortBy1stArg2(a1, a2, indexMid + 1, idxHi);
            
            mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     *
     @param a1
     @param a2
     @param idxLo
     @param idxHi
     */
    public static void sortBy1stArg(int[] a1, float[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortBy1stArg(a1, a2, idxLo, indexMid);
            
            sortBy1stArg(a1, a2, indexMid + 1, idxHi);
            
            mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     *
     @param a1
     @param a2
     @param idxLo
     @param idxHi
     */
    public static void sortBy1stArgDecrThen2ndIncr(int[] a1, int[] a2, 
        int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortBy1stArgDecrThen2ndIncr(a1, a2, idxLo, indexMid);
            
            sortBy1stArgDecrThen2ndIncr(a1, a2, indexMid + 1, idxHi);
            
            mergeBy1stArgDecrThen2ndIncr(a1, a2, idxLo, indexMid, idxHi);
        }
    }

    /**
     * use mergesort to sort by decreasing value a1 and apply 
       same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
       @param idxLo smallest index to participate in sort
       @param idxHi largest index to participate n sort
    */
    public static void sortByDecr(int[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }

    /**
     * use mergesort to sort by decreasing value a1 and apply 
       same changes to a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
       @param idxLo smallest index to participate in sort
       @param idxHi largest index to participate n sort
    */
    public static void sortByDecr(float[] a1, int[] a2, int idxLo, int idxHi) {

        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortByDecr(a1, a2, idxLo, indexMid);
            
            sortByDecr(a1, a2, indexMid + 1, idxHi);
            
            mergeByDecr(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeByDecr(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArg(int[] a1, float[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        float[] a2Left = new float[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        float[] a2Right = new float[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        float sentinel2 = Float.POSITIVE_INFINITY;
        a2Left[nLeft] = sentinel2;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel2;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l <= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArg(int[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l <= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeBy1stArgDecrThen2ndIncr(int[] a1, int[] a2, 
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;
        
        int leftPos = 0;
        int rightPos = 0;
        boolean lft = false;
        
        for (int k = idxLo; k <= idxHi; k++) {
            int l = a1Left[leftPos];
            int r = a1Right[rightPos];
            if (l > r) {
                lft = true;
            } else if (l == r) {
                int l2 = a2Left[leftPos];
                int r2 = a2Right[rightPos];
                if (l2 <= r2) {
                    lft = true;
                } else {
                    lft = false;
                }
            } else {
                lft = false;
            }
            
            if (lft) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    private static void mergeByDecr(float[] a1, int[] a2, int idxLo, 
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        float[] a1Left = new float[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        float[] a1Right = new float[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        float sentinel1 = Float.NEGATIVE_INFINITY;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel1;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel1;
        
        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            float l = a1Left[leftPos];
            float r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
    
    /**
     * sort by increasing value a1 for ties sort by a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        sortBy1stArgThen2nd(a1, a2, 0, a1.length - 1);
    }
    
    /**
     * sort by increasing value a1 for ties sort by a2.
     * Ties are further sorted by increasing values of a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param a3 array of points sorted by same order as a1 and a2, useful to pass in
     * the original index ordering for example.
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2, int[] a3) {
        if (a1 == null) {
            throw new IllegalArgumentException("a1 cannot be null");
        }
        if (a2 == null) {
            throw new IllegalArgumentException("a2 cannot be null");
        }
        if (a3 == null) {
            throw new IllegalArgumentException("a3 cannot be null");
        }
        if (a1.length != a2.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a2");
        }
        if (a1.length != a3.length) {
            throw new IllegalArgumentException(
            "number of items in a1 must be the same as in a3");
        }
        sortBy1stArgThen2nd(a1, a2, a3, 0, a1.length - 1);
    }
    
    /**
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also and break ties in a1
     @param indexes array of original indexes of a1 and a2, progressively sorted
     * by same swaps.
     @param idxLo starting index of sorting of a1, inclusive
     @param idxHi stopping index of sorting of a1, inclusive
     */
    private static void sortBy1stArgThen2nd(double[] a1, double[] a2, int[] indexes, int idxLo,
        int idxHi) {

        int indexMid = -1;
        
        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, indexes, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexes, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, indexes, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param idxLo starting index of sorting of a1, inclusive
     @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2, int idxLo,
        int idxHi) {

        int indexMid = -1;
        
        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, idxLo, indexMid, idxHi);
        }
    }
    
    /**
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also and act as tie breaker
     @param a3 array of points sorted by same order as a1 and a2, useful to pass in
     * the original index ordering for example.
     @param idxLo starting index of sorting of a1, inclusive
     @param idxHi stopping index of sorting of a1, inclusive
     */
    public static void sortBy1stArgThen2nd(int[] a1, int[] a2, int[] a3, int idxLo,
        int idxHi) {

        int indexMid = -1;
        
        if (idxLo < idxHi) {

            indexMid = (idxLo + idxHi) >> 1;
            sortBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid);
            sortBy1stArgThen2nd(a1, a2, a3, indexMid + 1, idxHi);
            mergeBy1stArgThen2nd(a1, a2, a3, idxLo, indexMid, idxHi);
        }
    }
    
    private static void mergeBy1stArgThen2nd(int[] a1, int[] a2, int idxLo,
        int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);

        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);

        int sentinel = Integer.MAX_VALUE;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;
        int l, r;

        for (int k = idxLo; k <= idxHi; k++) {
            l = a1Left[leftPos];
            r = a1Right[rightPos];

            if (l == r) {
                l = a2Left[leftPos];
                r = a2Right[rightPos];

                if (l <= r) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
 
    private static void mergeBy1stArgThen2nd(int[] a1, int[] a2, int[] a3, 
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a3Left = new int[nLeft + 1];
        int[] a2Left = new int[nLeft + 1];
        int[] a1Left = new int[nLeft + 1];

        int[] a3Right = new int[nRight + 1];
        int[] a2Right = new int[nRight + 1];
        int[] a1Right = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        System.arraycopy(a3, idxLo, a3Left, 0, nLeft);

        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        System.arraycopy(a3, idxMid + 1, a3Right, 0, nRight);

        int sentinel = Integer.MAX_VALUE;
        a3Left[nLeft] = sentinel;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a3Right[nRight] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;
        int l, r;

        for (int k = idxLo; k <= idxHi; k++) {
            l = a1Left[leftPos];
            r = a1Right[rightPos];

            if (l == r) {
                l = a2Left[leftPos];
                r = a2Right[rightPos];

                if (l <= r) {
                    a2[k] = a2Left[leftPos];
                    a3[k] = a3Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    leftPos++;
                } else {
                    a3[k] = a3Right[rightPos];
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a3[k] = a3Left[leftPos];
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a3[k] = a3Right[rightPos];
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
            }
        }
    }
 
    private static void mergeBy1stArgThen2nd(double[] a1, double[] a2, 
        int[] indexes, int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        double[] a2Left = new double[nLeft + 1];
        double[] a1Left = new double[nLeft + 1];

        double[] a2Right = new double[nRight + 1];
        double[] a1Right = new double[nRight + 1];
        
        int[] indexesLeft = new int[nLeft + 1];
        int[] indexesRight = new int[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);

        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        System.arraycopy(indexes, idxLo, indexesLeft, 0, nLeft);
        System.arraycopy(indexes, idxMid + 1, indexesRight, 0, nRight);

        double sentinel = Double.POSITIVE_INFINITY;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel;

        int leftPos = 0;
        int rightPos = 0;

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];

            if (l == r) {
                double lx = a2Left[leftPos];
                double rx = a2Right[rightPos];

                if (lx <= rx) {
                    a2[k] = a2Left[leftPos];
                    a1[k] = a1Left[leftPos];
                    indexes[k] = indexesLeft[leftPos];
                    leftPos++;
                } else {
                    a2[k] = a2Right[rightPos];
                    a1[k] = a1Right[rightPos];
                    indexes[k] = indexesRight[rightPos];
                    rightPos++;
                }
            } else if (l < r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                indexes[k] = indexesLeft[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                indexes[k] = indexesRight[rightPos];
                rightPos++;
            }
        }
    }

    /**
     * use mergesort to sort by decreasing value a1 and apply 
       same changes to a2.
     * runtime is O(N * log_2(N))
     *
     @param a1 array of points to be sorted
     @param a2 array of points to apply a1 sorting to also
     @param idxLo
     @param idxHi 
     @param ascendingSort
     @return number of inversions
     */
    private static int sortBy1stArg(double[] a1, int[] a2, int idxLo,
        int idxHi, boolean ascendingSort) {

        int nInv = 0;
        
        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;

            nInv += sortBy1stArg(a1, a2, idxLo, indexMid, ascendingSort);

            nInv += sortBy1stArg(a1, a2, indexMid + 1, idxHi, ascendingSort);

            nInv += mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi, ascendingSort);
        }
        return nInv;
    }
    
    private static int mergeBy1stArg(double[] a1, int[] a2,
        int idxLo, int idxMid, int idxHi, boolean ascendingSort) {
        if (ascendingSort) {
            return mergeBy1stArgIncr(a1, a2, idxLo, idxMid, idxHi);
        } else {
            return mergeBy1stArgDecr(a1, a2, idxLo, idxMid, idxHi);
        }
    }

    private static int mergeBy1stArgDecr(double[] a1, int[] a2,
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        double[] a1Left = new double[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        double[] a1Right = new double[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MIN_VALUE;
        double sentinel1 = Double.NEGATIVE_INFINITY;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel1;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel1;
        
        int leftPos = 0;
        int rightPos = 0;

        int nInv = 0;
        int leftLen = a1Left.length - 1;

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];
            if (l >= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;
                nInv += (leftLen - leftPos);
            }
        }
        return nInv;
    }
    
    private static int mergeBy1stArgIncr(double[] a1, int[] a2,
        int idxLo, int idxMid, int idxHi) {

        int nLeft = idxMid - idxLo + 1;
        int nRight = idxHi - idxMid;

        int[] a2Left = new int[nLeft + 1];
        double[] a1Left = new double[nLeft + 1];

        int[] a2Right = new int[nRight + 1];
        double[] a1Right = new double[nRight + 1];

        System.arraycopy(a1, idxLo, a1Left, 0, nLeft);
        System.arraycopy(a2, idxLo, a2Left, 0, nLeft);
        
        System.arraycopy(a1, idxMid + 1, a1Right, 0, nRight);
        System.arraycopy(a2, idxMid + 1, a2Right, 0, nRight);
        
        int sentinel = Integer.MAX_VALUE;
        double sentinel1 = Double.POSITIVE_INFINITY;
        a2Left[nLeft] = sentinel;
        a1Left[nLeft] = sentinel1;
        a2Right[nRight] = sentinel;
        a1Right[nRight] = sentinel1;
        
        int leftPos = 0;
        int rightPos = 0;

        int nInv = 0;
        int leftLen = a1Left.length - 1;

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];
            if (l <= r) {
                a2[k] = a2Left[leftPos];
                a1[k] = a1Left[leftPos];
                leftPos++;
            } else {
                //if (Math.abs(a1[k] - r) > 1e-7) System.out.printf("inversion pair (%.3e,%.3e)\n", a1[k], r);
                a2[k] = a2Right[rightPos];
                a1[k] = a1Right[rightPos];
                rightPos++;

                nInv += (leftLen - leftPos);
            }
        }
        return nInv;
    }

    /**
     * sort the given points in place by polar angle in counterclockwise order
     * around the first point x[0], y[0].
     * NOTE that the angles are rounded to degree integers for the "reduce to
     * unique" step.
     * The runtime complexity is O(n) where n = x.length.
     * @param x array of x points
     * @param y array of y points of same length as x.
     @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
      *                             from the first point in (x,y) will be kept. if false, all points are kept.
      * for reference, 0.1 degrees = 0.001745 radians.
      * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPointLinear(long[] x, long[] y, boolean reduceToUniqueAngles) {
        return sortCCWBy1stPointLinear(x, y, reduceToUniqueAngles, 1E-7);
    }

    /**
     * sort the given points in place by polar angle in counterclockwise order
     * around the first point x[0], y[0].
     * NOTE that the angles are rounded to degree integers for the "reduce to
     * unique" step.
     * The runtime complexity is O(n) where n = x.length.
     * @param x array of x points
     * @param y array of y points of same length as x.
     @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
      *                             from the first point in (x,y) will be kept. if false, all points are kept.
     * for reference, 0.1 degrees = 0.001745 radians.
     * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPointLinear(double[] x, double[] y, boolean reduceToUniqueAngles) {
        return sortCCWBy1stPointLinear(x, y, reduceToUniqueAngles, 1E-7);
    }

    /**
     * sort the given points in place by polar angle in counterclockwise order
     * around the first point x[0], y[0].
     * NOTE that the angles are rounded to degree integers for the "reduce to
     * unique" step.
     * The runtime complexity is O(n) where n = x.length.
     * @param x array of x points
     * @param y array of y points of same length as x.
     @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
      *                             from the first point in (x,y) will be kept. if false, all points are kept.
      * @param tolerance if reduceToUniqueAngles is true, tolerance is used to determing if two angles are equivalent.
     *                  The angles are in radians and so is the tolerance.
     * for reference, 0.1 degrees = 0.001745 radians.
     * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPointLinear(double[] x, double[] y, boolean reduceToUniqueAngles, double tolerance) {

        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }

        if (x.length < 3) {
            return x.length;
        }

        double x0 = x[0];
        double y0 = y[0];

        double[] x1 = Arrays.copyOfRange(x, 1, x.length);
        double[] y1 = Arrays.copyOfRange(y, 1, y.length);

        double degRadians;
        int[] polarAngleDegree1 = new int[x.length - 1];

        for (int i = 0; i < x1.length; i++) {
            degRadians = AngleUtil.polarAngleCCW((x1[i] - x0), (y1[i] - y0));
            polarAngleDegree1[i] = (int)Math.round(degRadians * (180./Math.PI));
        }

        int[] sortedIndexes1;

        int rangeOfData = MiscMath0.findMax(polarAngleDegree1) - MiscMath0.findMin(polarAngleDegree1);

        if ((rangeOfData < 360) && (polarAngleDegree1.length >= 360)) {
            // this will be O(360) if range of data is 0 to 359, inclusive
            sortedIndexes1 = CountingSort.sortAndReturnIndexes(polarAngleDegree1);
        } else {
            // this will be O(N) where N is length of input array
            sortedIndexes1 = BucketSort.sortAndReturnIndexes(polarAngleDegree1);
        }

        int i;
        double[] polarAngleDegree = new double[x.length];
        for (i = 0; i < x1.length; ++i) {
            x[i+1] = x1[sortedIndexes1[i]];
            y[i+1] = y1[sortedIndexes1[i]];
            polarAngleDegree[i+1] = polarAngleDegree1[i];
        }

        // for same polar angles, keep the one which is furthest from x0, p0.
        // assuming an angular resolution of 1 degree and using rounded integers for the angle degrees
        int nUsable;
        if (reduceToUniqueAngles) {
            nUsable = reduceToUniquePolarAngles(x, y, polarAngleDegree, tolerance);
        } else {
            nUsable = x.length;
        }

        return nUsable;
    }

    private static void rewriteBySortedIndexes(int[] sortedIndexes1, double[] x) {

        int n = x.length;

        double[] x2 = new double[n];
        for (int i = 0; i < n; ++i) {
            x2[i] = x[sortedIndexes1[i]];
        }

        System.arraycopy(x2, 0, x, 0, n);
    }


    /**
     * traverse the polar angles in pA, and if points have the same polar angles
     * within tolerance, only keep the one furthest from (x0, y0).
     * NOTE that x, y, deg should have been sorted by polar angle before using this method.
     * This method compacts arrays x, y, and pA are compacted so that the usable values are at the
     * top r indexes where r is returned by this method
     * @param x x coordinates of points
     * @param y y coordinates of points
     * @param polarAngle polar angles of points in x,y w.r.t to their first points x[0], y[0].
     *            x,y,deg should have been sorted by polar angle before using this method.
     * @return returns the number of indexes usable in each of x, y, and pA
     * after compacting the arrays to remove redundant polar angle degrees.
     */
    static int reduceToUniquePolarAngles(double[] x, double[] y, double[] polarAngle, double tolerance) {

        if (x.length != polarAngle.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and polarAngle must be same lengths");
        }

        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist;
        double dist;
        int nextI;
        int i2 = 1;
        int i = 1;
        assert(polarAngle[0] == 0);

        int n = x.length;

        // process any angle=0's first.  we want to keep the first point and not use distance argument on it
        while ((i < n) &&  (polarAngle[i] <= tolerance)) {
            ++i;
        }

        for (; i < x.length; i++) {
            // look ahead
            nextI = i + 1;
            iMaxDist = i;

            if ( (nextI < x.length)   && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= tolerance) ) {
                maxDist = relativeLengthOfLine(x[0], y[0], x[i], y[i]);
            }

            while ( (nextI < x.length)  && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= tolerance) ) {
                dist = relativeLengthOfLine(x[0], y[0], x[nextI], y[nextI]);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }

            x[i2] = x[iMaxDist];
            y[i2] = y[iMaxDist];
            polarAngle[i2] = polarAngle[iMaxDist];
            i = nextI - 1;
            ++i2;
        }

        return i2;
    }

    static double relativeLengthOfLine(double x1, double y1, double x2, double y2) {
        double dx2 = (x2 - x1);
        dx2 *= dx2;
        double dy2 = (y2 - y1);
        dy2 *= dy2;
        //double d = Math.sqrt(dx2 + dy2);
        return dx2 + dy2;
    }

    /**
     * sort the given points in place by polar angle in counterclockwise order
     * around the first point x[0], y[0].
     * NOTE that the angles are rounded to degree integers for the linear runtime.
     * The runtime complexity is O(n) where n = x.length.
     * @param x array of x points
     * @param y array of y points of same length as x.
     @param reduceToUniqueAngles if true, for any points having the same polar angle, only the point furthest
      *                             from the first point in (x,y) will be kept. if false, all points are kept.
      * @param tolerance if reduceToUniqueAngles is true, tolerance is used to determing if two angles are equivalent.
     *                  The angles are in radians and so is the tolerance.
     * for reference, 0.1 degrees = 0.001745 radians.
     * @return the number of usable points in the sorted arrays. if reduceToUniqueAngles is true,
     * then for points having same angle,
     * only the furthest from x[0], y[0] is kept, and so there may be points at the ends of arrays
     * x,y, and outPolarAngle that are unusable.
     */
    public static int sortCCWBy1stPointLinear(long[] x, long[] y, boolean reduceToUniqueAngles, double tolerance) {

        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }

        if (x.length < 3) {
            return x.length;
        }

        long x0 = x[0];
        long y0 = y[0];

        long[] x1 = Arrays.copyOfRange(x, 1, x.length);
        long[] y1 = Arrays.copyOfRange(y, 1, y.length);

        double degRadians;
        int[] polarAngleDegree1 = new int[x.length - 1];

        for (int i = 0; i < x1.length; i++) {
            degRadians = AngleUtil.polarAngleCCW((double)(x1[i] - x0), (double)(y1[i] - y0));
            polarAngleDegree1[i] = (int)Math.round(degRadians * (180./Math.PI));
        }

        int[] sortedIndexes1;

        int rangeOfData = MiscMath0.findMax(polarAngleDegree1) - MiscMath0.findMin(polarAngleDegree1);

        if ((rangeOfData < 360) && (polarAngleDegree1.length >= 360)) {
            // this will be O(360) if range of data is 0 to 359, inclusive
            sortedIndexes1 = CountingSort.sortAndReturnIndexes(polarAngleDegree1);
        } else {
            // this will be O(N) where N is length of input array
            sortedIndexes1 = BucketSort.sortAndReturnIndexes(polarAngleDegree1);
        }

        int i;
        double[] polarAngleDegree = new double[x.length];
        for (i = 0; i < x1.length; ++i) {
            x[i+1] = x1[sortedIndexes1[i]];
            y[i+1] = y1[sortedIndexes1[i]];
            polarAngleDegree[i+1] = polarAngleDegree1[i];
        }

        // for same polar angles, keep the one which is furthest from x0, p0.
        // assuming an angular resolution of 1 degree and using rounded integers for the angle degrees
        int nUsable;
        if (reduceToUniqueAngles) {
            nUsable = reduceToUniquePolarAngles(x, y, polarAngleDegree, tolerance);
        } else {
            nUsable = x.length;
        }

        return nUsable;
    }

    private static void rewriteBySortedIndexes(int[] sortedIndexes1, long[] x) {

        int n = x.length;

        long[] x2 = new long[n];
        for (int i = 0; i < n; ++i) {
            x2[i] = x[sortedIndexes1[i]];
        }

        System.arraycopy(x2, 0, x, 0, n);
    }


    /**
     * traverse the polar angles in pA, and if points have the same polar angles
     * within tolerance, only keep the one furthest from (x0, y0).
     * NOTE that x, y, deg should have been sorted by polar angle before using this method.
     * This method compacts arrays x, y, and pA are compacted so that the usable values are at the
     * top r indexes where r is returned by this method
     * @param x x coordinates of points
     * @param y y coordinates of points
     * @param polarAngle polar angles of points in x,y w.r.t to their first points x[0], y[0].
     *            x,y,deg should have been sorted by polar angle before using this method.
     * @return returns the number of indexes usable in each of x, y, and pA
     * after compacting the arrays to remove redundant polar angle degrees.
     */
    static int reduceToUniquePolarAngles(long[] x, long[] y, double[] polarAngle, double tolerance) {

        if (x.length != polarAngle.length || x.length != y.length) {
            throw new IllegalArgumentException("x, y, and polarAngle must be same lengths");
        }

        double maxDist = Double.NEGATIVE_INFINITY;
        int iMaxDist;
        double dist;
        int nextI;
        int i2 = 1;
        int i = 1;
        assert(polarAngle[0] == 0);

        int n = x.length;

        // process any angle=0's first.  we want to keep the first point and not use distance argument on it
        while ((i < n) &&  (polarAngle[i] <= tolerance)) {
            ++i;
        }

        for (; i < x.length; i++) {

            // look ahead
            nextI = i + 1;
            iMaxDist = i;

            if ( (nextI < x.length)   && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= tolerance) ) {
                maxDist = relativeLengthOfLine(x[0], y[0], x[i], y[i]);
            }

            while ( (nextI < x.length)  && (Math.abs( polarAngle[i] - polarAngle[nextI] ) <= tolerance) ) {
                dist = relativeLengthOfLine(x[0], y[0], x[nextI], y[nextI]);
                if (maxDist < dist) {
                    maxDist = dist;
                    iMaxDist = nextI;
                }
                nextI++;
            }

            x[i2] = x[iMaxDist];
            y[i2] = y[iMaxDist];
            polarAngle[i2] = polarAngle[iMaxDist];
            i = nextI - 1;
            ++i2;
        }

        return i2;
    }
}
