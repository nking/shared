package algorithms.sort;

import java.util.Arrays;

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
     */
    private static void sortBy1stArg(double[] a1, int[] a2, int idxLo, 
        int idxHi, boolean ascendingSort) {
        
        if (idxLo < idxHi) {

            int indexMid = (idxLo + idxHi) >> 1;
            
            sortBy1stArg(a1, a2, idxLo, indexMid, ascendingSort);
            
            sortBy1stArg(a1, a2, indexMid + 1, idxHi, ascendingSort);
            
            mergeBy1stArg(a1, a2, idxLo, indexMid, idxHi, ascendingSort);
        }
    }
    
    private static void mergeBy1stArg(double[] a1, int[] a2, 
        int idxLo, int idxMid, int idxHi, boolean ascendingSort) {
        if (ascendingSort) {
            mergeBy1stArgIncr(a1, a2, idxLo, idxMid, idxHi);
        } else {
            mergeBy1stArgDecr(a1, a2, idxLo, idxMid, idxHi);
        }
    }

    private static void mergeBy1stArgDecr(double[] a1, int[] a2, 
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
            }
        }
    }
    
    private static void mergeBy1stArgIncr(double[] a1, int[] a2, 
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

        for (int k = idxLo; k <= idxHi; k++) {
            double l = a1Left[leftPos];
            double r = a1Right[rightPos];
            if (l < r) {
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
}
