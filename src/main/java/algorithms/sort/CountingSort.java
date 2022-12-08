package algorithms.sort;

import java.util.Arrays;

/**
 * a sort for integers.  runtime complexity Math.max(n, max(a) - min(a))
 * <pre>
 * To use this algorithm: the range of the numbers in the array should probably not 
 * be much greater than 1e7 unless the jvm settings for maximum stack size
 * are increased.  An internal long array of size maximum of array values
 * is constructed and that consumes memory which also affects
 * performance for max &gt; 11e7.
 * 
 * implemented from Cormen, Leiserson, Rivest, and Stein "Introduction to Algorithms"

   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * </pre>
 * @author nichole
 */
public class CountingSort {
    
    /**
     * sort the members of a
     * <pre>
     * runtime complexity: Math.max(n, max(a) - min(a))
     * </pre>
     @param a
     @return  
     */
    public static int[] sort(int[] a) {

        int[] b = Arrays.copyOf(a, a.length);
        int[] idxs = sortAndReturnIndexes(b);
        
        return b;
    }
    
    /**
     * sort the members of a and return the original indexes.
     * <pre>
     * runtime complexity: Math.max(n, max(a) - min(a))
     * </pre>
     @param a input and output array a
     @return the original indexes of a, sorted
     */
    public static int[] sortAndReturnIndexes(final int[] a) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (a.length < 2) {
            return Arrays.copyOf(a, a.length);
        }
        
        // shift range of a to non-negative numbers, starting at 0.
        int max = Integer.MIN_VALUE;
        int min = Integer.MAX_VALUE;
        int i;
        for (i = 0; i < a.length; ++i) {
            if (a[i] > max) {
                max = a[i];
            }
            if (a[i] < min) {
                min = a[i];
            }
        }
        int[] a2 = Arrays.copyOf(a, a.length);
        for (i = 0; i < a2.length; ++i) {
           a2[i] -= min;
        }
        
        long[] c = new long[max - min + 1];

        // c holds frequency of each number by index, e.g. c[0] holds the number of 0's in a
        for (i = 0; i < a2.length; i++) {
            c[a2[i]]++;
        }
                
        // cumulative sum to end of array c.  the last item in c holds the 
        // total number of items in 'a' less than or equal to max
        for (i = 1; i < c.length; i++) {
            c[i] += c[i - 1];
        }
        
        int[] b = new int[a.length];
        
        int[] idxs = new int[a.length];
        Arrays.fill(idxs, -1);
                
        // use the order imposed by c to write the values of a into b.  c holds
        // frequency, that is place markers too so that is updated as b is written
        int aI;
        for (i = (a2.length - 1); i > -1; i--) {
            aI = a2[i];
            c[aI]--;
            b[(int)c[aI]] = aI + min;
            idxs[(int)c[aI]] = i;
        }
        
        System.arraycopy(b, 0, a, 0, a.length);
        
        return idxs;
    }
    
    /**
     * sort the members of a and apply the same
     * changes of item position to b.
     * runtime complexity: Math.max(n, max(a) - min(a))
     * 
     @param a
     @param b
     */
    public static void sort(final int[] a, final int[] b) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "the lengths of a and b must be the same");
        }
        if (a.length < 2) {
            return;
        }
        
        int[] a2 = Arrays.copyOf(a, a.length);
        int[] idxs = sortAndReturnIndexes(a2);
        
        int[] b2 = new int[a.length];
        for (int i = 0; i < b.length; ++i) {
            b2[i] = b[idxs[i]];
        }

        System.arraycopy(a2, 0, a, 0, a.length);
        System.arraycopy(b2, 0, b, 0, b.length);
    }
    
    /**
     * apply a descending sort to the members of a
     * and apply the same changes of item position to b.
     * runtime complexity: Math.max(n, max(a) - min(a))
     * 
     @param a
     @param b
     */
    public static void sortByDecr(final int[] a, final int[] b) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (b == null) {
            throw new IllegalArgumentException("b cannot be null");
        }
        if (a.length != b.length) {
            throw new IllegalArgumentException(
            "the lengths of a and b must be the same");
        }
        if (a.length < 2) {
            return;
        }
        
        int[] a2 = Arrays.copyOf(a, a.length);
        int[] idxs = sortAndReturnIndexes(a2);        
        int[] b2 = new int[a.length];
        int i;
        for (i = 0; i < b.length; ++i) {
            b2[i] = b[idxs[i]];
        }
        
        // reverse
        int n2 = a2.length >> 1;
        int i2 = a2.length - 1;
        int swap;
        for (i = 0; i < n2; ++i) {
            swap = a2[i];
            a2[i] = a2[i2];
            a2[i2] = swap;
            swap = b2[i];
            b2[i] = b2[i2];
            b2[i2] = swap;
            i2--;
        }

        System.arraycopy(a2, 0, a, 0, a.length);
        System.arraycopy(b2, 0, b, 0, b.length);
    }

}
