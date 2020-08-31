package algorithms.misc;

import java.util.Random;

/**
 * class to randomly select from the CDF of a discrete probability function.
 * 
 * @author nichole
 */
public class CDFRandomSelect {
    
    /**
     * choose k indexes from the cdf by randomly drawing from rand[0,1] and
     * using binary search to find the nearest cumulative probability in the cdf,
     * k times, returning the indexes.
     * 
     ** @param cdf cumulative distribution function
     * @param k number of random selects
     * @param rand random number generator
     * @return indexes from randomly selected cdf distribution values
     */
    public static int[] chooseKFromBinarySearch(double[] cdf, int k, Random rand) {
        return chooseKFromBinarySearch(cdf, k, rand, 1.e-15);
    }
    
    /**
     * choose k indexes from the cdf by randomly drawing from rand[0,1] and
     * using binary search to find the nearest cumulative probability in the cdf,
     * k times, returning the indexes.
     * 
     * @param cdf cumulative distribution function
     * @param k number of random selects
     * @param rand random number generator
     * @param tolerance tolerance for equality of doubles
     * @return indexes from randomly selected cdf distribution values
     */
    public static int[] chooseKFromBinarySearch(double[] cdf, int k, Random rand,
        double tolerance) {
        
        int[] selected = new int[k];
        
        for (int i = 0; i < k; ++i) {
            selected[i] = binarySearchForNearest(cdf, rand.nextDouble(), tolerance);
        }
        
        return selected;
    }
 
    private static int binarySearchForNearest(double[] cdf, final double srch, 
        final double tol) {
        
        int n = cdf.length;
           
        int lowIdx = 0;
        int highIdx = n - 1;
        int midIdx = (highIdx + lowIdx) >> 1;
        
        double v;
        int comp;
        
        while (lowIdx != highIdx) {

            midIdx = (highIdx + lowIdx) >> 1;
            
            v = cdf[midIdx];
            
            //-1, 0 or 1 when v is less than, equal to, or greater than value.
            comp = (Math.abs(v - srch) < tol) ? 0 : (v < srch) ? -1 : 1;

            if (comp > 0) {

                if (highIdx == midIdx) {
                    highIdx--;
                } else {
                    highIdx = midIdx;
                }

            } else if (comp < 0) {
                if (lowIdx == midIdx) {
                    lowIdx++;
                } else {
                    lowIdx = midIdx;
                }

            } else {
                // is equal
                return midIdx;
            }
        }

        //compare difference of midIdx to predecessor and successor
        int minIdx = midIdx;
        double diff = Math.abs(cdf[midIdx] - srch);
        
        if ((midIdx + 1) < n) {
            if (Math.abs(cdf[midIdx + 1] - srch) < diff) {
                minIdx = midIdx + 1;
                diff = Math.abs(cdf[midIdx] - srch);
            }
        }
        if ((midIdx - 1) >= 0) {
            if ( Math.abs(cdf[midIdx - 1] - srch) < diff) {
                minIdx = midIdx - 1;
            }
        }
        
        return minIdx;
    }
            
    /**
     * k indexes from randomly selected cdf distribution values.
     * Internally, the CDF is transformed from double values to integers
     * and stored in a YFastTrie.
     * <b>This method is a work in progress.</b>
     * The method is competitive when cdf.length is .gte. 2^15 limited by successor
     * method runtime roughly, and may be competitive for far smaller integer range
     * (depends upon: largest value to store in trie, how filled the trie is,
     * and whether the transformation to integers preserves enough significant 
     * digits for the distribution and for random numbers).
     * 
     * @param cdf cumulative distribution function
     * @param k number of random selects
     * @param rand random number generator
     * @return indexes from randomly selected cdf distribution values
     */
    public static int[] chooseKFromIntegerTransformAndTrie(double[] cdf, int k, Random rand) {       
        throw new UnsupportedOperationException("Not implemented yet."); 
    }

}
