package algorithms.statistics;

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
     *@param cdf cumulative distribution function.  NOTE that the values should
     * be normalized such that the last item is 1.
     @param k number of random selects
     @param rand random number generator
     @return indexes from randomly selected cdf distribution values
     */
    public static int[] chooseKFromBinarySearch(double[] cdf, int k, Random rand) {
        return chooseKFromBinarySearch(cdf, k, rand, 1.e-15);
    }

    /**
     * choose k indexes from the cdf by randomly drawing from rand[0,1] and
     * using binary search to find the nearest cumulative probability in the cdf,
     * k times, returning the indexes.
     * 
     @param cdf cumulative distribution function.  NOTE that the values should
     * be normalized such that the last item is 1.
     @param k number of random selects
     @param rand random number generator
     @param tolerance tolerance for equality of doubles
     @return indexes from randomly selected cdf distribution values
     */
    public static int[] chooseKFromBinarySearch(double[] cdf, int k, Random rand,
        double tolerance) {
        
        int[] selected = new int[k];
        
        for (int i = 0; i < k; ++i) {
            selected[i] = binarySearchForQuantile(cdf, rand.nextDouble());
        }
        
        return selected;
    }
 
    /**
     * Find the least index in the CDF where the probability is >= srch.
     * This is the quantile for the inverse mapping of the probability srcn.
     * 
     * The runtime is O(lg_2(cdf.length)).
     * 
     @param cdf cumulative distribution function.  NOTE that the values should
     * be normalized such that the last item is 1.
     @param srch the probability to search for in the CDF.
     @return index in the cdf distribution whose array value is nearest to
     * srch.
    */
    public static int binarySearchForQuantile(double[] cdf, final double srch) {

        /*
        if srch is equal to a point in the cdf,
        then we take the floor of that value to find
        the least index with that same value
        else if srch is not == to a point in the CDF,
        we pick the next highest value in the CDF.

        For a single objective instead of 1 for EQ and 1 for LT,
        we can use srch - very small number and perform a successor
        search for that.
         */
        
        if (cdf == null || cdf.length == 0) {
            throw new IllegalArgumentException("cdf cannot be null or length 0");
        }

        int n = cdf.length;

        // assuming machine precision is ~ 1E-11, use larger val for eps
        double srch2 = srch - 5E-10;
           
        int lowIdx = 0;
        int highIdx = n - 1;
        int midIdx;

        while (lowIdx <= highIdx) {

            midIdx = lowIdx + ((highIdx - lowIdx)/2);

            if (cdf[midIdx] <= srch2) {
                lowIdx = midIdx + 1;
            } else {
                highIdx = midIdx - 1;
            }

        }

        if (lowIdx == n) {
            lowIdx = n - 1;
        }
        
        return lowIdx;
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
     * NOTE: this method is considered for cases such as CURDecomposition
     * performed on a matrix with a very large number of parameters.
     * 
     @param cdf cumulative distribution function
     @param k number of random selects
     @param rand random number generator
     @return indexes from randomly selected cdf distribution values
     */
    /*public static int[] chooseKFromIntegerTransformAndTrie(double[] cdf, int k, Random rand) {       
        throw new UnsupportedOperationException("Not implemented yet."); 
    }*/

}
