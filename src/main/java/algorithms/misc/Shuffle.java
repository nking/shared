package algorithms.misc;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 * included is the Fisher-Yates algorithm which can be used to create unbiased random permutations.

   As the number of permutations of a sequence of integers of length n is n!, there are limits to
   calculating all permutations.  For n less than 20, could implement Heap's recursive algorithm if needed or
   Ives' iterative algorithm.

 * @author nichole
 */
public class Shuffle {
    
    /**
     * randomly shuffle the cards.  the distribution of permutations is gaussian (normal).
     * runtime complexity is O(n) where n = a.length.
     * @param a
     * @throws NoSuchAlgorithmException 
     */
    public static void fisherYates(double[] a) throws NoSuchAlgorithmException {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        fisherYates(a, rand);
    }
    
    /**
     * randomly shuffle the cards.  the distribution of permutations is gaussian (normal).
     * runtime complexity is O(n) where n = a.length.
     * @param a
     * @throws NoSuchAlgorithmException 
     */
    public static void fisherYates(int[] a) throws NoSuchAlgorithmException {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        fisherYates(a, rand);
    }
    
    /**
     * randomly shuffle the cards.  the distribution of permutations is gaussian (normal).
     * runtime complexity is O(n) where n = a.length.
     * @param a
     * @param rand 
     */
    public static void fisherYates(double[] a, SecureRandom rand) {
        
        int n = a.length;
        int j;
        double swap;
        for (int i = (n-1); i > 0; i--) {
            // 0 <= j <= i
            j = rand.nextInt(i + 1);
            swap = a[i];
            a[i] = a[j];
            a[j] = swap;
        }
    }
    
    /**
     * randomly shuffle the cards.  the distribution of permutations is gaussian (normal).
     * runtime complexity is O(n) where n = a.length.
     * @param a
     * @param rand 
     */
    public static void fisherYates(int[] a, SecureRandom rand) {
        
        int n = a.length;
        int j;
        int swap;
        for (int i = (n-1); i > 0; i--) {
            // 0 <= j <= i
            j = rand.nextInt(i + 1);// upper bound, exclusive, so need i+1
            swap = a[i];
            a[i] = a[j];
            a[j] = swap;
        }
    }
    
}
