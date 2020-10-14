package algorithms.misc;

import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;

/**
 *
 * @author nichole
 */
public class Shuffle {
    
    public static void fisherYates(double[] a) throws NoSuchAlgorithmException {
        
        SecureRandom rand = SecureRandom.getInstanceStrong();
        long seed = System.nanoTime();
        //System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        fisherYates(a, rand);
    }
    
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
    
}
