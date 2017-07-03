package algorithms.misc;

import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nichole
 */
public class Misc0 {
     
    /**
     * get an instance of SecureRandom, trying first
     * the algorithm SHA1PRNG, else the
     * default constructor.
     * @return 
     */
    public static SecureRandom getSecureRandom() {
        
        SecureRandom sr = null;
        
        try {
            sr = SecureRandom.getInstance("SHA1PRNG");
        } catch (NoSuchAlgorithmException ex) {
            Logger.getLogger(Misc0.class.getName()).log(Level.SEVERE, null, ex);
            sr = new SecureRandom();
        }
        
        return sr;
    }
    
    public static double[][] convertToBinary(TIntSet pixIdxs, int width, int height) {
        
        double[][] out = new double[width][height];
        for (int i = 0; i < width; ++i) {
            out[i] = new double[height];
        }
        
        PixelHelper ph = new PixelHelper();
        int[] xy = new int[2];
        
        TIntIterator iter = pixIdxs.iterator();
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            ph.toPixelCoords(pixIdx, width, xy);
        
            out[xy[0]][xy[1]] = 1;
        }
        
        return out;
    }
}
