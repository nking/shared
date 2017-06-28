package algorithms.misc;

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
    
}
