package algorithms.util;

import java.util.Arrays;

/**
 * a class to enclose a primitive int array and provide an
 * equals and hashcode which will have the same identity for two
 * instances that have the same ordered content.
 * @author nichole
 */
public class OneDIntArray {

    /**
     *
     */
    public int[] a;

    /**
     *
     @param t
     */
    public OneDIntArray(int[] t) {
        a = t;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof OneDIntArray)) {
            return false;    
        }
        OneDIntArray other = (OneDIntArray)obj;
        
        return Arrays.equals(a, other.a);
    }

    @Override
    public int hashCode() {
        return FNVHash.hash(a);
    }

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
}
