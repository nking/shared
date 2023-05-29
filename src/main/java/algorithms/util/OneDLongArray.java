package algorithms.util;

import java.util.Arrays;

/**
 * a class to enclose a primitive int array and provide an
 * equals and hashcode which will have the same identity for two
 * instances that have the same ordered content.
 * @author nichole
 */
public class OneDLongArray {

    /**
     *
     */
    public long[] a;

    /**
     *
     @param t
     */
    public OneDLongArray(long[] t) {
        a = t;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof OneDLongArray)) {
            return false;    
        }
        OneDLongArray other = (OneDLongArray)obj;
        
        return Arrays.equals(a, other.a);
    }

    @Override
    public int hashCode() {
        return FNVHash.hash64a(a);
    }

    @Override
    public String toString() {
        return Arrays.toString(a);
    }
}
