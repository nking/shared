package algorithms.util;

import java.util.Arrays;

/**
 class to hold a triple set of integers for use in Collections.
 2 different instances of TriInt(1,2,3) are equal to one another.
 * @author nichole
 */
public class PolyInt {

    private final int[] a;

    public PolyInt(int[] points) {
        this.a = Arrays.copyOf(points, points.length);
    }

    public int[] get() {
        return a;
    }

    /**
     *
     @return
     */
    public PolyInt copy() {
        return new PolyInt(Arrays.copyOf(a, a.length));
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PolyInt)) {
            return false;    
        }
        
        PolyInt other = (PolyInt)obj;

        if (a.length != other.get().length) {
            return false;
        }
        for (int i = 0; i < a.length; ++i) {
            if (a[i] != other.a[i]) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int hashCode() {
        return FNVHash.hash32a(a);
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Arrays.toString(a)).append(")");
        
        return sb.toString();
    }
    
}
