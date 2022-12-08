package algorithms.util;

/**
 *
 * @author nichole
 */
public class PairFloat {
    
    private float x;
    private float y;
    
    /**
     *
     */
    public PairFloat() {
    }

    /**
     *
     @param xPoint
     @param yPoint
     */
    public PairFloat(float xPoint, float yPoint) {
        x = xPoint;
        y = yPoint;
    }

    /**
     *
     @param xPoint
     */
    public void setX(float xPoint) {
        x = xPoint;
    }

    /**
     *
     @param yPoint
     */
    public void setY(float yPoint) {
        y = yPoint;
    }

    /**
     *
     @return
     */
    public float getX() {
        return x;
    }

    /**
     *
     @return
     */
    public float getY() {
        return y;
    }
    
    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairFloat)) {
            return false;    
        }
        
        PairFloat other = (PairFloat)obj;
        
        return (x == other.getX()) && (y == other.getY());
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(this.x, this.y);

        return hash;
    }

    /**
     *
     */
    protected final static int fnv321aInit = 0x811c9dc5;

    /**
     *
     */
    protected final static int fnv32Prime = 0x01000193;

    /**
     *
     @param i0
     @param i1
     @return
     */
    protected int fnvHashCode(float i0, float i1) {

        /*
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
         */

        int hash = 0;

        int sum = fnv321aInit;

        // xor the bottom with the current octet.
        sum ^= Float.floatToIntBits(i0);

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        sum ^= Float.floatToIntBits(i1);
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Float.toString(x)).append(",").append(Float.toString(y)).append(")");
        
        return sb.toString();
    }
}
