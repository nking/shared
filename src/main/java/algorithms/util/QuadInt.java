package algorithms.util;

/**
 *
 * @author nichole
 */
public class QuadInt {
    
    private int a = Integer.MIN_VALUE;
    private int b = Integer.MIN_VALUE;
    private int c = Integer.MIN_VALUE;
    private int d = Integer.MIN_VALUE;
    
    /**
     *
     */
    public QuadInt() {
    }

    /**
     *
     @param aPoint
     @param bPoint
     @param cPoint
     @param dPoint
     */
    public QuadInt(int aPoint, int bPoint, int cPoint, int dPoint) {
        a = aPoint;
        b = bPoint;
        c = cPoint;
        d = dPoint;
    }

    /**
     *
     @param ab
     @param cd
     */
    public QuadInt(PairInt ab, PairInt cd) {
        a = ab.getX();
        b = ab.getY();
        c = cd.getX();
        d = cd.getY();
    }

    /**
     *
     @param aPoint
     */
    public void setA(int aPoint) {
        a = aPoint;
    }

    /**
     *
     @param bPoint
     */
    public void setB(int bPoint) {
        b = bPoint;
    }

    /**
     *
     @param cPoint
     */
    public void setC(int cPoint) {
        c = cPoint;
    }

    /**
     *
     @param dPoint
     */
    public void setD(int dPoint) {
        d = dPoint;
    }

    /**
     *
     @return
     */
    public int getA() {
        return a;
    }

    /**
     *
     @return
     */
    public int getB() {
        return b;
    }

    /**
     *
     @return
     */
    public int getC() {
        return c;
    }

    /**
     *
     @return
     */
    public int getD() {
        return d;
    }
    
    /**
     *
     @return
     */
    public QuadInt copy() {
        QuadInt cp = new QuadInt(a, b, c, d);
        return cp;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof QuadInt)) {
            return false;    
        }
        
        QuadInt other = (QuadInt)obj;
        
        return (a == other.getA()) && (b == other.getB()) 
            && (c == other.getC()) && (d == other.getD());
    }

    @Override
    public int hashCode() {
        
        int hash = fnvHashCode(a, b, c, d);

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
     @param i2
     @param i3
     @return
     */
    protected int fnvHashCode(int i0, int i1, int i2, int i3) {

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
        sum ^= i0;

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        sum ^= i1;
        
        sum *= fnv32Prime;
        
        sum ^= i2;
        
        sum *= fnv32Prime;
        
        sum ^= i3;
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Integer.toString(a)).append(",")
            .append(Integer.toString(b)).append("), (")
            .append(Integer.toString(c)).append(",")
            .append(Integer.toString(d)).append(")");
        
        return sb.toString();
    }
    
}
