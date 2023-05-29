package algorithms.util;

/**
 *
    first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class PairInt {
    
    private int x = Integer.MIN_VALUE;
    private int y = Integer.MIN_VALUE;
    /**
     *
     */
    public PairInt() {
    }

    /**
     *
     @param xPoint
     @param yPoint
     */
    public PairInt(int xPoint, int yPoint) {
        x = xPoint;
        y = yPoint;
    }
    
    /**
     * constructor which rounds to the nearest integer
     @param xPoint
     @param yPoint 
     */
    public PairInt(double xPoint, double yPoint) {
        x = (int)Math.round(xPoint);
        y = (int)Math.round(yPoint);
    }
    
    /**
     *
     @param xPoint
     */
    public void setX(int xPoint) {
        x = xPoint;
    }

    /**
     *
     @param yPoint
     */
    public void setY(int yPoint) {
        y = yPoint;
    }

    /**
     *
     @return
     */
    public int getX() {
        return x;
    }

    /**
     *
     @return
     */
    public int getY() {
        return y;
    }
    
    /**
     *
     @return
     */
    public PairInt copy() {
        PairInt c = new PairInt(x, y);
        return c;
    }

    @Override
    public boolean equals(Object obj) {
        
        if (!(obj instanceof PairInt)) {
            return false;    
        }
        
        PairInt other = (PairInt)obj;

        return (x == other.getX()) && (y == other.getY());
    }

    @Override
    public int hashCode() {
        return FNVHash.hash32a(new int[]{this.x, this.y});
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Integer.toString(x)).append(",").append(Integer.toString(y)).append(")");
        
        return sb.toString();
    }
    
}
