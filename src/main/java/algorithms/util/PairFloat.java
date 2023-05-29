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
        return FNVHash.hash32a(new float[]{this.x, this.y});
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder("(");
        sb.append(Float.toString(x)).append(",").append(Float.toString(y)).append(")");
        
        return sb.toString();
    }
}
