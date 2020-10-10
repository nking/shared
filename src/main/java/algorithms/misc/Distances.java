package algorithms.misc;

/**
 *
 * @author nichole
 */
public class Distances {
    
    /**
     * calculate the square of the euclidean distance of p1 from p2, that
     * is sum of (p1[i]-p2[i])^2.
     * @param p1 a point with p1.length dimensions
     * @param p2 a point with p2.length dimensions
     * @return the euclidean difference between p1 and p2.
     */
    public static double calcEuclideanSquared(double[] p1, double[] p2) {
        
        if (p1.length != p2.length) {
            throw new IllegalArgumentException("p1.lenght must equal p2.length");
        }
        
        double d;
        double s = 0;
        for (int i = 0; i < p1.length; ++i) {
            d = p1[i] - p2[i];
            s += (d*d);
        }
        
        return s;
    }
}
