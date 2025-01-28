package algorithms.misc;

/**
 *
 * @author nichole
 */
public class Distances {
    
    /**
     * calculate the square of the euclidean distance of p1 from p2, that
     * is sum of (p1[i]-p2[i])^2.
     @param p1 a point with p1.length dimensions
     @param p2 a point with p2.length dimensions
     @return the euclidean difference between p1 and p2.
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

    public static double manhattan(double x1, double y1, double x2, double y2) {
        //return Math.abs(x1-x2) + Math.abs(y1-y2);
        return Math.max(x1-x2, -x1+x2) + Math.max(y1-y2, -y1+y2);
    }

    public static double manhattan2(double x1, double y1, double x2, double y2) {
        /*from Competitive Programmer's Handbook by Antti Laaksonen, Chap 29.4
        rotate the coordinates by 45 degrees.

        from https://cp-algorithms.com/geometry/manhattan-distance.html

        rotation by 45 degrees about origin (0,0) and scale by sqrt(2) for the diagonal.
        that transformation reduces to
           x' = x + y
           y' = x - y

        Note that the Chebyshev distance is between pairs of points transformed by rot=45, scale=sqrt(2).
        */
        double x1p = x1 + y1;
        double y1p = x1 - y1;
        double x2p = x2 + y2;
        double y2p = x2 - y2;
        return Math.max(Math.abs(x1p - x2p), Math.abs(y1p - y2p));
    }

    /**
     * given an array of points of dimension d, find the maximum distance between a pair of points.
     * The runtime complexity is O(n * (2^d) * d)).
     * e.g. for d=2, the runtime complexity is O(n*8).
     <pre>
     reference:
     https://cp-algorithms.com/geometry/manhattan-distance.html
     </pre>
     * @param points array of points, that is, each row is a point of d dimenstions.
     * @return indexes of the furthest pair, or indexes of 1 furthest pair if more than one.
     */
    public static int[] maxManhattan(int[][] points) {
        /*
        for points p and q:
           considering x points alone: abs(xp - xq) = max(xp-xq, -xp+xq)
           for all p,q in points: max(xp-xq) = (max(xp) for p in points)
               + (max(-xq) for q in points)

           for 2 dimensions:
              for all p,q in points: max((xp-xq) + (yp-yq))
                = (max(xp + yp) for p in points)
               + (max(-xq - yq) for q in points)
         */
        int d = points[0].length;
        int n = points.length;

        long ans = 0;

        int[] out = new int[2];

        for (int msk = 0; msk < (1 << d); msk++) {
            long mx = Long.MIN_VALUE;
            long mn = Long.MAX_VALUE;
            System.out.println("");
            for (int i = 0; i < n; i++) {
                long cur = 0;
                for (int j = 0; j < d; j++) {
                    if ((msk & (1 << j)) != 0) {
                        cur += points[i][j];
                    } else {
                        cur -= points[i][j];
                    }
                }

                if (mx < cur) {
                    //System.out.printf("mx %d, i=%d, msk=%d\n", cur, i, msk);
                    out[1] = i;
                    mx = cur;
                }
                if (mn > cur) {
                    //System.out.printf("mn %d, i=%d, msk=%d\n", cur, i, msk);
                    out[0] = i;
                    mn = cur;
                }
                //mx = Math.max(mx, cur);
                //mn = Math.min(mn, cur);
            }
            ans = Math.max(ans, mx - mn);
        }
        return out;
    }
}
