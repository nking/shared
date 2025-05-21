package algorithms.signalProcessing;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.NumberTheory;

public class CurveResampler {

    /**
     * resample the curve xy to the size n2.  An upsample is
     * performed to the Least common multiple of xy[0].length and n2,
     * and then a downsample is performed to length n2.
     *
     * @param xy an array of size [2 X nPoints].
     * @param n2 the final sampled number of points
     * @return
     */
    public static int[][] upDownSample(int[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        int n1 = xy[0].length;
        long lcm = NumberTheory.leastCommonMultiple(n1, n2);

        double[][] xyD = MatrixUtil.convertIntToDouble(xy);

        double[][] up = linearUpsample(xyD, (int)lcm);
        double[][] down = linearUpsample(up, n2);

        return MatrixUtil.convertDoubleToInt(down);
    }

    /**
     * resample the curve xy to the size n2.  An upsample is
     * performed to the Least common multiple of xy[0].length and n2,
     * and then a downsample is performed to length n2.
     *
     * @param xy an array of size [2 X nPoints].
     * @param n2 the final sampled number of points
     * @return
     */
    public static double[][] upDownSample(double[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        int n1 = xy[0].length;
        long lcm = NumberTheory.leastCommonMultiple(n1, n2);

        double[][] up = linearUpsample(xy, (int)lcm);
        // consider adding a smoothing filter after upsample
        double[][] down = linearDownsample(up, n2);
        return down;
    }

    /**
     * linear down-sample of curve xy from nPoints = xy[0].length to n2 points where n2 < nPoints.
     * @param xy array of size [2 x nPoints]
     * @param n2
     * @return
     */
    public static double[][] linearDownsample(double[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        int n1 = xy[0].length;
        if (n1 < n2) {
            throw new IllegalArgumentException("n2 must be less than xy[0].length.  Use upSample instead");
        }
        double frac = (double) n1 / (double) n2;
        /**
         *  0  1  2  3  4  5 6 7 8 9 n1=10
         *     *  *  *  *  *  *      n2=6
         *    (n1-1)/(n2-1) = 1.8
         *    sample curve at intervals of 1.8 indexes in x
         *    1.8 floor =1, 1.8 ceil=2
         *       vals[1] + 0.8*(diff between vals at ceil and floor)
         */
        double[][] xy2 = new double[2][n2];
        xy2[0][0] = xy[0][0];
        xy2[1][0] = xy[1][0];
        int idx = 1;
        while (idx < n2) {
            double x = idx * frac;
            int flIdx = (int)Math.floor(x);
            int clIdx = (int)Math.ceil(x);
            if (clIdx >= n1) {
                xy2[0][idx] = xy[0][n1-1];
                xy2[1][idx] = xy[1][n1-1];
                ++idx;
            } else if (flIdx == clIdx) {
                xy2[0][idx] = xy[0][flIdx];
                xy2[1][idx] = xy[1][flIdx];
                ++idx;
            } else {
                double xdiff = x - flIdx;
                xy2[0][idx] = xy[0][flIdx] + (xdiff * (xy[0][clIdx] - xy[0][flIdx]));
                xy2[1][idx] = xy[1][flIdx] + (xdiff * (xy[1][clIdx] - xy[1][flIdx]));
                ++idx;
            }
        }
        return xy2;
    }

    /**
     * linear upsample of curve xy from nPoints = xy[0].length to n2 points where n2 > nPoints.
     * @param xy array of size [2 x nPoints]
     * @param n2
     * @return
     */
    public static double[][] linearUpsample(double[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        int n1 = xy[0].length;
        if (n1 > n2) {
            throw new IllegalArgumentException("n2 must be larger than xy[0].length.  Use downSample instead");
        }
        double frac = (double)n1/(double)n2;

        double[][] xy2 = new double[2][n2];
        xy2[0][0] = xy[0][0];
        xy2[1][0] = xy[1][0];
        int idx = 1;

        for (int i = 0; i < n1 - 1; ++i) {
            double ax = xy[0][i];
            double bx = xy[0][i+1];
            double ay = xy[1][i];
            double by = xy[1][i+1];
            int j = 1;
            while (idx < n2) {
                double x = ax + j * (bx - ax) * frac;
                double y = ay + j * (by - ay) * frac;
                xy2[0][idx] = x;
                xy2[1][idx] = y;
                ++j;
                ++idx;
                if (ax < bx && !(x < bx)) {
                    break;
                } else if (ax >= bx && !(x >= ax)) {
                    break;
                }
            }
            int tt = 2;
        }

       return xy2;
    }
}
