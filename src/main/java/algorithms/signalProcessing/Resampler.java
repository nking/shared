package algorithms.signalProcessing;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.NumberTheory;
import algorithms.util.FormatArray;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class Resampler {

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

        double[][] up = resample(xyD, (int)lcm);
        double[][] down = resample(up, n2);

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

        double[][] up = resample(xy, (int)lcm);
        double[][] down = resample(up, n2);
        return down;
    }

    public static double[][] resample(double[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        PolynomialSplineFunction function = interpolator.interpolate(xy[0], xy[1]);

        // 2. Generate new x-values for the desired output length
        double minX = xy[0][0];
        double maxX = xy[0][xy[0].length - 1];
        double[] newX = new double[n2];
        for (int i = 0; i < n2; i++) {
            newX[i] = minX + (maxX - minX) * ((double) i / (n2 - 1));
        }

        // 3. Interpolate to get new y-values
        double[] newY = new double[n2];
        for (int i = 0; i < n2; i++) {
            newY[i] = function.value(newX[i]);
        }

        // 4. Combine new x and y into output array
        double[][] output = new double[2][n2];
        output[0] = newX;
        output[1] = newY;
        return output;
    }
}
