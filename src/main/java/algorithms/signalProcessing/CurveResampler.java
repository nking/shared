package algorithms.signalProcessing;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.NumberTheory;
import algorithms.util.FormatArray;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.*;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;


public class CurveResampler {

    public static int[][] resample(int[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        double[][] xyD = MatrixUtil.convertIntToDouble(xy);
        double[][] r = resample(xyD, n2);
        return MatrixUtil.convertDoubleToInt(r);
    }

    public static float[][] resample(float[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }
        double[][] xyD = MatrixUtil.convertToDouble(xy);
        double[][] r = resample(xyD, n2);
        return MatrixUtil.convertToFloat(r);
    }

    public static double[][] resample(double[][] xy, int n2) {
        if (xy.length < 2) {
            throw new IllegalArgumentException("xy length must be 2");
        }

        int n1 = xy[0].length;
        double[] xLine = new double[n1];
        for (int i = 0; i < n1; ++i) {
            xLine[i] = i;
        }

        double[] x2Line = new double[n2];
        double frac = (double) (n1 - 1.) / (double) (n2 - 1.);
        for (int i = 0; i < n2; ++i) {
            x2Line[i] = i * frac;
        }

        SplineInterpolator interpolator = new SplineInterpolator();
        PolynomialSplineFunction xInterp = interpolator.interpolate(xLine, xy[0]);
        PolynomialSplineFunction yInterp = interpolator.interpolate(xLine, xy[1]);

        double[][] xy2 = new double[2][n2];
        for (int i = 0; i < n2; i++) {
            xy2[0][i] = xInterp.value(x2Line[i]);
            xy2[1][i] = yInterp.value(x2Line[i]);
        }
        return xy2;
    }

}
