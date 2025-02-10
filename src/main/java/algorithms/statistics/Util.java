package algorithms.statistics;

public class Util {

    /**
     * given data points x and y representing a two dimensional
     * dataset, calculate the moments array.
     * The moments array is the sum of the followng combined terms:
     <pre>
     x    y    x^2    x^3    x^4    xy    yx^2  = moments
     0    1    2      3       4     5     6     = indexes in output
     </pre>
     * @param x array of 1st dimension of data
     * @param y array of 2nd dimension of data
     * @return moments
     <pre>
           x    y    x^2    x^3    x^4    xy    yx^2  = moments
           0    1    2      3       4     5     6     = indexes in output
     </pre>
     */
    public static double[] caldc2DMomentsX2Y4(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;
        // 0    1    2      3       4     5     6
        //x    y    x^2    x^3    x^4    xy    yx^2
        double[] moments = new double[7];
        for (int i = 0; i < n; ++i) {
            moments[0] += x[i];
            moments[1] += y[i];
            moments[2] += x[i] * x[i];
            moments[3] += moments[2] * x[i];
            moments[4] += moments[3] * x[i];
            moments[5] += x[i] * y[i];
            moments[6] += y[i] * x[i] * x[i];
        }
        return moments;
    }

    /**
     * given data points x and y representing a two dimensional
     * dataset, calculate the moments array.
     * The moments array is the sum of the followng combined terms:
     <pre>
     x    y    x^2    y^2    xy   = moments
     0    1    2      3       4   = indexes in output
     </pre>
     * @param x array of 1st dimension of data
     * @param y array of 2nd dimension of data
     * @return moments
    <pre>
    x    y    x^2    y^2    xy   = moments
    0    1    2      3       4   = indexes in output
    </pre>
     */
    public static double[] caldc2DMomentsX2Y2(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must be same length");
        }
        int n = x.length;
        //x    y    x^2    y^2    xy   = moments
        //0    1    2      3       4   = indexes in output
        double[] moments = new double[5];
        for (int i = 0; i < n; ++i) {
            moments[0] += x[i];
            moments[1] += y[i];
            moments[2] += x[i] * x[i];
            moments[3] += y[i] * y[i];
            moments[4] += x[i] * y[i];
        }
        return moments;
    }

}
