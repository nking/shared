package algorithms.statistics;

import algorithms.misc.MiscMath0;

public class GEVYFit implements IYFit {

    protected double[] yfit = null;
    protected double[] x = null;
    protected double xScale = 1;
    protected double yScale = 1;
    protected int xPeakIndex = -1;
    protected double k;
    protected double sigma;
    protected double mu;
    protected double chiSqSum = Double.POSITIVE_INFINITY;
    protected double chiSqStatistic = Double.POSITIVE_INFINITY;
    double kSolutionResolution;
    double sigmaSolutionResolution;
    double muSolutionResolution;
    protected double yDataErrSq;

    protected String[] parameterNames = new String[]{
        "mu", "sigma", "k"
    };

    public String[] getParameterNames() {
        return parameterNames;
    }

    public double[] getParameters() {
        return new double[]{mu, sigma, k};
    }

    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append(yfit.length).append(" points, mu=").append(mu).append(" sigma=").append(sigma)
            .append(" k=").append(k).append(" chiSqSum=").append(chiSqSum)
            .append(" chiSqStatistic=").append(chiSqStatistic);

        return sb.toString();
    }

    public double getXPeak() {
        if (xPeakIndex == -1) {
            xPeakIndex = MiscMath0.findYMaxIndex(yfit);
        }
        return xScale*x[xPeakIndex];
    }
    
    public int getXPeakIndex() {
        if (xPeakIndex == -1) {
            xPeakIndex = MiscMath0.findYMaxIndex(yfit);
        }
        return xPeakIndex;
    }
    
    public double[] getYFit() {
        return yfit;
    }

    public double getK() {
        return k;
    }

    public double getSigma() {
        return sigma;
    }

    public double getMu() {
        return mu;
    }

    public double getKResolution() {
        return kSolutionResolution;
    }

    public double getSigmaResolution() {
        return sigmaSolutionResolution;
    }
    
    public double getMuSolutionResolution() {
        return muSolutionResolution;
    }

    public double getChiSqSum() {
        return chiSqSum;
    }

    public double getYDataErrSq() {
        return yDataErrSq;
    }

    public double getChiSqStatistic() {
        // minus one if mean was computed from the data
        if (yfit == null) {
            return chiSqStatistic;
        }
        return chiSqSum / (yfit.length - 3 - 1);
    }

    protected double calculateArea(int index, boolean isStepFunction) {
        return calculateArea(x, yfit, index, isStepFunction, xScale, yScale);
    }

    protected static double calculateArea(double[] x, double[] y, int xyIndex,
         boolean isStepFunction, double xScaleFactor, double yScaleFactor) {

        if (x.length == 0) {
            return 0.0f;
        } else if (x.length == 1) {
            double w = xScaleFactor * x[0];
            double h = yScaleFactor * y[0];
            return w * h;
        }
        double x0, xmid, x1, y0, ymid, y1;
        /*                              *
         *                          .      .
         *   *         *         *  ........  *
         *     .     .              .      .
         *     ...*...              .      .
         *     .     .
         *     .     .
         *
         */

        xmid = x[xyIndex];
        ymid = y[xyIndex];

        if (xyIndex == 0) {
            double xDelta = ((x[xyIndex + 1] - x[xyIndex])/2.0f);
            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;
            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                double yDelta = (y[1] - y[0])/2.0f;
                y0 = y[0] - yDelta;
                y1 = y[0] + yDelta;
            }
        } else if (xyIndex == (y.length - 1)) {
            double xDelta = ((x[xyIndex] - x[xyIndex - 1])/2.0f);
            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;
            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                double yDelta = (y[xyIndex - 1] - y[xyIndex])/2.0f;
                y0 = ymid + yDelta;
                y1 = ymid - yDelta;
            }
        } else {
            x0 = ((x[xyIndex] + x[xyIndex - 1]) / 2.0f);
            x1 = ((x[xyIndex + 1] + x[xyIndex]) / 2.0f);
            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                y0 = (y[xyIndex - 1] + y[xyIndex]) / 2.0f;
                y1 = (y[xyIndex] + y[xyIndex + 1]) / 2.0f;
            }
        }
        double areaBase = xScaleFactor * (x1 - x0) * yScaleFactor * (ymid);
        if (isStepFunction) {
            return areaBase;
        }
        double areaTop0 = 0.5f * xScaleFactor * (xmid - x0) * yScaleFactor * (y0 - ymid);
        double areaTop1 = 0.5f * xScaleFactor * (x1 - xmid) * yScaleFactor * (y1 - ymid);
        double area = areaTop0 + areaTop1 + areaBase;
        return area;
    }

    public double getX(int index) {
        return xScale*x[index];
    }

    /**
     * @param yfit the yfit to set
     */
    public void setYFit(double[] yfit) {
        this.yfit = yfit;
    }

    public void setYScale(double scale) {
        this.yScale = scale;
    }

    /**
     * @return the x array of the fit
     */
    public double[] getX() {
        return x;
    }

    /**
     * @param x array of the fit
     */
    public void setX(double[] x) {
        this.x = x;
    }

    public void setXScale(double scale) {
        this.xScale = scale;
    }

    /**
     * @param k the k to set
     */
    public void setK(double k) {
        this.k = k;
    }

    /**
     * @param sigma the sigma to set
     */
    public void setSigma(double sigma) {
        this.sigma = sigma;
    }

    /**
     * @param mu the mu to set
     */
    public void setMu(double mu) {
        this.mu = mu;
    }

    /**
     * @param chiSq the chiSqSum to set
     */
    public void setChiSqSum(double chiSq) {
        this.chiSqSum = chiSq;
    }

    /**
     * @param chiSqStat the chiSqStatistic to set
     */
    public void setChiSqStatistic(double chiSqStat) {
        this.chiSqStatistic = chiSqStat;
    }

    /**
     * @param yErrSq the yDataErrSq to set
     */
    public void setYDataErrSq(double yErrSq) {
        this.yDataErrSq = yErrSq;
    }

    public double[] getOriginalScaleX() {
        if (x == null) {
            return null;
        }
        double[] xsc = new double[x.length];
        for (int i = 0; i < xsc.length; i++) {
            xsc[i] = x[i] * xScale;
        }
        return xsc;
    }

    public double[] getOriginalScaleYFit() {
        if (yfit == null) {
            return null;
        }
        double[] ysc = new double[yfit.length];
        for (int i = 0; i < ysc.length; i++) {
            ysc[i] = yfit[i] * yScale;
        }
        return ysc;
    }

    public double getXScale() {
        return xScale;
    }

    public double getYScale() {
        return yScale;
    }

}
