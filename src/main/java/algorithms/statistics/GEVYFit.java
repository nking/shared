package algorithms.statistics;

import algorithms.misc.MiscMath0;

public class GEVYFit implements IYFit {

    protected float[] yfit;
    protected float[] x;
    protected float xScale = 1;
    protected float yScale = 1;
    protected int xPeakIndex = -1;
    protected float k;
    protected float sigma;
    protected float mu;
    protected float chiSqSum = Float.MAX_VALUE;
    protected float chiSqStatistic = Float.MAX_VALUE;
    float kSolutionResolution;
    float sigmaSolutionResolution;
    float muSolutionResolution;
    protected float yDataErrSq;

    protected String[] parameterNames = new String[]{
        "k", "sigma", "mu"
    };

    public long approximateMemoryUsed() {

        String arch = System.getProperty("sun.arch.data.model");

        boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;

        int nbits = (is32Bit) ? 32 : 64;

        int overheadBytes = 16;

        int intBytes = (is32Bit) ? 4 : 8;
        int arrayBytes = 32/8;

        long sumBytes = 0;

        if (yfit != null) {
            sumBytes += (arrayBytes + (yfit.length*arrayBytes));
        }

        if (x != null) {
            sumBytes += (arrayBytes + (x.length*arrayBytes));
        }

        // 17 variables on the stack, each of size stack word size
        sumBytes += (17 * intBytes);

        // String size on the heap = reference size + content size?
        // parameterNames
        sumBytes += (arrayBytes + (3*(nbits/8) + (intBytes*1 + intBytes*5 + intBytes*2)));

        sumBytes += overheadBytes;

        // amount of padding needed to make it a round 8 bytes
        long padding = (sumBytes % 8);

        sumBytes += padding;

        return sumBytes;
    }

    public String[] getParameterNames() {
        return parameterNames;
    }

    public float[] getParameters() {
        return new float[]{k, sigma, mu};
    }

    public String toString() {

        StringBuffer sb = new StringBuffer();
        sb.append(yfit.length).append(" points, k=").append(k).append(" sigma=").append(sigma)
            .append(" mu=").append(mu).append(" chiSqSum=").append(chiSqSum)
            .append(" chiSqStatistic=").append(chiSqStatistic);

        return sb.toString();
    }

    public float getXPeak() {
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
    
    public float[] getYFit() {
        return yfit;
    }

    public float getK() {
        return k;
    }

    public float getSigma() {
        return sigma;
    }

    public float getMu() {
        return mu;
    }

    public float getKResolution() {
        return kSolutionResolution;
    }

    public float getSigmaResolution() {
        return sigmaSolutionResolution;
    }
    
    public float getMuSolutionResolution() {
        return muSolutionResolution;
    }

    public float getChiSqSum() {
        return chiSqSum;
    }

    public float getYDataErrSq() {
        return yDataErrSq;
    }

    public float getChiSqStatistic() {
        // minus one if mean was computed from the data
        if (yfit == null) {
            return chiSqStatistic;
        }
        return chiSqSum / (yfit.length - 3 - 1);
    }

    protected float calculateArea(int index, boolean isStepFunction) {
        return calculateArea(x, yfit, index, isStepFunction, xScale, yScale);
    }

    protected static float calculateArea(float[] x, float[] y, int xyIndex,
                                         boolean isStepFunction, float xScaleFactor, float yScaleFactor) {

        if (x.length == 0) {
            return 0.0f;
        } else if (x.length == 1) {
            float w = xScaleFactor * x[0];
            float h = yScaleFactor * y[0];
            return w * h;
        }
        float x0, xmid, x1, y0, ymid, y1;
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
            float xDelta = ((x[xyIndex + 1] - x[xyIndex])/2.0f);
            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;
            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                float yDelta = (y[1] - y[0])/2.0f;
                y0 = y[0] - yDelta;
                y1 = y[0] + yDelta;
            }
        } else if (xyIndex == (y.length - 1)) {
            float xDelta = ((x[xyIndex] - x[xyIndex - 1])/2.0f);
            x0 = x[xyIndex] - xDelta;
            x1 = x[xyIndex] + xDelta;
            if (isStepFunction) {
                y0 = ymid;
                y1 = ymid;
            } else {
                float yDelta = (y[xyIndex - 1] - y[xyIndex])/2.0f;
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
        float areaBase = xScaleFactor * (x1 - x0) * yScaleFactor * (ymid);
        if (isStepFunction) {
            return areaBase;
        }
        float areaTop0 = 0.5f * xScaleFactor * (xmid - x0) * yScaleFactor * (y0 - ymid);
        float areaTop1 = 0.5f * xScaleFactor * (x1 - xmid) * yScaleFactor * (y1 - ymid);
        float area = areaTop0 + areaTop1 + areaBase;
        return area;
    }

    public float getX(int index) {
        return xScale*x[index];
    }

    /**
     * @param yfit the yfit to set
     */
    public void setYFit(float[] yfit) {
        this.yfit = yfit;
    }

    public void setYScale(float scale) {
        this.yScale = scale;
    }

    /**
     * @return the x array of the fit
     */
    public float[] getX() {
        return x;
    }

    /**
     * @param x array of the fit
     */
    public void setX(float[] x) {
        this.x = x;
    }

    public void setXScale(float scale) {
        this.xScale = scale;
    }

    /**
     * @param k the k to set
     */
    public void setK(float k) {
        this.k = k;
    }

    /**
     * @param sigma the sigma to set
     */
    public void setSigma(float sigma) {
        this.sigma = sigma;
    }

    /**
     * @param mu the mu to set
     */
    public void setMu(float mu) {
        this.mu = mu;
    }

    /**
     * @param chiSq the chiSqSum to set
     */
    public void setChiSqSum(float chiSq) {
        this.chiSqSum = chiSq;
    }

    /**
     * @param chiSqStat the chiSqStatistic to set
     */
    public void setChiSqStatistic(float chiSqStat) {
        this.chiSqStatistic = chiSqStat;
    }

    /**
     * @param yErrSq the yDataErrSq to set
     */
    public void setYDataErrSq(float yErrSq) {
        this.yDataErrSq = yErrSq;
    }

    public float[] getOriginalScaleX() {
        if (x == null) {
            return null;
        }
        float[] xsc = new float[x.length];
        for (int i = 0; i < xsc.length; i++) {
            xsc[i] = x[i] * xScale;
        }
        return xsc;
    }

    public float[] getOriginalScaleYFit() {
        if (yfit == null) {
            return null;
        }
        float[] ysc = new float[yfit.length];
        for (int i = 0; i < ysc.length; i++) {
            ysc[i] = yfit[i] * yScale;
        }
        return ysc;
    }

    public float getXScale() {
        return xScale;
    }

    public float getYScale() {
        return yScale;
    }

}
