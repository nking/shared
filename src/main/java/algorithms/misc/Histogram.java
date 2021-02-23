package algorithms.misc;

import algorithms.sort.CountingSort;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TDoubleIntMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TIntIntHashMap;
import gnu.trove.set.TIntSet;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Logger;
import algorithms.sort.MiscSorter;

/**
 *  TODO:  improve this class...
 *
    first implemented in projects
     https://github.com/nking/two-point-correlation
     w/ Copyright (c) 2013-2015 Nichole King
     http://nking.github.io/two-point-correlation/
     using The MIT License (MIT)
     and
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class Histogram {

    protected static Logger log = Logger.getLogger(Histogram.class.getName());

    /**
     * create a histogram from the data that has little or no adjustment
     * for min and max.
     *
     * @param a
     * @param nBins
     * @param xHist
     * @param yHist
     */
    public static void createHistogram(float[] a, int nBins, float[] xHist, int[] yHist) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }

        float[] minMax = MiscMath0.calculateOuterRoundedMinAndMax(a);

        createHistogram(a, nBins, minMax[0], minMax[1], xHist, yHist);
    }

    protected static float calculateBinWidth(float minValue, float maxValue, int nBins) {

        float xInterval = (maxValue - minValue)/(float)nBins;

        // expand interval if necessary to make sure the last point is in the last bin
        if ((int) ((maxValue - minValue)/xInterval) != (nBins - 1)) {
            float t = (maxValue + minValue)/2.0f;
            int powDelta = MiscMath0.findPowerOf10(t);
            float pow10 = (float)Math.pow(10, powDelta);
            xInterval = (maxValue - minValue + pow10)/(float)nBins;
        }

        return xInterval;
    }

    public static void createHistogram(float[] a, int nBins,
        float aMin, float aMax, float[] xHist, int[] yHist, float binWidth) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (xHist == null || xHist.length != nBins) {
            throw new IllegalArgumentException("xHist has to be of size nBins and initialized");
        }
        if (yHist == null || yHist.length != nBins) {
            throw new IllegalArgumentException("yHist has to be of size nBins and initialized");
        }

        Arrays.fill(yHist, 0);

        for (int i = 0; i < nBins; i++) {
            xHist[i] = aMin + (float)i*binWidth + (binWidth/2.f);
        }

        for (int i = 0; i < a.length; i++) {
            int bin = (int) ((a[i] - aMin)/binWidth);
            if ((bin > -1) && (bin < nBins)) {
                yHist[bin]++;
            }
        }
    }
    
    public static void createHistogram(float[] a, int nBins,
        float aMin, float aMax, float[] xHist, int[] yHist) {

        if (a == null) {
            throw new IllegalArgumentException("a cannot be null");
        }
        if (xHist == null || xHist.length != nBins) {
            throw new IllegalArgumentException("xHist has to be of size nBins and initialized");
        }
        if (yHist == null || yHist.length != nBins) {
            throw new IllegalArgumentException("yHist has to be of size nBins and initialized");
        }

        Arrays.fill(yHist, 0);

        float xInterval = calculateBinWidth(aMin, aMax, nBins);

        createHistogram(a, nBins, aMin, aMax, xHist, yHist, xInterval);
    }
    
     public static HistogramHolder createSimpleHistogram(List<Integer> values) {

        if (values == null) {
            throw new IllegalArgumentException(
            "values and valueErrors cannot be null and must be the same length");
        }
        if (values.isEmpty()) {
            return null;
        }
        float[] v = new float[values.size()];
        for (int i = 0; i < values.size(); ++i) {
            v[i] = values.get(i).intValue();
        }
        float[] ve = Errors.populateYErrorsBySqrt(v);
        
        int nBins = (int)(2*Math.pow(v.length, 0.3333));
        
        if (v.length == 1) {
            nBins = 1;
        }
        
        return createSimpleHistogram(nBins, v, ve);
    }
   
    public static HistogramHolder createSimpleHistogram(float[] values, 
        float[] valueErrors) {

        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException(
            "values and valueErrors cannot be null and must be the same length");
        }
        
        int nBins = (int)(2*Math.pow(values.length, 0.3333));
        
        if (values.length == 1) {
            nBins = 1;
        }
        
        return createSimpleHistogram(nBins, values, valueErrors);
    }
     
    /**
     * WARNING:this is an incomplete method that is lossy too.
     * 
     * @param a
     * @param binSize
     * @param aMin
     * @param aMax
     * @return 
     */
    public static HistogramHolder createSimpleHistogram(double[] a, 
        double binSize, double aMin, double aMax) {
        if (a == null) {
            throw new IllegalArgumentException(
            "a cannot be null and must be the same length");
        }
        
        int nBins = (int)Math.ceil(((aMax - aMin))/binSize);
        if (nBins < 0) {
            nBins *= -1;
        }
        
        float[] afloat = new float[a.length];
        for (int i = 0; i < a.length; ++i) {
            afloat[i] = (float)a[i];
        }
        
        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(afloat, nBins, (float)aMin, (float)aMax, 
            xHist, yHist, (float)binSize);
        
        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHist);
        
        return histogram;
    }
      
    public static HistogramHolder createSimpleHistogram(int nBins, 
        float[] values, float[] valueErrors) {

        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }

        float[] minMax = MiscMath0.calculateOuterRoundedMinAndMax(values);
                
        float binWidth = calculateBinWidth(minMax[0], minMax[1], nBins);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, minMax[0], minMax[1], 
            xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }
    
    public static HistogramHolder createSimpleHistogram(float binWidth, 
        float[] values, float[] valueErrors) {

        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        
        float[] minMax = MiscMath0.calculateOuterRoundedMinAndMax(values);
                
        int nBins = (int)Math.ceil(((minMax[1] - minMax[0]))/binWidth);
        if (nBins < 0) {
            nBins *= -1;
        }

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, minMax[0], minMax[1], xHist, 
            yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }

    public static HistogramHolder createSimpleHistogram(float binWidth, 
        float[] values) {

        if (values == null) {
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        float[] minMax = MiscMath0.getMinMax(values);
        float xmax = MiscMath0.roundUpByLargestPower((float)minMax[1]);
        float xmin = MiscMath0.roundDownByLargestPower((float)minMax[0]);
        // xmax > 1 and xmin is between 0 and 1, round xmin down
        if ((xmax > 1) && (xmin > 0) && (xmin < 1.0)) {
            xmin = 0;
        }
        
        int nBins = (int)Math.ceil(((xmax - xmin))/binWidth);
        if (nBins < 0) {
            nBins *= -1;
        }

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, xmin, xmax, xHist, 
            yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        
        return histogram;
    }
    
    public static HistogramHolder createSimpleHistogram(float binWidth, 
        float[] values, float min, float max) {

        if (values == null) {
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        
        int nBins = (int)Math.ceil(((max - min))/binWidth);
        if (nBins < 0) {
            nBins *= -1;
        }

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, min, max, xHist, 
            yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        
        return histogram;
    }

    
    public static HistogramHolder createSimpleHistogram(float minX, float maxX,
        float binWidth, 
        float[] values, float[] valueErrors) {

        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        
        int nBins = (int)Math.ceil(((maxX - minX))/binWidth);
        if (nBins < 0) {
            nBins *= -1;
        }

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, minX, maxX, xHist, 
            yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }
    
    public static HistogramHolder createSimpleHistogram(int binWidth, 
        List<Integer> theValues) {

        if (theValues == null || theValues.isEmpty()) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        
        float[] values = new float[theValues.size()];
        for (int i = 0; i < theValues.size(); ++i) {
            int v = theValues.get(i).intValue();
            values[i] = v;
        }
        
        float[] valueErrors = Errors.populateYErrorsBySqrt(values);
        
        float[] minMax = MiscMath0.calculateOuterRoundedMinAndMax(values);
                
        int nBins = (int)Math.ceil(((minMax[1] - minMax[0]))/binWidth);
        if (nBins < 0) {
            nBins *= -1;
        }
        
        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, minMax[0], minMax[1], 
            xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, 
            yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }

    public static HistogramHolder createSimpleHistogram(
        final float xMin, final float xMax, int nBins, 
        float[] values, float[] valueErrors) {

        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
                
        float binWidth = calculateBinWidth(xMin, xMax, nBins);

        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
        
        Histogram.createHistogram(values, nBins, xMin, xMax,
            xHist, yHist, binWidth);

        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }

    /**
     * Populate the arrays xHistErrorsOutput and yHistErrorsOutput with errors
     * for the x and y bins calculated from valueErrors.
     *
     *
     * Trying to adjust errors to use the fact that a bin with a value of zero is
     * actually determined from all bins, so the error in a bin with a value of
     * zero is a positive number less than infinity and is due to all errors of
     * points that go into making the histogram.
     *
     * For each bin:
     *   (sigma)^2 = (ave sigma_from_all)^2 + (sigma from all points in bin, in function that reduces it by the number of points)
     *
     * The total of all bins should be equal to the sigma_from_all;
     *
     * N = 3 bins:
     *    [0] (sigma_bin0)^2 = (sigma_from_all/nBins)^2 + (sigma_bin0 * Fnc(bin0))^2
     *    [1] (sigma_bin1)^2 = (sigma_from_all/nBins)^2 + (sigma_bin1 * Fnc(bin1))^2
     *    [2] (sigma_bin2)^2 = (sigma_from_all/nBins)^2 + (sigma_bin2 * Fnc(bin2))^2
     * Where (sigma_bin0)^2 is the sum of all points that went into bin0 added in quadrature.
     *
     * The total of all 3 should equal the total from all points added in quadrature and no more.
     *
     *    (sigma_from_all_points)^2 = nBins*(sigma_from_all/nBins)^2 + (sigma_bin0 * Fnc(bin0))^2
     *                                + (sigma_bin1 * Fnc(bin1))^2 + (sigma_bin2 * Fnc(bin2))^2
     *
     *    (sigma_from_all_points)^2 * (1 - (1/nBins)) = (sigma_bin0 * Fnc(bin0))^2 + (sigma_bin1 * Fnc(bin1))^2 + (sigma_bin2 * Fnc(bin2))^2
     *
     * Not easy to see how to solve, so returning to first principles:
     *     Each bin should have a contribution from all points and then a contribution from its own.
     *     In order for the total to not exceed the sum in quadrature of all points, a bin's sigma squared should
     *     be the ave from all + it's own times (1-1/N) roughly, which is what the equation above suggests.
     *
     * @param xHist
     * @param yHist
     * @param values
     * @param valueErrors errors that are on the same scale as the values, that is, these are
     *   NOT percent errors
     * @param xHistErrorsOutput
     * @param yHistErrorsOutput
     */
    public static void calulateHistogramBinErrors(float[] xHist, int[] yHist,
        float[] values, float[] valueErrors, float[] xHistErrorsOutput, 
        float[] yHistErrorsOutput) {
        
        if ((xHist == null) || (xHist.length == 0)) {
            return;
        }

        float xInterval = (xHist.length > 1) ? xHist[1] - xHist[0] : 0;
        float xmin = xHist[0] - (xInterval/2.0f);

        float[] sumErrorPerBin = new float[xHist.length];

        float sumErrorAllPoints = 0;
        //float sumPercentErrorAllPoints = 0;

        for (int i = 0; i < valueErrors.length; i++) {

            int bin = (int) ((values[i] - xmin)/xInterval);

            // in units of y
            float a = valueErrors[i];
            a *= a;

            if ((bin > -1) && (bin < xHist.length)) {

                sumErrorPerBin[bin] += a;

                sumErrorAllPoints += a;

                //float b = valueErrors[i]/values[i];
                //sumPercentErrorAllPoints += (b*b);
            }
        }

        // the x value was determined from all points, so the error should be taken
        //   as the average error of all points
        float aveErrorOfAllPoints = (float)Math.sqrt(sumErrorAllPoints)/yHist.length;

        // the percent errors, divided over all bins, which is what was done to learn the binwdith
        //float avePercentErrorOfAllPoints = (float)Math.sqrt(sumPercentErrorAllPoints)/yHist.length;

        float sumAlt = 0;
        float sumOrigSquared = 0;
        for (int i = 0; i < yHist.length; i++) {

            float sumErrorBinSquared = sumErrorPerBin[i];

            // compare to aveErrorOfAllPoints... should be similar
            //float be = (float) Math.sqrt(sumErrorBinSquared);

            float c = aveErrorOfAllPoints;

            // estimating it as binWidth/2. any point has x value = bin center +- binWidth/2
            xHistErrorsOutput[i] = xInterval/2.0f; //c;

            //float a = (float)Math.sqrt(sumErrorBinSquared * (1.0f - (1.0f/yHist.length));
            float ai = sumErrorBinSquared;
            float af = sumErrorBinSquared/yHist.length;
            float a = (float)Math.sqrt(ai - af);

            float yBinError = (yHist[i] == 0) ? 0 : a/yHist[i];

            yHistErrorsOutput[i] = yBinError;

            sumAlt += yBinError;


            float yBinErrorOrig = (yHist[i] == 0) ? 0 : (float)Math.sqrt(sumErrorPerBin[i])/yHist[i];

            sumOrigSquared += yBinErrorOrig;
        }

        float contribFromAllToEachBin = (sumOrigSquared - sumAlt)/yHist.length;

        for (int i = 0; i < yHist.length; i++) {

            yHistErrorsOutput[i] += contribFromAllToEachBin;
        }
    }
    
    /**
     * determine the errors in determining the width of the histogram for points 
     * with y above yLimit.  This is meant to determine the error in 
     * calculations of things like fwhm.
     * 
     * @param xHist
     * @param yHist
     * @param xErrors
     * @param yErrors
     * @param yMaxFactor
     * @return
     */
    public static float calculateHistogramWidthYLimitError(float[] xHist, 
        float[] yHist, float[] xErrors, float[] yErrors, float yMaxFactor) {

        /* Errors in histogram:
         *     error in Y is sqrt(Y) and that is already in standard units.
         *     error in X is resolvability, which is bin size = (xHist[1] - xHist[0])/2.
         * 
         *                                | df |^2               | df |^2         df   df
         *      (sigma_f)^2 =  (sigma_x)^2|----|   +  (sigma_y)^2|----|    +  2 * -- * -- * cov_ab
         *                                | dx |                 | dy |           dx   dy
         * 
         *      For uncorrelated variables the covariance terms are zero.
         * 
         *      If f = XY, and X and Y are not correlated, we have:
         *          sigma^2  =  xError^2*(Y^2)  +  yError^2*(X^2) 
         *          
         *      For sum defined as a integrated area divided by Y:
         *      
         *                                       X_i*Y_i
         *          f = sum_over_i_to_yLimitIdx( ------- ) = sum_over_i_to_yLimitIdx( X_i )
         *                                         Y_i
         *       
         *          sigma^2  =  xError^2*(1)
         *          
         */
        int yPeakIdx = MiscMath0.findYMaxIndex(yHist);
        float yLimit = yMaxFactor * yHist[yPeakIdx];
        int yLimitIdx = -1;
        for (int i = 0; i < xHist.length; i++) {
            if (i > yPeakIdx) {
                if (yHist[i] > yLimit) {
                    yLimitIdx = i;
                } else {
                    break;
                }
            } else {
                yLimitIdx = i;
            }
        }
        
        float sum = 0.0f;
        for (int i = 0; i <= yLimitIdx; i++) {
            float xe = xErrors[i];
            sum += (xe * xe);
        }

        sum = (float) Math.sqrt(sum);

        return sum;
    }

    public static HistogramHolder defaultHistogramCreator(float[] values, 
        float[] valueErrors) {
        
        if (values == null || valueErrors == null || values.length != valueErrors.length) {
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }
        
        if (values.length < 15) {
            float[] minMax = MiscMath0.calculateOuterRoundedMinAndMax(values);
            return calculateSturgesHistogram(minMax[0], minMax[1], values, 
                valueErrors);
        } else if (values.length < 100) {
            return createSimpleHistogram(values, valueErrors);
        }
        
        HistogramHolder hist = calculateSturgesHistogramRemoveZeroTail(values, 
            valueErrors);
                
        return hist;
    }
    
    /**
     * binWidth = 2*IQR * n^(−1/3)
     * @param values
     * @return 
     */
    public static HistogramHolder calculateFreedmanDiaconisHistogram(
        double[] values) {
        
        // binWidth = 2*IQR * n^(−1/3)
        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(values);
        
        float binWidth = (float)(2.*medianAndIQR[1]*Math.pow(values.length, -1./3.));
        
        float[] vf = new float[values.length];
        for (int i = 0; i < vf.length; ++i) {
            vf[i] = (float)values[i];
        }
       
        return createSimpleHistogram(binWidth, vf);
    }
    
    /**
     * warning: casts all to float
     * binWidth = 3.49 * stDev * n^(−1/3)
     * @param values
     * @return 
     */
    public static HistogramHolder calculateScottsHistogram(
        double[] values, double min, double max) {
        
        // binWidth = 3.49 * stDev * n^(−1/3)
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(values);
        
        float binWidth = (float)(3.39 * avgAndStDev[1] * Math.pow(values.length, -1./3.));
        
        float[] vf = new float[values.length];
        for (int i = 0; i < vf.length; ++i) {
            vf[i] = (float)values[i];
        }
        
        return createSimpleHistogram(binWidth, vf, (float)min, (float)max);
    }
    
    public static HistogramHolder calculateScottsHistogram(
        double[] values) {
        
        // binWidth = 3.49 * stDev * n^(−1/3)
        double[] avgAndStDev = MiscMath0.getAvgAndStDev(values);
        
        float binWidth = (float)(3.39 * avgAndStDev[1] * Math.pow(values.length, -1./3.));
        
        float[] vf = new float[values.length];
        for (int i = 0; i < vf.length; ++i) {
            vf[i] = (float)values[i];
        }
        
        return createSimpleHistogram(binWidth, vf);
    }
    
    /**
     * nBins = log_2(n) + 1:
     * @param xMin
     * @param xMax
     * @param values
     * @param valueErrors
     * @return 
     */
    public static HistogramHolder calculateSturgesHistogram(
        final float xMin, final float xMax,
        float[] values, float[] valueErrors) {
    
        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }

        /*
        nBins = log_2(n) + 1:  note log2(x) = Math.log(x)/Math.log(2.)
        ----------------------
        
        log_2(values.length) + 1 = (xMax - xMin)/binWidth
        
        ==> binWidth = (xMax - xMin) / (log_2(values.length) + 1)
        */
        
        float binWidth = (float) ((xMax - xMin)/((Math.log(values.length)/Math.log(2.)) + 1));
        
        int nBins = (int)((Math.log(values.length)/Math.log(2.)) + 1);
        
        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
       
        Histogram.createHistogram(values, nBins, xMin, xMax, xHist, yHist, 
            binWidth);
        
        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, 
            yErrors);
        
        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }
    
    public static HistogramHolder calculateSturgesHistogramRemoveZeroTail(
        float[] values, float[] valueErrors) {
    
        if (values == null || valueErrors == null || 
            values.length != valueErrors.length) {
            
            throw new IllegalArgumentException(
                "values and valueErrors cannot be null and must be the same length");
        }

        int nIntervalsSturges = (int)Math.ceil( Math.log(values.length)/Math.log(2));
        
        //int nItervalsRice = (int)(2*Math.pow(values.length, 0.3333));
        
        int nBins = 25;
        
        if (values.length > 10000) {
            nBins = 40;
        }

        nBins = Math.max(nIntervalsSturges, nBins);
        
        float[] xHist = new float[nBins];
        int[] yHist = new int[nBins];
       
        float minx = MiscMath0.findMin(values);
        float maxx = MiscMath0.findMax(values);

        float binWidth = calculateBinWidth(minx, maxx, nBins);

        Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
        
        float maxy = MiscMath0.findMax(yHist);

        int minCountsLimit = (int)Math.max(5, 0.03f*maxy);
        int countsBelowMinAtTail = 0;
        int lastLowCountIdx = yHist.length - 1;
        for (int i = (yHist.length - 1); i > -1; i--) {
            if (yHist[i] < minCountsLimit) {
                countsBelowMinAtTail++;
                lastLowCountIdx = i;
            } else {
                break;
            }
        }
        
        if (countsBelowMinAtTail > 0) {
            
            maxx = xHist[lastLowCountIdx];
            
            // keep nbins the same?            
            binWidth = calculateBinWidth(minx, maxx, nBins);

            Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
            
            if (countsBelowMinAtTail > (nBins >> 1)) {
                // one more round of trimming
                maxy = MiscMath0.findMax(yHist);
                minCountsLimit = (int)Math.max(5, 0.03f*maxy);
                countsBelowMinAtTail = 0;
                lastLowCountIdx = yHist.length - 1;
                for (int i = (yHist.length - 1); i > -1; i--) {
                    if (yHist[i] < minCountsLimit) {
                        countsBelowMinAtTail++;
                        lastLowCountIdx = i;
                    } else {
                        break;
                    }
                }
                if (countsBelowMinAtTail > 0) {
                    maxx = xHist[lastLowCountIdx];
                    binWidth = calculateBinWidth(minx, maxx, nBins);
                    Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
                }
            }
        }
        
        if (values.length > 100) {
            // if there are a large number of points, we'd like to increase the 
            // resolution of the peak if needed
            int nLeftOfPeak = MiscMath0.findYMaxIndex(yHist);
            int nIter = 0;
            while (nIter < 30 && nLeftOfPeak < 3 && (yHist[nLeftOfPeak] > 100)) {
                binWidth *= 0.8f;
                Histogram.createHistogram(values, nBins, minx, maxx, xHist, yHist, binWidth);
                nLeftOfPeak = MiscMath0.findYMaxIndex(yHist);
                nIter++;
            }
        }
          
        float[] yHistFloat = new float[yHist.length];
        for (int i = 0; i < yHist.length; i++) {
            yHistFloat[i] = (float) yHist[i];
        }

        float[] yErrors = new float[xHist.length];
        float[] xErrors = new float[xHist.length];

        calulateHistogramBinErrors(xHist, yHist, values, valueErrors, xErrors, yErrors);

        HistogramHolder histogram = new HistogramHolder();
        histogram.setXHist(xHist);
        histogram.setYHist(yHist);
        histogram.setYHistFloat(yHistFloat);
        histogram.setYErrors(yErrors);
        histogram.setXErrors(xErrors);
        
        return histogram;
    }
    
    /**
     * if there is more than one peak in the histogram, reduce the histogram
     * to only that peak, else leave unaltered.
     * 
     * @param hist
     * @param values
     * @param valueErrors
     * @return 
     */
    public static HistogramHolder reduceHistogramToFirstPeak(HistogramHolder 
        hist, float[] values, float[] valueErrors) {
        
        int yPeakIdx = findFirstPeakIndex(hist);
        
        if (yPeakIdx == -1) {
            return hist;
        }
      
        int yMinPeakIdx = findFirstMinimaFollowingPeak(hist, yPeakIdx);
        
        if (yMinPeakIdx == -1) {
            return hist;
        }
        
        int n = yMinPeakIdx + 1;
        
        HistogramHolder tmp = new HistogramHolder();
        tmp.setXHist(Arrays.copyOfRange(hist.getXHist(), 0, n));
        tmp.setYHistFloat(Arrays.copyOfRange(hist.getYHistFloat(), 0, n));
        tmp.setXErrors(Arrays.copyOfRange(hist.getXErrors(), 0, n));
        tmp.setYErrors(Arrays.copyOfRange(hist.getYErrors(), 0, n));
        
        return tmp;
    }

    public static int findFirstPeakIndex(HistogramHolder hist) {
        
        float yPeak = Float.NEGATIVE_INFINITY;
        int yPeakIdx = -1;
        
        // specific to use here, find max within first half of histogram
        for (int i = 0; i < (hist.getXHist().length >> 1); i++) {
            
            float y = hist.getYHistFloat()[i];
            
            if (y > yPeak) {
                yPeak = y;
                yPeakIdx = i;
            }
        }
        
        return yPeakIdx;
    }

    public static int findFirstMinimaFollowingPeak(HistogramHolder hist, 
        int yPeakIdx) {
    
        //TODO:  could be improved to smooth over noise to find true minimum
        
        float yPeakMinimum = Float.MAX_VALUE;
        int yPeakMinIdx = -1;
        
        // find min within first half of histogram, after peak
        int n = (int)(0.5f * hist.getXHist().length);
        
        if ((n - yPeakIdx) < 3) {
            n = hist.getXHist().length;
        }
        
        for (int i = (yPeakIdx + 1); i < n; i++) {
            
            float y = hist.getYHistFloat()[i];
            
            if (y < yPeakMinimum) {
                yPeakMinimum = y;
                yPeakMinIdx = i;
            }
        }
        
        return yPeakMinIdx;
    }

    /**
     * 
     * @param indexes
     * @param imageValues these should be non negative numbers
     * @return 
     */
    public static PairIntArray createADescendingSortByKeyArray(
        TIntSet indexes, int[] imageValues) {
        
        TIntIterator iter = indexes.iterator();
        TIntIntMap freqMap = new TIntIntHashMap();
        while (iter.hasNext()) {
            int idx = iter.next();
            int v = imageValues[idx];
            if (freqMap.containsKey(v)) {
                freqMap.put(v, freqMap.get(v) + 1);
            } else {
                freqMap.put(v, 1);
            }
        }
       
        int[] v = new int[freqMap.size()];
        int[] c = new int[freqMap.size()];
                
        int vMax = Integer.MIN_VALUE;
        
        TIntIntIterator iter2 = freqMap.iterator();
        for (int i = 0; i < freqMap.size(); ++i) {
            iter2.advance();          
            v[i] = iter2.key();
            c[i] = iter2.value();
            if (v[i] > vMax) {
                vMax = v[i];
            } 
        }
        
        double nlg2n = v.length * Math.log(v.length)/Math.log(2);
        
        if (nlg2n < vMax || (vMax > 100000)) {
            MiscSorter.sortByDecr(v, c);
        } else if (v.length > 1) {
            CountingSort.sortByDecr(v, c, vMax);
        }
        
        PairIntArray p = new PairIntArray();
        
        for (int i = 0; i < c.length; i++) {
            p.add(v[i], c[i]);
        }
        
        return p;
    }
    
    public static float measureFWHMOfStrongestPeak(HistogramHolder hist) {
        
        if (hist == null) {
            throw new IllegalArgumentException("hist cannot be null");
        }
        
        int yMaxIdx = MiscMath0.findYMaxIndex(hist.getYHist());
        
        return measureFWHM(hist, yMaxIdx);
    }
    
    public static float[] measureFWHMOfAllPeaks(HistogramHolder hist, float frac) {
        
        if (hist == null) {
            throw new IllegalArgumentException("hist cannot be null");
        }
        
        List<Integer> yPeakIndexes = MiscMath0.findStrongPeakIndexes(hist, frac);

        float[] fwhms = new float[yPeakIndexes.size()];
        
        for (int i = 0; i < fwhms.length; i++) {
            
            int yPeakIdx = yPeakIndexes.get(i).intValue();
            
            fwhms[i] = measureFWHM(hist, yPeakIdx);
        }
        
        return fwhms;
    }
    
    public static float measureFWHM(HistogramHolder hist, int yPeakIndex) {
        
        if (hist == null) {
            throw new IllegalArgumentException("hist cannot be null");
        }
                
        if ((yPeakIndex == -1) || (yPeakIndex > (hist.getXHist().length - 1))) {
            return 0;
        }
        
        float halfPeak = hist.getYHist()[yPeakIndex]/2.f;
        
        float x0 = Float.NEGATIVE_INFINITY;
        if (yPeakIndex == 0) {
            x0 = hist.getXHist()[0];
        } else {
            for (int i = 0; i <= yPeakIndex; i++) {
                float y = hist.getYHistFloat()[i];
                if (y == halfPeak) {
                    x0 = hist.getXHist()[i];
                } else if (y > halfPeak) {
                    if (i == 0) {
                        x0 = hist.getXHist()[0];
                    } else {
                        // interpret between i and i-1
                        float dx01 = hist.getXHist()[i] - hist.getXHist()[i - 1];
                        float dy01 = y - hist.getYHistFloat()[i - 1];
                        float dy0h = y - halfPeak;
                        float ratio = dy0h/dy01;
                        x0 = hist.getXHist()[i] - (dx01 * ratio);
                    }
                    break;
                }
            }
        }
        
        float x1 = Float.NEGATIVE_INFINITY;
        if (yPeakIndex == (hist.getYHist().length - 1)) {
            x1 = hist.getXHist()[(hist.getYHist().length - 1)];
        } else {
            for (int i = (yPeakIndex + 1); i < hist.getYHist().length; i++) {    
                float y = hist.getYHist()[i];
                if (y == halfPeak) {
                    x1 = hist.getXHist()[i];
                } else if (y < halfPeak) {
                    if (i == (hist.getYHist().length - 1)) {
                        x1 = hist.getXHist()[(hist.getYHist().length - 1)];
                    } else {
                        // interpret between i and i-1
                        float dx01 = hist.getXHist()[i] - hist.getXHist()[i - 1];
                        float dy01 = y - hist.getYHistFloat()[i - 1];
                        float dy0h = y - halfPeak;
                        float ratio = dy0h/dy01;
                        x1 = hist.getXHist()[i] - (dx01 * ratio);
                    }
                    break;
                }
            }
        }
        
        if ((x0 == Float.NEGATIVE_INFINITY) || (x1 == Float.NEGATIVE_INFINITY)) {
            return 0;
        }
        
        float fwhm = x1 - x0;
            
        return fwhm;
    }
    
    /**
     * 
     * @param pointValues
     * @param minBin first bin's pixel value, inclusive
     * @param maxBin last bin's pixel value, inclusive.
     * @param nBins
     * @return 
     */
    public static int[] createHistogram(Map<PairInt, Integer> pointValues,
        int minBin, int maxBin, int nBins) {
        
        int[] h = new int[nBins];
        
        // (255 - 0 + 1)/256
        int binWidth = (maxBin - minBin + 1)/nBins;
        
        for (Entry<PairInt, Integer> entry : pointValues.entrySet()) {
            
            int v = entry.getValue().intValue();
            
            int binNumber = (v - minBin)/binWidth;
            
            //assert(binNumber >= 0);
            //assert(binNumber < maxBin);
            
            h[binNumber]++;
        }
        
        return h;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findMax(int[] a) {
        int max = Integer.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }

    public static int findMin(int[] a) {
        int min = Integer.MAX_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i] < min) {
                min = a[i];
            }
        }
        return min;
    }
}
