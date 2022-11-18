package algorithms.misc;

import algorithms.matrix.MatrixUtil;
import algorithms.sort.CountingSort;
import algorithms.util.Errors;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import gnu.trove.iterator.TIntIntIterator;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.map.TDoubleIntMap;
import gnu.trove.map.TIntDoubleMap;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TDoubleIntHashMap;
import gnu.trove.map.hash.TIntDoubleHashMap;
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
     * given data points and the number of bins to use, return the histogram as a 2 dimensional
     * array with row 0 being the bin centers, and row 1 being the histogram counts.
     * the range of data used are the maximum - minimum.
     * @param data data points
     * @param nBins the number of bins to use
     * @return 2 dimensional
     *      * array with row 0 being the bin centers, and row 1 being the histogram counts
     */
    public static double[][] createHistogram(double[] data, int nBins) {

        double[] minMax = MiscMath0.getMinMax(data);

        return createHistogram(data, nBins, minMax[1], minMax[0]);
    }

    /**
     * given data points and the number of bins to use, return the histogram as a 2 dimensional
     * array with row 0 being the bin centers, and row 1 being the histogram counts.
     * the range of data used are the maximum - minimum.
     * @param data data points
     * @param nBins the number of bins to use
     * @param max maximum data value that the histogram will hold.  this is the end of the last bin.
     * @param min minimum data value that thehistogram will hold.  this is the beginning of the first bin.
     * @return 2 dimensional
     *      * array with row 0 being the bin centers, and row 1 being the histogram counts
     */
    public static double[][] createHistogram(double[] data, int nBins, double max, double min) {

        int n = data.length;

        double binWidth = (max - min)/nBins;

        double[][] hist = new double[2][];
        hist[0] = new double[nBins];
        hist[1] = new double[nBins];

        for (int i = 0; i < nBins; i++) {
            hist[0][i] = min + binWidth*i + (binWidth/2.);
        }

        int bin;
        for (int i = 0; i < n; i++) {
            bin = (int) ((data[i] - min)/binWidth);
            if ((bin > -1) && (bin < nBins)) {
                hist[1][bin]++;
            }
        }

        return hist;
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
        System.out.println("nBins=" + nBins);
        
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
       TODO: revisit this.  have added wasserman statistics in other methods.
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
            CountingSort.sortByDecr(v, c);
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

    /**
     * calculate the risk of a histogram where risk = bias^2 + variance.
     * <pre>
     * risk estimator:
     *     JEst(h) = (2/((n-1)*h)) - ((n+1)/(n-1)) * sum over j=[0:n-1]( pEst_j^2 )
     *     where h is the bin width,
     *     n is the number of points used to generate the histogram,
     *     pEst_j is the histogram count in the jth bin divided by n.
     *
     *  Note that the algorithm expects that the data range of X was normalized
     *  to be between 0 and 1, inclusive.  Also, that h=1/m where m is pEst.length.
     *
     *  reference:
     *  Wasserman's "All of Statistics" eqn (20.14).
     *  and
     *  https://en.m.wikipedia.org/wiki/Histogram
     * </pre>
     * One can estimate the best binwidth h by minimizing the risk over
     * histograms generated with small to increasingly larger binwidths.
     * @param n the number of points used when constructing the histogram.
     * @param h the binwidth.
     * @param pEst array of the histogram counts, where each bin should have been divided by n.
     *             Note that the algorithm expects that the data range of X was normalized
     *      * to be between 0 and 1, inclusive.  Also, that h=1/m where m is pEst.length.
     * @return the estimated risk
     */
    public static double crossValidationRiskEstimator(int n, double h, double[] pEst) {
        int m = pEst.length;

        // assert that h == 1./m
        //assert(Math.abs(h - (1./m)) < 1e-2);

        double sum = 0;
        for (int j = 0; j < m; ++j) {
            sum += (pEst[j] * pEst[j]);
        }
        // from wikipedia, https://en.m.wikipedia.org/wiki/Histogram
        // reference to Stone 1984
        // AN ASYMPTOTICALLY OPTIMAL HISTOGRAM SELECTION RULE
        double t1 = 2./((n-1.)*h);
        //double t1 = 2./(n-1.);
        double t2 = (n+1)/(n*n*h*(n-1.));

        double jEst = t1 - t2*sum;
        return jEst;
    }

    /**
     * calculate the 1-alpha confidence envelope (a.k.a. confidence band) for the histogram.
     <pre>
     reference:
         Wasserman's "All of Statistics" eqn (20.17), (20.18).
     </pre>
     * @param n the number of points used when constructing the histogram.
     * @param h the binwidth. Note that the algorithm expects that the data range of X was normalized
     * to be between 0 and 1, inclusive.  Also, that h=1/m where m is pEst.length.
     * @param pEst array of is the histogram counts, where each bin should have been divided by n.
     * @param zAlpha the value from the z-table for the confidence level.  e.g. for 95% confidence level,
     *               zAlpha is 1.96.
     *  One can roughly estimate zAlpha using zAlpha = CDFStandardNormal.approxInverseShort(p).
     *  or use these:
     *    for 90%, zAlpha=1.645, for 95% zAalpha=1.96, for 98% zAlpha=2.326, for 99% zAlpha=2.576.
     * @return the lower and upper confidence envelope as 2 rows of length pEst (2 x pEst.length)
     */
    public static double[][] confidenceEnvelope(int n, double h, double[] pEst, double zAlpha) {
        int m = pEst.length;
        // assert that h == 1./m
        double c = ((zAlpha/(2.*m))/2.) * Math.sqrt((double)m/(double)n);

        double[][] lu = new double[2][m];
        lu[0] = new double[m];
        lu[1] = new double[m];
        int i;
        double fEstSR;
        for (i = 0; i < m; ++i) {
            fEstSR = Math.sqrt(pEst[i]/h);
            lu[0][i] = Math.max(fEstSR - c, 0);
            lu[0][i] *= lu[0][i];
            lu[1][i] = fEstSR + c;
            lu[1][i] *= lu[1][i];
        }
        return lu;
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

    /**
     * calculate the multidimensional histogram of the data using a simple brute force approach.
     * Only the populated bins are returned in a hash map with key = bin index, value = counts in the bin
     * where index is calculated in a telescoping manner.  e.g.:
     <pre>
      for 3D to one dimensional index:
          index = (((bin0Idx * n1) + bin1Idx) * n2) + bin2Idx
          where bin0Idx = x[row=some number][col=0] / h
          and bin1Idx = x[row=same row number][col=1] / h

     and for the same example, to extract the bin numbers from the one-dimensional index:
         bin2Idx = index % n1
         bin1Idx = index % (n1*n2)
         bin0Idx = index / (n1*n2)
     </pre>

     TODO: upon use, will create a way to read off the populated bins of a histogram as a slice
     of the returned sparse multidimensional histogram.  currently, the user has to iterate
     through the 1-D index to get the counts in the bins of interest.

     <pre>
     There is mention of choosing the bandwidth h in the lecture notes of
     Larry Wasserman, CMU
     36-708 Statistical Methods for Machine Learning by
     https://www.stat.cmu.edu/~larry/=sml/densityestimation.pdf
     after equation (8), but c2 isn't defined.
     One can estimate the bandwidth by using a range of h values to find the minimum result of
     crossValidationRiskEstimator() in this class, but the estimator is for 1-D.
     </pre>
     * @param h the binwidth to use.  note that using the same binwidth for all dimensions is not
     *          ideal, but it does allow for a n-Dimensional to 1-Dimensional bin index to be made
     *          in which the individual components can be extracted if needed... putting the results
     *          into a sparse form was a goal of this method.
     * @param x n x m data array where the n rows are the data points and the m columns are the
                dimensions. each column of data must have already been normalized to the range [0, 1].
                a 2-dimensional sample of 4 data points as an example:
                <pre>
                2 x 4 array:
                       dim 0    dim 1
                x[0] = (0,      0.5)
                x[1] = (0.25,   0.45)
                x[2] = (0.1,    0.55)
                x[3] = (0.35,   0.46)
                </pre>
     * @return
     */
    public static TIntIntMap createMultidimensionalHistogram(double h, double[][] x) {
        if (x.length == 0 || x[0].length == 0) {
            throw new IllegalArgumentException("x must have > 0 rows and > 0 columns");
        }
        int n0 = x[0].length;

        TIntIntMap hist = new TIntIntHashMap();
        int i;
        int d;
        int idx;
        int b;
        for (i = 0; i < x.length; ++i) {
            idx = 0;
            for (d = 0; d < x[i].length; ++d) {
                b = (int)(x[i][d]/h);
                //for 2D: (bin0 * n1) + bin1.  for 3D:  (((bin0 * n1) + bin1) * n2) + bin2...
                idx = (idx * n0) + b;
            }
            if (hist.containsKey(idx)) {
                hist.put(idx, hist.get(idx) + 1);
            } else {
                hist.put(idx, 1);
            }
        }
        return hist;
    }

    public static int numberOfBinsSqrt(int dataLength) {
        //https://en.m.wikipedia.org/wiki/Histogram
        return (int)Math.ceil(Math.sqrt(dataLength));
    }

    /**
     * assumes an approximately normal distribution.
     * for best results, dataLength >= 30.
     * reference: https://en.m.wikipedia.org/wiki/Histogram
     * @param dataLength
     * @return
     */
    public static int numberOfBinsSturges(int dataLength) {
        //https://en.m.wikipedia.org/wiki/Histogram
        return 1 + (int)Math.ceil(Math.log(dataLength)/Math.log(2));
    }
    /**
     * reference: https://en.m.wikipedia.org/wiki/Histogram
     * @param dataLength
     * @return
     */
    public static int numberOfBinsRice(int dataLength) {
        //https://en.m.wikipedia.org/wiki/Histogram
        return (int)Math.ceil(2.*Math.pow(dataLength, 1./3));
    }

    /**
     * improvement of sturges formula for non-normal data.
     * reference: https://en.m.wikipedia.org/wiki/Histogram
     * @return
     */
    public static int numberOfBinsDoane(double[] data) {

        double g1 = MiscMath0.calcSampleSkewness(data);

        int n = data.length;

        double sigmaG1 = Math.sqrt(  (6.*(n-2.))/((n+1.)*(n+3.)));

        double t2 = Math.log(n)/Math.log(2);
        double t3 = Math.log(1. + (Math.abs(g1)/sigmaG1))/Math.log(2);
        double nBins = 1 + t2 + t3;
        return (int)nBins;
    }

    /**
     * this is less sensitive to outliers in data because it uses IQR.
     * reference: https://en.m.wikipedia.org/wiki/Histogram
     * @return
     */
    public static int numberOfBinsFreedmanDiaconis(double[] data) {

        double[] medianAndIQR = MiscMath0.calcMedianAndIQR(data);
        double d = Math.pow(data.length, 1./3.);
        return (int) (2. * medianAndIQR[1]/d);
    }
}
