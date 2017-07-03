package algorithms.misc;

import algorithms.CountingSort;
import algorithms.util.PairInt;
import algorithms.util.PairIntArray;
import java.util.List;
import java.util.ArrayList;
import java.util.Collection;
import gnu.trove.set.TIntSet;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.list.TIntList;

/**
    miscellaneous math methods. some could probably be improved.

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
public class MiscMath0 {

    public static float[] calculateOuterRoundedMinAndMax(float[] a) {

        // find the powers of 10 for the data min and max
        float xmin = Float.MAX_VALUE;
        float xmax = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > xmax) {
                xmax = a[i];
            }
            if (a[i] < xmin) {
                xmin = a[i];
            }
        }

        xmax = MiscMath0.roundUpByLargestPower(xmax);

        xmin = MiscMath0.roundDownByLargestPower(xmin);

        // xmax > 1 and xmin is between 0 and 1, round xmin down
        if ((xmax > 1) && (xmin > 0) && (xmin < 1.0)) {
            xmin = 0;
        }

        return new float[]{xmin, xmax};
    }

    /**
     Round number up to the next digit in the largest power.

        MiscMath0.roundUpByLargestPower(31.1f) == 40.0f;

        MiscMath0.roundUpByLargestPower(0.11f) == 1.0f;

        MiscMath0.roundUpByLargestPower(-0.011f) == -0.02f;

        MiscMath0.roundUpByLargestPower(310.1f) == 400.0f;

        MiscMath0.roundUpByLargestPower(-3.1) == -4.0f;

        MiscMath0.roundUpByLargestPower(3.1) == 4.0f;

    @param f
    @return
    */
    public static float roundUpByLargestPower(float f) {

        int pow = findPowerOf10(f);

        double pow10 = Math.pow(10, pow);

        int d = (int)(f/pow10);
        int m;
        if (f == pow10) {
            m = d;
        } else if (f > 0) {
            // residual ?
            if (f > (d*pow10)) {
                m = d + 1;
            } else {
                m = d;
            }
        } else {
            // decimals
            m = d - 1;
        }

        float r = (float) (m * pow10);

        return r;
    }

    /**
     * Round number down to in the largest power.
     * For example,
     *     roundDownByLargestPower(3.1) returns 3.0
     *     roundDownByLargestPower(-3.1) returns -4.0
     *     roundDownByLargestPower(31.1) returns 30.0
     *
     * @param f
     * @return
     */
    public static float roundDownByLargestPower(float f) {

        if (f == 0) {
            return 0;
        }

        int pow = findPowerOf10(f);
        double pow10 = Math.pow(10, pow);

        double r;
        if (f > 0) {
            double m = f % pow10;
            r = f - m;
        } else {
            int m = (int)(f/pow10);
            r = pow10 * (m - 1);
        }

        return (float)r;
    }
    
    public static int findPowerOf10(float a) {

        if (a == 0) {
            return 0;
        }
        if (a < 0.f) {
            a *= -1.0f;
        }
        double b = Math.log10(a);
        if (b > 0) {
            return (int)b;
        } else {
            if (b >= 1) {
                return (int)Math.round(b);
            } else if (b > -1) {
                // fractions between -1 and +1
                 return findPowerOf10_2(a);
            } else {
                return (int)Math.round(b);
            }
        }
    }

    public static int findPowerOf10_2(float a) {

        if (a == 0) {
            return 0;
        }

        int power = 0;

        if (a <= 1.0f) {
            while (a < 1.0) {
                a *=  10.0;
                power--;
            }
        } else {
            // precision errors in multiplication here are trouble for non base2 numbers such as powers of 10
            while (a >= 1.0) {
                a /= 10.0;
                power++;
            }
            power--;
        }

        return power;
    }

    /**
     * find max but ignore values such as FLOAT.MAX_VALUE, infinity, and NAN
     * @param a
     * @return
     */
    public static int findYMaxIndex(float[] a) {
        if (a == null || a.length == 0) {
            return -1;
        }
        float max = Float.NEGATIVE_INFINITY;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] > max) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i]) 
                && (a[i] < Float.MAX_VALUE)) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(int[] a) {
        if (a == null || a.length == 0) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(TIntList a) {
        if (a == null || a.size() == 0) {
            return -1;
        }
        int max = Integer.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.size(); i++) {
            if (a.get(i) > max) {
                max = a.get(i);
                index = i;
            }
        }
        return index;
    }

    public static int findYMaxIndex(long[] a) {
    
        if (a == null || a.length == 0) {
            return -1;
        }
        long max = Long.MIN_VALUE;
        int index = 0;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
                index = i;
            }
        }
        return index;
    }
 
    public static float findMin(float[] a) {
        float min = Float.POSITIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if ((a[i] < min) && !Float.isInfinite(a[i]) && !Float.isNaN(a[i])) {
                min = a[i];
            }
        }
        return min;
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
     * find max 
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

    /**
     * find max 
     * @param a
     * @return
     */
    public static float findMax(float[] a) {
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < a.length; i++) {
            if (a[i] > max) {
                max = a[i];
            }
        }
        return max;
    }

    /**
     * for the given histogram, returns the indexes of the primary peak and
     * any peaks which are larger than frac*maxPeak above their neighboring
     * values.
     * @param h
     * @param frac the fraction of the y value of the maximum peak that is used
     * as a critical limit that any other peaks must have in excess of their
     * neighboring points.  For example, 0.1.
     * @return 
     */
    public static List<Integer> findStrongPeakIndexes(HistogramHolder h, float frac) {
        
        float[] x = h.getXHist();
        int[] y = h.getYHist();
        
        if (x == null || y == null || x.length == 0 || y.length == 0) {
            return null;
        }
        
        if (x.length == 1) {
            List<Integer> minMaxIndexes = new ArrayList<Integer>();
            if (y[0] > 0) {
                minMaxIndexes.add(Integer.valueOf(0));
            }
            return minMaxIndexes;
        }
        
        int yPeakIdx = MiscMath0.findYMaxIndex(y);
        
        if (yPeakIdx == -1) {
            return null;
        }
        
        int yMaxPeak = y[yPeakIdx];
        
        /*
        storing the minima and maxima in the same array list.
        the minima have -1*index within k
        and the maxima keep their positive values of the index within k.
        */
        List<Integer> minMaxIndexes = new ArrayList<Integer>();
        
        float lastY = y[0];
        boolean incr = true;
        for (int ii = 1; ii < y.length; ii++) {
            if ((y[ii] < lastY) && incr) {
                minMaxIndexes.add(Integer.valueOf(ii - 1));
                incr = false;
            } else if ((y[ii] > lastY) && !incr) {
                minMaxIndexes.add(Integer.valueOf(-1*(ii - 1)));
                incr = true;
            }
            lastY = y[ii];
        }
        
        if (incr) {
            // add the last point
             minMaxIndexes.add(Integer.valueOf(y.length - 1));
        }
        
        // for the histograms of euclidean point combination differences,
        // this should usually be singly peaked
        if (minMaxIndexes.size() == 1) {
            return minMaxIndexes;
        }
        
        float limit = frac * yMaxPeak;
        
        // find peaks where y[ii] is > limit above adjacent local minima
        
        List<Integer> peaks = new ArrayList<Integer>();

        for (int ii = 0; ii < minMaxIndexes.size(); ii++) {

            int idx = minMaxIndexes.get(ii).intValue();

            if (idx > -1) {
                // this is a maxima
                
                boolean found = false;
                
                // compare to preceding minimum
                for (int iii = (ii - 1); iii > -1; iii--) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                            found = true;
                        }
                        break;
                    }
                }
                if (found) {
                    continue;
                }

                //compare to proceeding minimum
                for (int iii = (ii + 1); iii < minMaxIndexes.size(); iii++) {
                    int idx2 = minMaxIndexes.get(iii).intValue();
                    if (idx2 < 0) {
                        float compare = y[-1*idx2];
                        float diff = y[idx] - compare;
                        if (diff >= limit) {
                            peaks.add(Integer.valueOf(idx));
                        }
                        break;
                    }
                }
            }
        }
       
        return peaks;
    }

    public static List<Integer> findStrongPeakIndexesDescSort(
        HistogramHolder hist, float fracMax) {

        List<Integer> indexes = findStrongPeakIndexes(hist, fracMax);
        
        if (indexes == null || indexes.size() < 2) {
            return indexes;
        }
        
        int[] idxs = new int[indexes.size()];
        int[] c = new int[idxs.length];
        
        int maxC = Integer.MIN_VALUE;
        for (int i = 0; i < indexes.size(); ++i) {
            idxs[i] = indexes.get(i).intValue();
            c[i] = hist.getYHist()[idxs[i]];
            if (c[i] > maxC) {
                maxC = c[i];
            }
        }
        
        CountingSort.sortByDecr(c, idxs, maxC);
        
        indexes.clear();
        
        for (int i = 0; i < idxs.length; ++i) {
            indexes.add(Integer.valueOf(idxs[i]));
        }
        
        return indexes;
    }

    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x) {
        
        return getAvgAndStDev(x, x.length);
    }

    /**
     * given an array of points, return the average and standard deviation from
     * the average
     * @param x
     * @return 
     */
    public static double[] getAvgAndStDev(int[] x, int length) {
        
        long sumX = 0;
        for (int i = 0; i < length; i++) {
            sumX += x[i];
        }
        
        double avgX = (double)sumX/(double)length;
        
        sumX = 0;
        for (int i = 0; i < length; i++) {
            double diffX = x[i] - avgX;
            sumX += (diffX * diffX);
        }
        double stdDevX = Math.sqrt(sumX/(length - 1.0f));
        
        return new double[]{avgX, stdDevX};
    }    

    public static PairIntArray get20NeighborOffsets() {
        
        PairIntArray r2Offsets = new PairIntArray();
        
                              r2Offsets.add(-1, 2); r2Offsets.add(0, 2); r2Offsets.add(1, 2);
        r2Offsets.add(-2, 1); r2Offsets.add(-1, 1); r2Offsets.add(0, 1); r2Offsets.add(1, 1); r2Offsets.add(2, 1);
        r2Offsets.add(-2, 0); r2Offsets.add(-1, 0);                       r2Offsets.add(1, 0); r2Offsets.add(2, 0);
        r2Offsets.add(-2, -1); r2Offsets.add(-1, -1); r2Offsets.add(0, -1); r2Offsets.add(1, -1); r2Offsets.add(2, -1);
                              r2Offsets.add(-1, -2); r2Offsets.add(0, -2); r2Offsets.add(1, -2);
        
        return r2Offsets;
    }

     /**
     *
     * @param a
     * @return
     */
    public static int[] findMinMaxValues(int[][] a) {
        
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                int v = a[i][j];
                if (v < min) {
                    min = v;
                }
                if (v > max) {
                    max = v;
                }
            }
        }
        
        return new int[]{min, max};
    }

    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(Collection<PairInt> points) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (PairInt p : points) {
            int x = p.getX();
            int y = p.getY();
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }
    
    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(TIntSet pixelIdxs, int imgWidth) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        TIntIterator iter = pixelIdxs.iterator();
        
        while (iter.hasNext()) {
            
            int pixIdx = iter.next();
            int y = pixIdx/imgWidth;
            int x = pixIdx - (y * imgWidth);
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }

    /**
     * find the minima and maxima of x and y and return them as
     * int[]{xMin, xMax, yMin, yMax}
     * @param points
     * @return minMaxXY int[]{xMin, xMax, yMin, yMax}
     */
    public static int[] findMinMaxXY(PairIntArray points) {
        
        int xMin = Integer.MAX_VALUE;
        int xMax = Integer.MIN_VALUE;
        int yMin = Integer.MAX_VALUE;
        int yMax = Integer.MIN_VALUE;
        
        for (int i = 0; i < points.getN(); ++i) {
            int x = points.getX(i);
            int y = points.getY(i);
            if (x < xMin) {
                xMin = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y > yMax) {
                yMax = y;
            }
        }
        return new int[]{xMin, xMax, yMin, yMax};
    }
    
    public static boolean isAPowerOf2(int n) {
        // n XOR n-1
        return ((n == 0) || ((n & (n - 1)) == 0));
    }
    
    /**
     * determine the number of bits, that is, the msb position + 1.
     * Note that a value of 0 returns a bit length of 1.
     * @param v
     * @return 
     */
    public static int numberOfBits(int v) {
        
        if (v < 0) {
            v *= -1;
        } else if (v == 0) {
            return 1;
        }
        return 32 - Integer.numberOfLeadingZeros(v);
    }
    
    /**
     * 
     * @param v
     * @return 
     */
    public static int bitReverse(int v, int nBits) {

        int r = v;
                
        int rev = 0;
        for (int i = 0; i < nBits; i++) {
            rev = (rev << 1) | (r & 1);
            r >>= 1;
        }
        
        return rev;
        
    }
}
