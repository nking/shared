
package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import java.util.Arrays;

/**
 * various methods for performing statistical standardization and de-normalization
 * on a data.
 * 
 * see https://en.wikipedia.org/wiki/Normalization_(statistics) for more.
 * 
 * TODO: add a non-isotropic scaling (see Hartley 1997, end of Section 5):
     transform points so that
     1) Their centroid is at the origin.
     2) The principal moments are both equal to unity
  
 * @author nichole
 */
public class Standardization {
    
    /**
     * performs "standard unit normalization" on the points, that is,
     * transforms the data to have a mean of 0 and a standard deviation of 1 (unit variance)
     * for each dimension.
     * 
     @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     @param nDimensions the number of dimensions of a point
     @param outputMean the mean of the input data per dimension
     @param outputStandardDeviation the standard deviation of the mean subtracted
     * data per dimension.
     @return an array holding the standardized data in format of data array, 
     * e.g. point_0_0, point_0_1, point_0_2...point_0_{nDimensions-1}].
     * Note that the results should have a mean of 0 and a standard deviation
     * of 1 or 0.
     */
    public static double[] standardUnitNormalization(double[] data, 
        int nDimensions, double[] outputMean, double[] outputStandardDeviation) {
        
        if (outputMean.length != nDimensions) {
            throw new IllegalArgumentException("outputMean.length must equal nDimensions");
        }
        if (outputStandardDeviation.length != nDimensions) {
            throw new IllegalArgumentException("outputStandardDeviation.length must equal nDimensions");
        }
        
        // subtract the mean from the data,
        // then calculate the standard deviation and divide the data by it.
        double[] out = Arrays.copyOf(data, data.length);
        
        double[] c = MiscMath0.mean(data, nDimensions);
        System.arraycopy(c, 0, outputMean, 0, c.length);
        
        int nData = data.length/nDimensions;
        
        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                out[j] -= c[d];
            }
        }
        
        double[][] mAndSD = MiscMath0.standardDeviation(out, nDimensions);
        assert(mAndSD[0].length == nDimensions);
        assert(mAndSD[1].length == nDimensions);
        
        System.arraycopy(mAndSD[1], 0, outputStandardDeviation, 0, nDimensions);
        
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                if (mAndSD[1][d] > 0.) {
                    out[j] /= mAndSD[1][d];
                }
            }
        }
        
        return out;
    }
    
    /**
     * performs de-normalization per dimension on data that has been
     * standardized with "standard unit normalization", that is,
     * multiplies by the given standardDeviation and adds the given mean.
     * for each dimension.
     * 
     @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     @param nDimensions the number of dimensions of a point
     @param mean the mean of the unnormalized input data per dimension
     @param standardDeviation the standard deviation of the mean subtracted
     * unnormalized data per dimension.
     @return an array holding the de-normalized data in format similar to data array.
     */
    public static double[] standardUnitDenormalization(final double[] data, 
        final int nDimensions, final double[] mean, final double[] standardDeviation) {
        
        if (mean.length != nDimensions) {
            throw new IllegalArgumentException("mean.length must equal nDimensions");
        }
        if (standardDeviation.length != nDimensions) {
            throw new IllegalArgumentException("standardDeviation.length must equal nDimensions");
        }
        
        double[] out = Arrays.copyOf(data, data.length);
        
        int nData = data.length/nDimensions;
        
        int i, j, d;
        for (i = 0; i < nData; ++i) {
            for (d = 0; d < nDimensions; ++d) {
                j = i * nDimensions + d;
                if (standardDeviation[d] > 0){
                    out[j] *= standardDeviation[d];
                }
                out[j] += mean[d];
            }
        }
        
        return out;
    }

    /**
     * unit standard normalize data and return the mean and standard deviation of the zero-centered data.
     * @param im
     * @return the mean and standard deviation of zero-centered data
     */
    public static double[] unitStandardNormalizeImage(double[][] im) {
        double mean = 0;
        int np = 0;
        for (double[] rows : im) {
            for (double val : rows) {
                mean += val;
                ++np;
            }
        }
        mean /= (double)(np);

        double std = 0;
        int col;
        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                im[row][col] -= mean;
                // stdev of 0 centered mean.  diff from 0:
                std += (im[row][col] * im[row][col]);
            }
        }
        // remove 1 degree of freedom for having to calculate the mean:
        std = Math.sqrt(std/(np - 1.));

        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                im[row][col] /= std;
            }
        }
        return new double[]{mean, std};
    }

    /**
     * apply min max normalization to image and return the min and max.
     * The operation changes the data values so that they are all in range [0, 1].
     * @param im
     * @return the mean and standard deviation of zero-centered data
     */
    public static double[] minMaxNormalizeImage(double[][] im) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        int col;
        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                min = Math.min(min, im[row][col]);
                max = Math.max(max, im[row][col]);
            }
        }
        double range = max - min;
        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                im[row][col] = (im[row][col] - min)/range;
            }
        }
        return new double[]{min, max};
    }

    /**
     * apply min max normalization to image, given the min and max.  useful if need to
     * apply same normalization to different dataset
     * @param im
     * @return the mean and standard deviation of zero-centered data
     */
    public static void minMaxNormalizeImage(double[][] im, double min, double max) {
        double range = max - min;
        int col;
        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                im[row][col] = (im[row][col] - min)/range;
            }
        }
    }


    /**
     * subtract mean from data then divide by std.  used to perform same normalization on different images.
     * @param im
     * @param the mean to subtract from data
     * @param the standard deviation to divide data by
     */
    public static void applyMeanAndStdToImage(double[][] im, double mean, double std) {
        int col;
        for (int row = 0; row < im.length; ++row) {
            for (col = 0; col < im[row].length; ++col) {
                im[row][col] -= mean;
                im[row][col] /= std;
            }
        }
    }
    
     /**
     * performs "standard unit normalization" on the points, that is,
     * transforms the data to have a mean of 0 and a standard deviation of 1 (unit variance)
     * for each dimension.
     * the input data must have format [nSamples][nVariables], e.g.
     * data[0] = [10, 100, 1000];
     * data[1] = [9,  101, 999];
     *   for nSamples=2 and nVariables = 3;
     * 
     @param data the input data must have format [nSamples][nVariables], e.g.
     * data[0] = [10, 100, 1000];
     * data[1] = [9,  101, 999];
     *   for nSamples=2 and nVariables = 3;
     @param outputMean the mean of the input data per variable (a.k.a. dimension)
     @param outputStandardDeviation the standard deviation of the mean subtracted
     * data per dimension.
     @return an array holding (data - mean)/stdev
     */
    public static double[][] standardUnitNormalization(double[][] data, 
        double[] outputMean, double[] outputStandardDeviation) {
        
        int m = data.length;
        int n = data[0].length;
        
        if (outputMean.length != n) {
            throw new IllegalArgumentException("outputMean.length must equal data[0].length");
        }
        if (outputStandardDeviation.length != n) {
            throw new IllegalArgumentException("outputStandardDeviation.length must equal data[0].length");
        }
        
        // subtract the mean from the data,
        // then calculate the standard deviation and divide the data by it.
        double[][] out = MatrixUtil.copy(data);
        
        double[] mean = MatrixUtil.columnMeans(out);
        System.arraycopy(mean, 0, outputMean, 0, mean.length);

        int i;
        int j;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                out[i][j] = data[i][j] - mean[j];
            }
        }
        
        Arrays.fill(outputStandardDeviation, 0);

        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                //mean has already been subtracted:
                outputStandardDeviation[j] += (out[i][j]*out[i][j]);
            }
        }
        for (j = 0; j < n; ++j) {
            outputStandardDeviation[j] = Math.sqrt(outputStandardDeviation[j]/(m - 1.0));
        }

        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                if (outputStandardDeviation[j] > 0.) {
                    out[i][j] /= outputStandardDeviation[j];
                }
            }
        }
        
        return out;
    }

    /**
     * calculate the mean of each column of data and subtract it from each value in the respective column of data.
     @param data
     @return
     */
    public static double[][] zeroCenterMean(double[][] data) {
        int m = data.length;
        int n = data[0].length;

        double[][] out = MatrixUtil.zeros(m, n);

        double[] mean = MatrixUtil.columnMeans(data);

        int i;
        int j;
        for (j = 0; j < n; ++j) {
            for (i = 0; i < m; ++i) {
                out[i][j] = data[i][j] - mean[j];
            }
        }
        return out;
    }
}
