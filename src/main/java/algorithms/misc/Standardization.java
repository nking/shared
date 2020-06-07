
package algorithms.misc;

import java.util.Arrays;

/**
 * various methods for performing statistical standardization and de-normalization
 * on a data.
 * 
 * see https://en.wikipedia.org/wiki/Normalization_(statistics) for more.
 * 
 * @author nichole
 */
public class Standardization {
    
    /**
     * performs "standard unit normalization" on the points, that is,
     * re-scales data to have a mean of 0 and a standard deviation of 1 (unit variance)
     * for each dimension.
     * 
     * @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     * @param nDimensions the number of dimensions of a point
     * @param outputMean the mean of the input data per dimension
     * @param outputStandardDeviation the standard deviation of the mean subtracted
     * data per dimension.
     * @return an array holding the standardized data in format of data array, 
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
        double diff;
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
     * @param data nDimensional data points in format [ point_0 in all dimensions,
     *   point_1 in all dimensions, ... point_{n-1} in all dimensions[
     * @param nDimensions the number of dimensions of a point
     * @param mean the mean of the unnormalized input data per dimension
     * @param standardDeviation the standard deviation of the mean subtracted
     * unnormalized data per dimension.
     * @return an array holding the de-normalized data in format similar to data array.
     */
    public static double[] standardUnitDenormalization(double[] data, 
        int nDimensions, double[] mean, double[] standardDeviation) {
        
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
}
