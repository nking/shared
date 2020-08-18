
package algorithms.misc;

import algorithms.matrix.MatrixUtil;
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
     * performs "standard unit normalization" on the points, that is,
     * re-scales data to have a mean of 0 and a standard deviation of 1 (unit variance)
     * for each dimension.
     * the input data must have format [nSamples][nVariables], e.g.
     * data[0] = [10, 100, 1000];
     * data[1] = [9,  101, 999];
     *   for nSamples=2 and nVariables = 3;
     * 
     * @param data the input data must have format [nSamples][nVariables], e.g.
     * data[0] = [10, 100, 1000];
     * data[1] = [9,  101, 999];
     *   for nSamples=2 and nVariables = 3;
     * @param outputMean the mean of the input data per variable (a.k.a. dimension)
     * @param outputStandardDeviation the standard deviation of the mean subtracted
     * data per dimension.
     * @return an array holding (data - mean)/stdev
     */
    public static double[][] standardUnitNormalization(double[][] data, 
        double[] outputMean, double[] outputStandardDeviation) {
        
        int nSamples = data.length;
        int nVars = data[0].length;
        
        if (outputMean.length != nVars) {
            throw new IllegalArgumentException("outputMean.length must equal data[0].length");
        }
        if (outputStandardDeviation.length != nVars) {
            throw new IllegalArgumentException("outputStandardDeviation.length must equal data[0].length");
        }
        
        // subtract the mean from the data,
        // then calculate the standard deviation and divide the data by it.
        double[][] out = MatrixUtil.copy(data);
        
        double[] c = MatrixUtil.mean(out);
        System.arraycopy(c, 0, outputMean, 0, c.length);
                
        int i, d;
        for (i = 0; i < out.length; ++i) {
            for (d = 0; d < out[i].length; ++d) {
                out[i][d] -= c[d];
            }
        }
        
        Arrays.fill(outputStandardDeviation, 0);
                
        for (i = 0; i < nSamples; ++i) {
            for (d = 0; d < nVars; ++d) {
                //mean has already been subtracted:
                outputStandardDeviation[d] += (out[i][d]*out[i][d]);
            }
        }
        for (d = 0; d < nVars; ++d) {
            outputStandardDeviation[d] = Math.sqrt(outputStandardDeviation[d]/(nSamples - 1.0)); 
        }
        
        for (i = 0; i < nSamples; ++i) {
            for (d = 0; d < nVars; ++d) {
                if (outputStandardDeviation[d] > 0.) {
                    out[i][d] /= outputStandardDeviation[d];
                }
            }
        }
        
        return out;
    }
}
