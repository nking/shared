package algorithms.signalProcessing;

import algorithms.util.OneDFloatArray;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform1D {
   
    /**
      The a-trous algorithm is a fast implementation of a wavelet transform 
       with no downsampling.   It is non-orthogonal, semi-linear runtime
       complexity, is invariant under translation, and the transform is 
       isotropic.
       Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
       The scaling function used is the lower resolution choice, the triangle
       function.
       <pre>
       The method uses recursive convolution operations, including previous
       result to make next.
       </pre>
     @param input
     @param outputTransformed
     @param outputCoeff 
     */
    public void calculateWithTriangleScalingFunction(float[] input,
        List<OneDFloatArray> outputTransformed, List<OneDFloatArray> outputCoeff) {
        
        int imgDimen = input.length;

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        TriangleFunction1D scalingFunction = new TriangleFunction1D();
        
        OneDFloatArray cp = new OneDFloatArray(Arrays.copyOf(input, input.length));
        outputTransformed.add(cp);
        
        OneDFloatArray cf = new OneDFloatArray(new float[input.length]);
        outputCoeff.add(cf);

        for (int j = 0; j < nr; ++j) {
            
            OneDFloatArray cJ = outputTransformed.get(j);
 
            OneDFloatArray cJPlus1 = new OneDFloatArray(
                scalingFunction.calculateNextLevel(cJ.a, j));
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            float[] s = Arrays.copyOf(cJ.a, cJ.a.length);
            for (int jj = 0; jj < s.length; ++jj) {
                s[jj] -= cJPlus1.a[jj];
            }
            OneDFloatArray wJPlus1 = new OneDFloatArray(s);
          
            outputCoeff.add(wJPlus1);
        }
    }

    /**
     * The a trous algorithm is a fast implementation of a wavelet transform 
     * with no downsampling.   It is non-orthogonal, has semi-linear runtime
     * complexity, is invariant under translation, and the transform is 
     * isotropic.
     * Implemented from pseudocode in http://www.multiresolution.com/svbook.pdf
     * The scaling function used is the higher resolution choice, the 3rd order 
     * B Spline function.
     * <pre>
     * The method uses recursive convolution operations, including previous
       * result to make next.
       * Each convolution uses two passes of one dimensional binomial kernels,
       * starting with the equivalent of sigma=1.
       * For each step, the equivalent resulting sigma is from 
       * sigma^2 = sigma_1^2 + sigma_2^2.
       * 
       * outputTransformed[1] = sigma = 1 convolution
       * outputTransformed[2] = sqrt( (1)^2 + (1)^2 ) = sqrt(2) convolution
       * outputTransformed[3] = sqrt( 2 + 1 ) = sqrt(3) convolution
       * outputTransformed[4] = sqrt( 3 + 1 ) = sqrt(4) = 2 convolution
       * outputTransformed[5] = sqrt( 4 + 1 ) = sqrt(5) convolution
       * outputTransformed[6] = sqrt( 5 + 1 ) = sqrt(6) convolution
       * ...
       * outputTransformed[8] = sqrt( 8 + 1 ) = sqrt(9) = 3 convolution
       * </pre>
     @param input
     @param outputTransformed
     @param outputCoeff 
     */
    public void calculateWithB3SplineScalingFunction(float[] input,
        List<OneDFloatArray> outputTransformed, List<OneDFloatArray> outputCoeff) {

        int imgDimen = input.length;

        int nr = (int)(Math.log(imgDimen)/Math.log(2));

        B3SplineFunction1D scalingFunction = new B3SplineFunction1D();
        
        OneDFloatArray cp = new OneDFloatArray(Arrays.copyOf(input, input.length));
        outputTransformed.add(cp);
        
        OneDFloatArray cf = new OneDFloatArray(new float[input.length]);
        outputCoeff.add(cf);

        for (int j = 0; j < nr; ++j) {
            
            OneDFloatArray cJ = outputTransformed.get(j);
 
            OneDFloatArray cJPlus1 = new OneDFloatArray(
                scalingFunction.calculate(cJ.a));
           
            outputTransformed.add(cJPlus1);
            
            // c_(j,k) − c_(j+1,k)
            float[] s = Arrays.copyOf(cJ.a, cJ.a.length);
            for (int jj = 0; jj < s.length; ++jj) {
                s[jj] -= cJPlus1.a[jj];
            }
            OneDFloatArray wJPlus1 = new OneDFloatArray(s);
            
            outputCoeff.add(wJPlus1);
        }
        
    }
    
    /**
     *
     @param c0
     @param mmCoeff
     @return
     */
    public float[] reconstruct(float[] c0, List<OneDFloatArray> mmCoeff) {

        int nr = mmCoeff.size();

        OneDFloatArray output = new OneDFloatArray(Arrays.copyOf(c0, c0.length));

        for (int j = 0; j < nr; ++j) {
           
            float[] add = mmCoeff.get(j).a;
            for (int jj = 0; jj < add.length; ++jj) {
                output.a[jj] += add[jj];
            }
        }

        return output.a;
    }
    
    /**
     * Following
     * "Edge-Optimized À-Trous Wavelets for Local Contrast Enhancement with 
     * Robust Denoising" by Johannes Hanika, Holger Dammertz, and Hendrik Lensch
     * https://jo.dreggn.org/home/2011_atrous.pdf
     * to estimate del c_i_jj as part of creating an error image.
     * The authors calculate gradient c_i_jj using Cranley Patterson rotation 
       sampling within the A Trous B3Spline window (which is 25 pixels).
       
       This looks a little like calculating auto-correlation, except not wanting 
       the center pixel as the fixed first pixel of the difference.
       
       If gradient c_i_jj is meant to be a measure of local image noise, would 
       presumably want to select only differences between adjacent pixel pairs.
       So the use of Cranley Patterson rotation must be in selecting the second
       point using an offset chosen from the vector U of values.
       That offset is applied uniformly to the set to help choose the 2nd point.
       The universe of offsets U can only be the offsets to result in the 8
       neighbor region.
        
       Not sure, but I think that is what the authors implemented.
        
       Given to this method are the center pixel index for the A Trous window
       and the offsets as dx and dy chosen from the universe U of 8 neighbor
       offsets.
        
       For each pixel in the window, will determine its intensity difference 
       from the pixel and the pixel that is it's coordinates plus the offsets.
       The result returned will be the average of those.
       
       Note that another paper 
       ("Efficient Multidimensional Sampling" by Kollig and Keller, 
       http://www.uni-kl.de/AG-Heinrich/EMS.pdf)
       suggests different sampling methods, so may change this in the future.
        
     @return 
     */
    /*
    private float estimateLocalNoise(float[] data, int pixIdx, int xOffset) {
        
        int len = data.length;
        
        float v = data[pixIdx];
        
        int count = 0;
        float diff = 0;
        // iterate within window to find first pixel
        for (int dx = -2; dx <= 2; ++dx) {
            int x1 = pixIdx + dx;
            if (x1 < 0 || x1 > (len - 1)) {
                continue;
            }
                          
            int x2 = x1 + xOffset;
            if ((x2 < 0) || (x2 > (len - 1))) {
                continue;
            }

            diff += Math.abs(data[x1] - data[x2]);

            count++;
        }
        assert(count > 0);
        
        diff /= (float)count;
        
        return diff;
    }
    */
}
