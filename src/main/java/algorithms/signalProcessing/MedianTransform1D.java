package algorithms.signalProcessing;

import algorithms.util.OneDFloatArray;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author nichole
 */
public class MedianTransform1D {
    
    /**
     * pyramidal median transform (faster than multiscalePyramidalMedianTransform
     * but reconstruction from coefficients is not exact;
     * following pseudocode in http://www.multiresolution.com/svbook.pdf
     * "Handbook of Astronomical Data Analysis" by
     * Jean-Luc Starck and Fionn Murtagh, pg 121.
     * 
     * Note that the length of outputCoeff is one less than 
     * outputTransformed because the iteration stops and cannot calculate the
     * next difference.
     @param input
     @param outputTransformed
     @param outputCoeff 
     */
    public void multiscalePyramidalMedianTransform2(
        float[] input,
        List<OneDFloatArray> outputTransformed, 
        List<OneDFloatArray> outputCoeff) {

        int imgDimen = input.length;

        int nr = (int)(Math.log(imgDimen)/Math.log(2));
        int s = 1;
        int winL = 2*s + 1;
        
        Interp interp = new Interp();
        
        MedianSmooth1D med = new MedianSmooth1D();
        
        outputTransformed.add(
            new OneDFloatArray(Arrays.copyOf(input, input.length)));
        
        outputCoeff.add(
            new OneDFloatArray(new float[input.length]));

        for (int j = 0; j < (nr - 1); ++j) {
                       
            OneDFloatArray cJ = outputTransformed.get(j);
            
            if (cJ.a.length <= winL) {
                break;
            }
            
            float[] cJPlus1Ast = med.calculate(cJ.a, winL);   
            
            assert(cJ.a.length == cJPlus1Ast.length);
                        
            // decimation:
            float[] cJPlus1;
            if ((cJPlus1Ast.length & 1) == 1) {
                int outLength = cJPlus1Ast.length/2;
                cJPlus1 = interp.linearInterp(
                    cJPlus1Ast, outLength);//, -256, 255);
            } else {
                cJPlus1 = interp.bin(cJPlus1Ast, 2);
            }
            outputTransformed.add(new OneDFloatArray(cJPlus1));
            
            OneDFloatArray wJPlus1 = new OneDFloatArray(new float[cJ.a.length]);
            for (int ii = 0; ii < cJPlus1Ast.length; ++ii) {
                wJPlus1.a[ii] = cJ.a[ii] - cJPlus1Ast[ii];
            }            
            
            outputCoeff.add(wJPlus1);
            
            assert(cJ.a.length == wJPlus1.a.length);
        }
        outputCoeff.remove(0);
    }
    
     /**
     * reconstruct image from products of pyramidal median transform.
     * following pseudocode in http://www.multiresolution.com/svbook.pdf
     * "Handbook of Astronomical Data Analysis" by
     * Jean-Luc Starck and Fionn Murtagh
     * 
     @param c0
     @param mmCoeff
     @return 
     */
    public float[] reconstructPyramidalMultiscaleMedianTransform(
        OneDFloatArray c0, List<OneDFloatArray> mmCoeff) {

        int nr = mmCoeff.size();

        Interp interp = new Interp();
        
        float[] output = Arrays.copyOf(c0.a, c0.a.length);

        for (int j = (nr - 1); j > -1; --j) {

            OneDFloatArray wJ = mmCoeff.get(j);
            
            //up-sample wJ to size output
            float[] cJPrime;
            if (output.length * 2 == wJ.a.length) {
                
                cJPrime = interp.unbin(output, 2);
               
            } else {
                
                cJPrime = interp.linearInterp(
                    output, wJ.a.length);//, -256, 255);
            }
            
            output = Arrays.copyOf(cJPrime, cJPrime.length);
            for (int ii = 0; ii < wJ.a.length; ++ii) {
                output[ii] += wJ.a[ii];
            }
            
        }

        return output;
    }

}
