package algorithms.signalProcessing;

import algorithms.misc.MiscMath0;
import algorithms.util.OneDFloatArray;
import algorithms.util.PolygonAndPointPlotter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ATrousWaveletTransform1DTest extends TestCase {
    
    public ATrousWaveletTransform1DTest(String testName) {
        super(testName);
    }

    public void testTriangleFunction() {
        
        ATrousWaveletTransform1D wave = new ATrousWaveletTransform1D();
        float fwhm;
        
        float sigma = 1.f;
        
        float[] a = Util.createGaussian(sigma, 0);
        
        //FWHM is approx 2.35 * sigma
        fwhm = Util.measureFWHM(a);
        System.out.println("fwhm=" + fwhm);
        assertTrue(Math.abs(fwhm - 2.35*sigma) <= 1);
        
        List<OneDFloatArray> outputTransformed = new 
            ArrayList<OneDFloatArray>();
        List<OneDFloatArray> outputCoeff = new 
            ArrayList<OneDFloatArray>();
        wave.calculateWithTriangleScalingFunction(a, outputTransformed, 
            outputCoeff);
      
        for (OneDFloatArray b : outputTransformed) {
            //System.out.println("b=" + Arrays.toString(b.a));
            fwhm = Util.measureFWHM(b.a);
            System.out.println("fwhm=" + fwhm);
            
            float sigma2 = fwhm/2.35f;
            double s = Math.sqrt(sigma2*sigma2 - sigma*sigma);
            if (!Double.isNaN(s)) {
                System.out.println("implied sigma for level convolution=" + 
                    s + " resulting in " + sigma2);
            }
        }
        
        /*
        s^2 = sigma^2 + 0.707^2 --> s2-s02
        s = sqrt(sigma^2 + 0.707^2)
        */
    //    System.out.println("a=" + Arrays.toString(a));
    //    System.out.println("last=" + Arrays.toString(
    //        outputTransformed.get(outputTransformed.size() - 1).a));
    
    }

    public void testB3SplineScalingFunction() throws IOException {
    
        ATrousWaveletTransform1D wave = new ATrousWaveletTransform1D();
               
        float[] input = createCurve1();
        
        List<OneDFloatArray> outputTransformed = new 
            ArrayList<OneDFloatArray>();
        List<OneDFloatArray> outputCoeff = new 
            ArrayList<OneDFloatArray>();
        wave.calculateWithB3SplineScalingFunction(input, outputTransformed, 
            outputCoeff);
      
        
        //System.out.println("input=" + Arrays.toString(input));
        //System.out.println("last=" + Arrays.toString(
        //    outputTransformed.get(outputTransformed.size() - 1).a));
          
        float[] x = new float[input.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i;
        }
        
        // errors:
        float[] diffSq = Arrays.copyOf(outputCoeff.get(outputTransformed.size() - 1).a,
            outputCoeff.get(outputTransformed.size() - 1).a.length);
        float[] diffX = new float[diffSq.length];
        for (int i = 0; i < diffSq.length; ++i) {
            //TODO: should add these differences in quadrature for all 
            //   transformations up to the index of transformed, 
            //   which is the last index here
            diffSq[i] *= diffSq[i];
            diffSq[i] *= outputTransformed.get(outputTransformed.size() - 1).a[i];
        }
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(0.f, input.length, 0.f, 10, 
            x, input, x, input, 
            "input");
        plotter.addPlot(0.f, input.length, 0.f, 10, 
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            diffX, diffSq,
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            "transformed");
        plotter.addPlot(0.f, input.length, 
            0.9f*MiscMath0.findMin(outputCoeff.get(outputTransformed.size() - 1).a), 
            1.2f*MiscMath0.findMax(outputCoeff.get(outputTransformed.size() - 1).a), 
            x, 
            outputCoeff.get(outputTransformed.size() - 1).a, 
            x, 
            outputCoeff.get(outputTransformed.size() - 1).a, 
            "coeff");
        
        System.out.println(plotter.writeFile());
        
        
        int min0 = MiscMath0.findYMaxIndex(input);
        int min1 = MiscMath0.findYMaxIndex(
            outputTransformed.get(outputTransformed.size() - 1).a);
        assertEquals(min0, min1);
    }

    public void testReconstruct() throws IOException {
        
        ATrousWaveletTransform1D wave = new ATrousWaveletTransform1D();
        
        float[] input = createCurve1();
        
        List<OneDFloatArray> outputTransformed = new 
            ArrayList<OneDFloatArray>();
        List<OneDFloatArray> outputCoeff = new 
            ArrayList<OneDFloatArray>();
        wave.calculateWithB3SplineScalingFunction(input, outputTransformed, 
            outputCoeff);
        
        float[] r = wave.reconstruct(
            outputTransformed.get(outputTransformed.size() - 1).a, 
            outputCoeff);
    
        float[] x = new float[input.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i;
        }
        
        /*
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(0.f, input.length, 0.f, 10.f, 
            x, input, x, input, 
            "input");
        plotter.addPlot(0.f, input.length, 0.f, 10.f, 
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            "transformed");
        plotter.addPlot(0.f, input.length, 0.f, 10.f, 
            x, r, 
            x, r, 
            "reconstructed");
        
        System.out.println(plotter.writeFile2());
        */
            
        for (int i = 0; i < input.length; ++i) {
            float v0 = input[i];
            float v1 = r[i];
            float diff = Math.abs(v0 - v1);
            //System.out.println("diff=" + diff);
            assertEquals(0.f, diff);
        }
        
    }
    
    private float[] createCurve1() {
        float sigma = 1.f;
        
        float[] a = Util.createGaussian(sigma, 0);
       
        int n = a.length;
        float bias = 5.f;
        
        float[] input = new float[n*3];
        for (int i = 0; i < n; ++i) {
            input[i] = bias + 4.f*a[i];
        }
        int idx = n;
        for (int i = 0; i < n; ++i) {
            input[idx] = bias - 4.f*a[i];
            idx++;
        }
        for (int i = 0; i < n; ++i) {
            input[idx] = bias + 4.f*a[i];
            idx++;
        }
        return input;
    }
    
}
