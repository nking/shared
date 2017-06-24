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
        
        float[] a = createGaussian(sigma, 0);
        
        //FWHM is approx 2.35 * sigma
        fwhm = measureFWHM(a);
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
            fwhm = measureFWHM(b.a);
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
      
        /*
        System.out.println("input=" + Arrays.toString(input));
        System.out.println("last=" + Arrays.toString(
            outputTransformed.get(outputTransformed.size() - 1).a));
          
        float[] x = new float[input.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i;
        }
        
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(0.f, input.length, 0.f, 2.f*bias, 
            x, input, x, input, 
            "input");
        plotter.addPlot(0.f, input.length, 0.f, 2.f*bias, 
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            x, 
            outputTransformed.get(outputTransformed.size() - 1).a, 
            "transformed");
        System.out.println(plotter.writeFile());
        */
        
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
        
        float[] a = createGaussian(sigma, 0);
       
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
    
    /*
     The code below is from project
     http://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)
     */
    
    private float[] createGaussian(float sigma, float mu) {
        
        float normalization = (float)(sigma * Math.sqrt(2.f * Math.PI));
        
        float hwi = estimateHWZI(sigma, 0.001f);
        if (hwi < 0) {
            hwi *= -1.f;
        }
        int halfWidthInPixels = (int)Math.ceil(hwi);
        
        int start = -1*halfWidthInPixels;
        int stopExcl = halfWidthInPixels + 1;
        
        float d, dsq;
       
        int nPoints = stopExcl - start;
        
        float[] yPoints = new float[nPoints];
        int count = 0;
        for (int i = start; i < stopExcl; i++) {
            
            float x = i;
            d = (x - mu);
            dsq = d*d;
            
            float y = (float) Math.exp(-1.f * dsq/(2.f * sigma * sigma));
            
            yPoints[count] = y/normalization;
                        
            count++;
        }
                                
        return yPoints;
    }
    
    private float estimateHWZI(float sigma, float fractionMax) {
        
        /*
           exp( (-(x0 - mu)^2)/2o~^2) / exp( (-(xcenter - mu)^2)/2o~^2) = fractionMax
        
           where xcenter = mu is the center of the gaussian
        
           exp( (-(x0 - mu)^2)/2o~^2) / 1 = fractionMax
        
           (-(x0 - mu)^2)/2o~^2) = ln(fractionMax)
        
           -(x0 - mu)^2 = -1*(2o~^2 * ln(fractionMax))
           
           x0 - mu = math.sqrt(-1*(2o~^2 * ln(fractionMax)))
        
           x0 = mu + o~ * math.sqrt(-1*(2 * ln(fractionMax)))
        
        */

        float x0 = (float)(sigma * Math.sqrt(-2. * Math.log(fractionMax)));
                
        return x0;
    }
    
    private float measureFWHM(float[] x) {
                
        int xMaxIdx = MiscMath0.findYMaxIndex(x);
        float xMax = x[xMaxIdx];
                
        int xHalfMax0Idx = -1;
        int xHalfMax1Idx = -1;
        for (int i = xMaxIdx; i > -1; i--) {
            if (x[i] > (xMax / 2.)) {
                xHalfMax0Idx = i;
            } else {
                break;
            }
        }
        for (int i = xMaxIdx - 1; i < x.length; i++) {
            if (x[i] > (xMax / 2.)) {
                xHalfMax1Idx = i;
            } else {
                break;
            }
        }
        
        //System.out.println("max=" + xMax + "\n" +
        //   " " + xHalfMax0Idx + " : " + xMaxIdx + " : " + xHalfMax1Idx);
        
        if ((xHalfMax1Idx == -1) || (xHalfMax0Idx == -1)) {
            return -1;
        }
        return xHalfMax1Idx - xHalfMax0Idx;
    }
}
