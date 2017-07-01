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
public class MedianTransform1DTest extends TestCase {
    
    public MedianTransform1DTest(String testName) {
        super(testName);
    }

    public void test() throws IOException {
    
        MedianTransform1D wave = new MedianTransform1D();
               
        float[] input = createCurve1();
        
        List<OneDFloatArray> outputTransformed = new 
            ArrayList<OneDFloatArray>();
        List<OneDFloatArray> outputCoeff = new 
            ArrayList<OneDFloatArray>();
        
        wave.multiscalePyramidalMedianTransform2(
            input, outputTransformed, outputCoeff);
      
        /*
        System.out.println("input=" + Arrays.toString(input));
        
        float[] x = new float[input.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i;
        }
         
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter();
        plotter.addPlot(0.f, input.length, 0.f, 2.f*5, 
            x, input, x, input, 
            "input");
        
        for (int i = 0; i < outputTransformed.size(); ++i) {
            
            //System.out.println("last=" + Arrays.toString(
            //    outputTransformed.get(i).a));
            
            x = new float[outputTransformed.get(i).a.length];
            for (int ii = 0; ii < x.length; ++ii) {
                x[ii] = ii;
            }
         
            plotter.addPlot(0.f, outputTransformed.get(i).a.length, 
                0.f, 2.f*5, 
            x, outputTransformed.get(i).a, 
            x, outputTransformed.get(i).a, 
            "transformed");
        }
               
        System.out.println(plotter.writeFile());
        */
        
        /*
        plotter = new PolygonAndPointPlotter();
        
        for (int i = 0; i < outputCoeff.size(); ++i) {
            
            //System.out.println("last=" + Arrays.toString(
            //    outputTransformed.get(i).a));
            
            x = new float[outputCoeff.get(i).a.length];
            for (int ii = 0; ii < x.length; ++ii) {
                x[ii] = ii;
            }
            
            float yMin = MiscMath0.findMin(outputCoeff.get(i).a);
            float yMax = MiscMath0.findMax(outputCoeff.get(i).a);
         
            plotter.addPlot(0.f, outputCoeff.get(i).a.length, 
                yMin, yMax, 
            x, outputCoeff.get(i).a, 
            x, outputCoeff.get(i).a, 
            "coeff");
        }
               
        System.out.println(plotter.writeFile2());
        */
        
        
        //plotter = new PolygonAndPointPlotter();
        
        float[] r = wave.reconstructPyramidalMultiscaleMedianTransform(
            outputTransformed.get(outputTransformed.size() - 1),
            outputCoeff);
        
        /*
        x = new float[input.length];
        for (int i = 0; i < x.length; ++i) {
            x[i] = i;
        }
        
        plotter.addPlot(0.f, input.length, 0.f, 10.f, 
            x, input, x, input, 
            "input");
        
        plotter.addPlot(0.f, input.length, 0.f, 10.f, 
            x, r, 
            x, r, 
            "reconstructed");
        
        System.out.println(plotter.writeFile3());
        */
            
        assertEquals(r.length, input.length);
        
        for (int i = 0; i < input.length; ++i) {
            float v0 = input[i];
            float v1 = r[i];
            float diff = Math.abs(v0 - v1);
            //System.out.println("diff=" + diff);
            assertTrue(Math.abs(diff) < 1);
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
