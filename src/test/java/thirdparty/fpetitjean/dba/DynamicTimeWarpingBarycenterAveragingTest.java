package thirdparty.fpetitjean.dba;

import algorithms.util.FormatArray;
import java.util.Arrays;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 */
public class DynamicTimeWarpingBarycenterAveragingTest extends TestCase {
    
    public DynamicTimeWarpingBarycenterAveragingTest(String testName) {
        super(testName);
    }
    
    /*
    characters to numbers: can use ascii code, can use java's string caracter codes
    jaccard similarity, hamming distance, 
    differences in images of the characters as font normalized correlation coefficients
    of summed pixel differences.
    
    related topic, but for determining classes is building a data dictionary:
    "Time Series Classification under More Realistic Assumptions" by Hu, chen, & Keogh
    useful for the UCR datasets
    www.cs.ucr.edu/~eamonn/time_series_data/
    */
    
    public void est0() {
        int nSeries = 20;
        int length = 10;
        Random r = new Random(3071980);

        double[][] sequences = new double[nSeries][];
        for (int i = 0; i < sequences.length; i++) {
            sequences[i] = new double[length];
            for (int j = 0; j < sequences[i].length; j++) {
                sequences[i][j] = Math.cos(r.nextDouble() * j / 20.0 * Math.PI);
            }
        }
        System.out.printf("sequence[=%s\n", FormatArray.toString(sequences, "%.3f"));

        double[] averageSequence = DynamicTimeWarpingBarycenterAveraging.performDBA(sequences, 10);

        System.out.print("average sequenc=[");
        for (int j = 0; j < averageSequence.length; j++) {
            System.out.print(averageSequence[j] + " ");
        }
        System.out.println("]");
    }
    
    public void test1() {
        // from fig 2 of Petitjea & Gancarski, 2011
        // "Summarizing a set of time series by averaging: From Steiner sequence to
        // compact multiple alignment"
        
        int nIter = 1;
        
        double[][] sequences = new double[3][];
        sequences[0] = new double[]{1,10,0,0,4};
        sequences[1] = new double[]{0,2,10,0,0};
        sequences[2] = new double[]{0, 0, 0, 10, 0};
        
        double[] averageSequence = DynamicTimeWarpingBarycenterAveraging.performDBA(sequences, nIter);
        System.out.print("average sequenc=[");
        for (int j = 0; j < averageSequence.length; j++) {
            System.out.print(averageSequence[j] + " ");
        }
        System.out.println("]");
    }
    
    /*
    converts a to lower case ascii and returns as integers.
    for more complex conversions, could use the unicode multiple character codes,
    or jaccard similarity.
    */
    public static int[] convertToAscii(String a) {
        a = a.toLowerCase();
        int[] out = new int[a.length()];
        for (int i = 0; i < a.length(); ++i) {
            out[i] = (int)a.charAt(i);
        }
        return out;
    }
    public static double[] convertToAsciiReal(String a) {
        int[] aAscii = convertToAscii(a);
        double[] out = new double[a.length()];
        for (int i = 0; i < a.length(); ++i) {
            out[i] = aAscii[i];
        }
        return out;
    }
}
