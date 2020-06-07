package algorithms.misc;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class StandardizationTest extends TestCase {
    
    public StandardizationTest(String testName) {
        super(testName);
    }

    public void testStandardUnitNormalization() {
        
        double[] data;
        int nDimensions;
        double[] outputMean, expectedMean;
        double[] outputStDev, expectdStDev;
        double[] mean2;
        double[][] stDev2;
        double[] st, expectedSt;
        double tol = 0.001;
        int nc0, nc1, nc;
        
        //(1, 1), (3, 3)                           --> (2,2)        
        nDimensions = 2;
        data = new double[]{1, 1,  3, 3}; // -1, 1
        expectedMean = new double[]{2, 2};
        expectdStDev = new double[]{1.414, 1.414};
        expectedSt = new double[]{-0.707, -0.707, 0.707, 0.707};
        outputMean = new double[nDimensions];
        outputStDev = new double[nDimensions];
        
        st = Standardization.standardUnitNormalization(data, 
            nDimensions, outputMean, outputStDev);
        
        assertEquals(st.length, expectedSt.length);
        
        nc0 = 0;  nc1=0;
        for (int i = 0; i < outputMean.length; ++i) {
            if (Math.abs(outputMean[i] - expectedMean[i]) < tol) {
                nc0++;
            }
            if (Math.abs(outputStDev[i] - expectdStDev[i]) < tol) {
                nc1++;
            }
        }
        assertEquals(expectedMean.length, nc0);
        assertEquals(expectdStDev.length, nc1);
        nc = 0;
        for (int i = 0; i < expectedSt.length; ++i) {
            if (Math.abs(st[i] - expectedSt[i]) < tol) {
                nc++;
            }
        }
        assertEquals(st.length, nc);
        
        mean2 = MiscMath0.mean(st, nDimensions);
        stDev2 = MiscMath0.standardDeviation(st, nDimensions);
        assertEquals(nDimensions, mean2.length);
        assertEquals(2, stDev2.length);
        assertEquals(nDimensions, stDev2[0].length);
        assertEquals(nDimensions, stDev2[1].length);
        nc0 = 0;  nc1=0; nc = 0;
        for (int i = 0; i < nDimensions; ++i) {
            if (Math.abs(mean2[i] - 0.) < tol) {
                nc0++;
            }
            if (Math.abs(stDev2[0][i] - 0) < tol) {
                nc1++;
            }
            if ((Math.abs(stDev2[1][i] - 1.) < tol) || (Math.abs(stDev2[1][i] - 0.) < tol)) {
                nc++;
            }
        }
        assertEquals(nDimensions, nc0);
        assertEquals(nDimensions, nc1);
        assertEquals(nDimensions, nc);
        
        
        //------------------------------------------
        //(0, 0), (0, 0), (0, 12) 
        nDimensions = 2;
        data = new double[]{0, 0, 0, 0, 0, 12};
        expectedMean = new double[]{0, 4};
        expectdStDev = new double[]{0, 6.928};
        
        expectedSt = new double[]{0, -0.5774, 0, -0.577, 0, 1.1547};
        outputMean = new double[nDimensions];
        outputStDev = new double[nDimensions];
        
        st = Standardization.standardUnitNormalization(data, 
            nDimensions, outputMean, outputStDev);
        
        assertEquals(st.length, expectedSt.length);
        
        nc0 = 0;  nc1=0;
        for (int i = 0; i < outputMean.length; ++i) {
            if (Math.abs(outputMean[i] - expectedMean[i]) < tol) {
                nc0++;
            }
            if (Math.abs(outputStDev[i] - expectdStDev[i]) < tol) {
                nc1++;
            }
        }
        assertEquals(expectedMean.length, nc0);
        assertEquals(expectdStDev.length, nc1);
        nc = 0;
        for (int i = 0; i < expectedSt.length; ++i) {
            if (Math.abs(st[i] - expectedSt[i]) < tol) {
                nc++;
            }
        }
        assertEquals(st.length, nc);
        
        mean2 = MiscMath0.mean(st, nDimensions);
        stDev2 = MiscMath0.standardDeviation(st, nDimensions);
        assertEquals(nDimensions, mean2.length);
        assertEquals(2, stDev2.length);
        assertEquals(nDimensions, stDev2[0].length);
        assertEquals(nDimensions, stDev2[1].length);
        nc0 = 0;  nc1=0; nc = 0;
        for (int i = 0; i < nDimensions; ++i) {
            if (Math.abs(mean2[i] - 0.) < tol) {
                nc0++;
            }
            if (Math.abs(stDev2[0][i] - 0) < tol) {
                nc1++;
            }
            if ((Math.abs(stDev2[1][i] - 1.) < tol) || (Math.abs(stDev2[1][i] - 0.) < tol)) {
                nc++;
            }
        }
        assertEquals(nDimensions, nc0);
        assertEquals(nDimensions, nc1);
        assertEquals(nDimensions, nc);
        
        // ------------
        //(-20, 48), (-20, -48), (20, 0), (59, 0)
        nDimensions = 2;
        data = new double[]{-20, 48, -20, -48, 20, 0, 59, 0};
        
    }

    public void testStandardUnitDenormalization() {
    }
    
}
