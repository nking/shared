package algorithms.statistics;

import algorithms.matrix.MatrixUtil;
import algorithms.misc.MiscMath0;
import java.security.SecureRandom;
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
        double[] st, expectedSt, dSt;
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
        
        dSt = Standardization.standardUnitDenormalization(st, nDimensions, 
            outputMean, outputStDev);
        
        assertEquals(data.length, dSt.length);
        nc = 0;
        for (int i = 0; i < data.length; ++i) {
            if ((Math.abs(dSt[i] - data[i]) < tol)) {
                nc++;
            }
        }
        assertEquals(data.length, nc);
        
        
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
        
        dSt = Standardization.standardUnitDenormalization(st, nDimensions, 
            outputMean, outputStDev);
        
        assertEquals(data.length, dSt.length);
        nc = 0;
        for (int i = 0; i < data.length; ++i) {
            if ((Math.abs(dSt[i] - data[i]) < tol)) {
                nc++;
            }
        }
        assertEquals(data.length, nc);
        
        // ------------
        //(-20, 48), (-20, -48), (20, 0), (59, 0)
        nDimensions = 2;
        data = new double[]{-20, 48, -20, -48, 20, 0, 59, 0};
        expectedMean = new double[]{9.75, 0};
        expectdStDev = new double[]{37.863, 39.192};
        expectedSt = new double[]{-0.7857, 1.225, -0.786, -1.225, 0.271, 0, 1.301, 0};
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
        
        dSt = Standardization.standardUnitDenormalization(st, nDimensions, 
            outputMean, outputStDev);
        
        assertEquals(data.length, dSt.length);
        nc = 0;
        for (int i = 0; i < data.length; ++i) {
            if ((Math.abs(dSt[i] - data[i]) < tol)) {
                nc++;
            }
        }
        assertEquals(data.length, nc);
    }

    public void test2() {
        SecureRandom sr = new SecureRandom();
        long seed = System.nanoTime();
        System.out.println("test2: seed=" + seed);
        sr.setSeed(seed);
        
        int nIter = 10;
        double tol = 0.001;
        
        int nData, nDimensions;
        double[] data, st, outputM, outputStDev, dSt;
        double[] mean2;
        double[][] stDev2;
        int i, it, nc0, nc1, nc;
        
        for (it = 0; it < nIter; ++it) {
            nData = sr.nextInt(1023) + 5;
            nDimensions = sr.nextInt(nData / 4) + 1;

            data = new double[nData * nDimensions];
            outputM = new double[nDimensions];
            outputStDev = new double[nDimensions];
            
            for (i = 0; i < data.length; ++i) {
                data[i] = sr.nextDouble();
                if (sr.nextBoolean()) {
                    data[i] *= -1.;
                }
                data[i] *= sr.nextInt();
            }
            
            st = Standardization.standardUnitNormalization(data, nDimensions, 
                outputM, outputStDev);
            
            assertEquals(st.length, data.length);
            
            mean2 = MiscMath0.mean(st, nDimensions);
            stDev2 = MiscMath0.standardDeviation(st, nDimensions);
            assertEquals(nDimensions, mean2.length);
            assertEquals(2, stDev2.length);
            assertEquals(nDimensions, stDev2[0].length);
            assertEquals(nDimensions, stDev2[1].length);
            nc0 = 0;
            nc1 = 0;
            nc = 0;
            for (i = 0; i < nDimensions; ++i) {
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

            dSt = Standardization.standardUnitDenormalization(st, nDimensions,
                    outputM, outputStDev);

            assertEquals(data.length, dSt.length);
            nc = 0;
            for (i = 0; i < data.length; ++i) {
                if ((Math.abs(dSt[i] - data[i]) < tol)) {
                    nc++;
                }
            }
            assertEquals(data.length, nc);
        }
    }
    
    public void test3() {
        SecureRandom sr = new SecureRandom();
        long seed = System.nanoTime();
        System.out.println("test3: seed=" + seed);
        sr.setSeed(seed);
        
        int nIter = 10;
        double tol = 0.001;
        
        int nData, nDimensions;
        double[] outputM, outputStDev, dSt;
        double[] mean2, stDev2;
        double[][] data, st;
        int i, j, it, nc0, nc1, nc;
        
        for (it = 0; it < nIter; ++it) {
            nData = sr.nextInt(1023) + 5;
            nDimensions = sr.nextInt(nData / 4) + 1;

            data = new double[nData][nDimensions];
            outputM = new double[nDimensions];
            outputStDev = new double[nDimensions];
            
            for (i = 0; i < data.length; ++i) {
                data[i] = new double[nDimensions];
                for (j = 0; j < data[i].length; ++j) {
                    data[i][j] = sr.nextDouble();
                    if (sr.nextBoolean()) {
                        data[i][j] *= -1.;
                    }
                    data[i][j] *= sr.nextInt();
                }
            }
            
            st = Standardization.standardUnitNormalization(data, 
                outputM, outputStDev);
            
            assertEquals(st.length, data.length);
            
            mean2 = MatrixUtil.columnMeans(st);
            stDev2 = MatrixUtil.standardDeviation(st);
            assertEquals(nDimensions, mean2.length);
            assertEquals(nDimensions, stDev2.length);
            nc0 = 0;
            nc = 0;
            for (i = 0; i < nDimensions; ++i) {
                if (Math.abs(mean2[i] - 0.) < tol) {
                    nc0++;
                }
                if ((Math.abs(stDev2[i] - 1.) < tol) || (Math.abs(stDev2[i] - 0.) < tol)) {
                    nc++;
                }
            }
            assertEquals(nDimensions, nc0);
            assertEquals(nDimensions, nc);

            /*
            dSt = Standardization.standardUnitDenormalization(st, nDimensions,
                    outputM, outputStDev);

            assertEquals(data.length, dSt.length);
            nc = 0;
            for (i = 0; i < data.length; ++i) {
                if ((Math.abs(dSt[i] - data[i]) < tol)) {
                    nc++;
                }
            }
            assertEquals(data.length, nc);
            */
        }
    }

    public void testUnitStandardNormalizeImage() {

        double[][] im = new double[][]{{3,9,12}};
        double[] meanStd = Standardization.unitStandardNormalizeImage(im);
        assertTrue(Math.abs(meanStd[0] - 8) < 1E-7);
        assertTrue(Math.abs(meanStd[1] - 4.58257569) < 1E-7);

        double[] res = MiscMath0.getAvgAndStDev(im[0]);
        assertTrue(Math.abs(res[0]) < 1E-7);
        assertTrue(Math.abs(res[1] - 1.) < 1E-7);

        double[][] im2 = new double[][]{{3,9,12}};
        Standardization.applyMeanAndStdToImage(im2, meanStd[0], meanStd[1]);

        for (int i = 0; i < im[0].length; ++i) {
            assertTrue(Math.abs(im[0][i] - im2[0][i]) < 1E-7);
        }

    }

    public void testMinMaxNormalizeImage() {

        double[][] im = new double[][]{{3,9,12}};
        double[] minMax = Standardization.minMaxNormalizeImage(im);
        assertTrue(Math.abs(minMax[0] - 3) < 1E-7);
        assertTrue(Math.abs(minMax[1] - 12) < 1E-7);

        for (int i = 0; i < im.length; ++i) {
            for (int j = 0; j < im[i].length; ++j) {
                assertTrue(im[i][j] >= 0.);
                assertTrue(im[i][j] <= 1.);
            }
        }

        double[][] im2 = new double[][]{{3,9,12}};
        Standardization.minMaxNormalizeImage(im2, minMax[0], minMax[1]);

        for (int i = 0; i < im[0].length; ++i) {
            assertTrue(Math.abs(im[0][i] - im2[0][i]) < 1E-7);
        }

    }
    
}
