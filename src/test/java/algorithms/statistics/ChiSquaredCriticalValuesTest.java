package algorithms.statistics;

import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ChiSquaredCriticalValuesTest extends TestCase {
    
    public ChiSquaredCriticalValuesTest(String testName) {
        super(testName);
    }
    
    public void testLin() {
        
        double p, diff, chiSquared2;
        int dOF;
        
        double[] x2 = new double[]{2.706, 3.841, 5.024, 6.635, 10.828};
        double[] ep = new double[]{0.90, 0.95,  0.975,  0.99,  0.999};
        dOF = 1;
                
        int c=0;
        for (double chiSquared : x2) {  
            p = ChiSquaredCriticalValues.approxPValueLin(chiSquared, dOF);
            System.out.println("  expected 1-p=" + ep[c] + " => " + p);
            diff = Math.abs(p - ep[c]);
            assertTrue(diff < 0.08);
            
            chiSquared2 = ChiSquaredCriticalValues.approxChiSqStatLin(1.-p, dOF);
            System.out.println("  expected x^2=" + chiSquared + " => " + chiSquared2);
            diff = Math.abs(chiSquared - chiSquared2);
            assertTrue(diff < 0.002);
            
            c++;
        }
        
        x2 = new double[]{9.236, 11.070, 12.833, 15.086, 20.515};
        dOF = 5;
        c=0;
        for (double chiSquared : x2) {  
            p = ChiSquaredCriticalValues.approxPValueLin(chiSquared, dOF);
            System.out.println("  expected 1-p=" + ep[c] + " => " + p);
            diff = Math.abs(p - ep[c]);
            assertTrue(diff < 0.02);
            
            chiSquared2 = ChiSquaredCriticalValues.approxChiSqStatLin(1.-p, dOF);
            System.out.println("  expected x^2=" + chiSquared + " => " + chiSquared2);
            diff = Math.abs(chiSquared - chiSquared2);
            assertTrue(diff < 0.002);
            
            c++;
        }
        
        x2 = new double[]{118.498, 124.342, 129.561, 135.807, 149.449};
        dOF = 100;
        c=0;
        for (double chiSquared : x2) {  
            p = ChiSquaredCriticalValues.approxPValueLin(chiSquared, dOF);
            System.out.println("  expected 1-p=" + ep[c] + " => " + p);
            diff = Math.abs(p - ep[c]);
            assertTrue(diff < 0.01);
            
            chiSquared2 = ChiSquaredCriticalValues.approxChiSqStatLin(1.-p, dOF);
            System.out.println("  expected x^2=" + chiSquared + " => " + chiSquared2);
            diff = Math.abs(chiSquared - chiSquared2);
            assertTrue(diff < 0.002);
            
            c++;
        }
        
        
        // lower tail:
        x2 = new double[]{1.610, 1.145, .831, .554, .210};
        ep = new double[]{0.10, 0.05, 0.025, 0.01, 0.001};
        dOF = 5;
        c=0;
        for (double chiSquared : x2) {  
            p = ChiSquaredCriticalValues.approxPValueLin(chiSquared, dOF);
            System.out.println("  expected p=" + ep[c] + " => " + p);
            diff = Math.abs(p - ep[c]);
            assertTrue(diff < 0.02);
            
            chiSquared2 = ChiSquaredCriticalValues.approxChiSqStatLin(1.-p, dOF);
            System.out.println("  expected x^2=" + chiSquared + " => " + chiSquared2);
            diff = Math.abs(chiSquared - chiSquared2);
            assertTrue(diff < 0.002);
            
            c++;
        }
        
        // lower tail:
        x2 = new double[]{82.358, 77.929, 74.222, 70.065, 61.91};
        ep = new double[]{0.10, 0.05, 0.025, 0.01, 0.001};
        dOF = 100;
        c=0;
        for (double chiSquared : x2) {  
            p = ChiSquaredCriticalValues.approxPValueLin(chiSquared, dOF);
            System.out.println("  expected p=" + ep[c] + " =>" + p);
            diff = Math.abs(p - ep[c]);
            assertTrue(diff < 0.01);
            
            chiSquared2 = ChiSquaredCriticalValues.approxChiSqStatLin(1.-p, dOF);
            System.out.println("  expected x^2=" + chiSquared + " => " + chiSquared2);
            diff = Math.abs(chiSquared - chiSquared2);
            assertTrue(diff < 0.002);
            
            c++;
        }
    }
    
    /*public void test0() {
                
        double p;
        int dOF;
        
        double[] x2 = new double[]{2.706, 3.841, 5.024, 6.635, 10.828};
        dOF = 5;
                
        for (double chiSquared : x2) {
            
            p = ChiSquaredCriticalValues.approxUpperTailPValueBeh(chiSquared, dOF);
            
        }
    }*/
}
