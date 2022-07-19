package algorithms.statistics;

import java.util.logging.Logger;
import junit.framework.TestCase;

public class GeneralizedExtremeValueTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getName());
    
    public void testGenerateCurve() throws Exception {
        
        double[] xPoints;
        double[] yPoints;
        double[] dXPoints;
        double[] dYPoints;
        GeneralizedExtremeValue gev;
        double[] curve;
        double sigma, k, mu;
        yPoints = new double[0];
        dXPoints = new double[0];
        dYPoints = new double[0];
        
        // TypeII
        xPoints = new double[11];
        for (int i = 0; i < 11; i++) {
            xPoints[i] = -4+i;
        }
        k = 0.5f;
        sigma = 1.0f;
        mu = 0.0f;
        gev = new GeneralizedExtremeValue(xPoints, yPoints, dXPoints, dYPoints);
        curve = gev.generateCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.15f) < 0.01f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.05f);
        assertTrue(Math.abs(curve[5] - 0.19f) < 0.05f);
        assertTrue(Math.abs(curve[6] - 0.097f) < 0.01f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        assertTrue(Math.abs(curve[8] - 0.03f) < 0.01f);
        
        curve = GeneralizedExtremeValue.genCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.15f) < 0.01f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.05f);
        assertTrue(Math.abs(curve[5] - 0.19f) < 0.05f);
        assertTrue(Math.abs(curve[6] - 0.097f) < 0.01f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        assertTrue(Math.abs(curve[8] - 0.03f) < 0.01f);
        
        for (int i = 0; i < xPoints.length; i++) {
            Double a = GeneralizedExtremeValue.generateYGEV(xPoints[i], mu, sigma, k);
            switch(i) {
                case 3:
                    assertTrue(Math.abs(a - 0.15f) < 0.01f);
                    break;
                case 4:
                    assertTrue(Math.abs(a - 0.37f) < 0.05f);
                    break;
                case 5:
                    assertTrue(Math.abs(a - 0.19f) < 0.05f);
                    break;
                case 6:
                    assertTrue(Math.abs(a - 0.097f) < 0.01f);
                    break;
                case 7:
                    assertTrue(Math.abs(a - 0.05f) < 0.01f);
                    break;
                case 8:
                    assertTrue(Math.abs(a - 0.03f) < 0.01f);
                    break;
                default:
                    break;
            }
        }
        
        // TypeI
        k = 0.0f;
        sigma = 1.0f;
        mu = 0.0f;
        gev = new GeneralizedExtremeValue(xPoints, yPoints, dXPoints, dYPoints);
        curve = gev.generateCurve(xPoints, mu, sigma, k);
        assertNotNull(curve);        
        assertTrue(Math.abs(curve[3] - 0.18f) < 0.1f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.01f);
        assertTrue(Math.abs(curve[5] - 0.25f) < 0.01f);
        assertTrue(Math.abs(curve[6] - 0.12f) < 0.05f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        
        curve = GeneralizedExtremeValue.genCurve(xPoints, mu, sigma, k);
        assertTrue(Math.abs(curve[3] - 0.18f) < 0.1f);
        assertTrue(Math.abs(curve[4] - 0.37f) < 0.01f);
        assertTrue(Math.abs(curve[5] - 0.25f) < 0.01f);
        assertTrue(Math.abs(curve[6] - 0.12f) < 0.05f);
        assertTrue(Math.abs(curve[7] - 0.05f) < 0.01f);
        
        for (int i = 0; i < xPoints.length; i++) {
            Double a = GeneralizedExtremeValue.generateYEVTypeI(xPoints[i], mu, sigma);
            switch(i) {
                case 3:
                    assertTrue(Math.abs(a - 0.18f) < 0.1f);
                    break;
                case 4:
                    assertTrue(Math.abs(a - 0.37f) < 0.01f);
                    break;
                case 5:
                    assertTrue(Math.abs(a - 0.25f) < 0.01f);
                    break;
                case 6:
                    assertTrue(Math.abs(a - 0.12f) < 0.05f);
                    break;
                case 7:
                    assertTrue(Math.abs(a - 0.05f) < 0.01f);
                    break;
                default:
                    break;
            }
        }
    }
    
}
