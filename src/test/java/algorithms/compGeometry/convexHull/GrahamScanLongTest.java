package algorithms.compGeometry.convexHull;

import algorithms.compGeometry.convexHull.GrahamScan.CHL;
import algorithms.compGeometry.MiscellaneousCurveHelper;
import algorithms.util.PairIntArray;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Random;
import java.util.logging.Logger;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.Test;
import junit.framework.TestSuite;
import junit.framework.TestCase;

public class GrahamScanLongTest extends TestCase {

    protected static Logger log = Logger.getLogger(GrahamScanLongTest.class.getName());
    
    public void estScan() throws Exception {

        /*            2,6*
         *
         *
         *     0,2*   2,2 3,2
         *                      7,1*
         *            2,0*
         *         
         */
        int n = 6;
        long[] x = new long[]{0, 2, 7, 2, 2, 3};
        long[] y = new long[]{2, 2, 1, 6, 0, 2};

        boolean useLinearSort = false;
        CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
        
        //System.out.printf("ch=%s\n", ch.toString());
        
        assertEquals(5, ch.getXH().length);

        // clockwise order
        double[] expectedxx = new double[]{0, 2, 7, 2, 0};
        double[] expectedyy = new double[]{2, 6, 1, 0, 2};

        double tol = 1E-7;
        for (int i = 0; i < expectedxx.length; i++) {
            assertTrue(Math.abs(expectedxx[i] - ch.getXH()[i]) < tol);
            assertTrue(Math.abs(expectedyy[i] - ch.getYH()[i]) < tol);
        }

        useLinearSort = true;
        x = new long[]{0, 2, 7, 2, 2, 3};
        y = new long[]{2, 2, 1, 6, 0, 2};

        ch = GrahamScan.computeHull(x, y, useLinearSort);

        assertEquals(5, ch.getXH().length);

        for (int i = 0; i < expectedxx.length; i++) {
            assertTrue(Math.abs(expectedxx[i] - ch.getXH()[i]) < tol);
            assertTrue(Math.abs(expectedyy[i] - ch.getYH()[i]) < tol);
        }
    }
    
    public void estScanExceptions() throws Exception {

        boolean threwException = false;
        boolean useLinearSort = false;
        long[] x = null;
        long[] y = null;
        try {
            CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        x = new long[]{0, 2};
        y = new long[]{2, 2};
        
        try {
            CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
        } catch (IllegalArgumentException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
        x = new long[]{0, 2, 2};
        y = new long[]{2, 2, 2};
        
        try {
            CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
        } catch (GrahamScanTooFewPointsException e) {
            threwException = true;
        }
        assertTrue(threwException);
        
    }

    public void testCalculateConvexHull6() throws Exception {

        int ntries = 10;

        long xMin = 10;
        long yMin = 10;
        long xMax = 1000;
        long yMax = 1000;
            
        PolygonAndPointPlotter plotter = new PolygonAndPointPlotter(xMin, xMax, yMin, yMax);
        
        for (int i = 0; i < ntries; i++) {
            long seed = System.nanoTime();
            //seed = 171920883721844L;
            System.out.println("seed=" + seed);
            Random sr = new Random(seed);

            int n = 1000;
            long[] x = new long[n];
            long[] y = new long[n];
            
            float[] xF = new float[n];
            float[] yF = new float[n];

            float maxRadius = 200;

            for (int ii = 0; ii < n; ii++) {

                float radius = maxRadius * sr.nextFloat();
                double angle = 360. * sr.nextDouble();

                float[] xy = calculateXAndYFromXcYcAndRadius(600.f, 400.f, radius, angle);

                xF[ii] = xy[0];
                yF[ii] = xy[1];
                x[ii] = (long)(Math.round(xy[0]));
                y[ii] = (long)(Math.round(xy[1]));
            }
            boolean useLinearSort = false;
            CHL ch = GrahamScan.computeHull(x, y, useLinearSort);
            
            int nH = ch.getXH().length;

            PairIntArray curve = new PairIntArray(nH - 1);
            float[] xHull = new float[nH];
            float[] yHull = new float[nH];
            for (int ii = 0; ii < nH; ++ii) {
                xHull[ii] = ch.getXH()[ii];
                yHull[ii] = ch.getYH()[ii];
                if (ii < (nH - 1)) {
                    curve.add((int)xHull[ii], (int)yHull[ii]);
                }
            }
            plotter.addPlot(xF, yF, xHull, yHull, "gs");
            plotter.writeFile();

        }
        
    }
    
    protected void createRandomPointsAroundCenter(int maxRadius,
        int numberOfPoints, int xc, int yc, PairIntArray xyPoints) throws 
        NoSuchAlgorithmException {

        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        sr.setSeed(System.nanoTime());

        for (int i = 0; i < numberOfPoints; i++) {

            float radius = maxRadius * sr.nextFloat();
            double angle = 360. * sr.nextDouble();

            float[] xy = calculateXAndYFromXcYcAndRadius((float)xc, (float)yc, 
                radius, angle);

            xyPoints.add(Math.round(xy[0]), Math.round(xy[1]));
        }
    }

    static float[] calculateXAndYFromXcYcAndRadiusCCW(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        return calculateXAndYFromXcYcAndRadius(xc, yc, radius, 360 - angleInDegreesFromYEQ0XGT0);
    }
    
    public static float[] calculateXAndYFromXcYcAndRadius(float xc, float yc, float radius, double angleInDegreesFromYEQ0XGT0) {

        double dx = radius * Math.cos(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));
        double dy = radius * Math.sin(angleInDegreesFromYEQ0XGT0 * (Math.PI/180.f));

        float x = (float) (xc + dx);
        float y = (float) (yc - dy);
        return new float[]{x, y};
    }

    /**
     * Test suite
     *
     * @return static Test
     */
    public static Test suite() {
        log.fine("Creating a TestSuite for GrahamScan");
        return new TestSuite(GrahamScanLongTest.class);
    }

    /**
     * Set up a Junit test runner
     *
     * @param args Not used.
     */
    public static void main(String[] args) {
        junit.textui.TestRunner.run(suite());
    }
}
