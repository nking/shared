package algorithms.compGeometry;

import java.util.Random;
import java.util.logging.Logger;
import junit.framework.TestCase;

public class LinesAndAnglesTest extends TestCase {

    protected Logger log = Logger.getLogger(this.getClass().getSimpleName());

    public LinesAndAnglesTest(String testName) {
        super(testName);
    }

    @Override
    protected void setUp() throws Exception {
        super.setUp();
    }

    @Override
    protected void tearDown() throws Exception {
        super.tearDown();
    }

    public void testCrossProduct() {

        log.info("crossProduct");

        /*         |
         *       4 -
         *         |
         *       3 -  o2
         *         |
         * 1 o   2 -
         *         |
         *       1 -
         *         |
         *       0 o--|--|--|--|--|--|
         *         0  1  2  3  4
         *
         */
        double x1 = -3;
        double y1 = 2;
        double x2 = 1;
        double y2 = 3;
        double result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result < 0); // clockwise
        assertTrue(result == -11);

        x1 = 3;
        result = LinesAndAngles.crossProduct(x1, y1, x2, y2);
        assertTrue(result >= 0); // counter clockwise
        assertTrue(result == 7);
    }
    
    public void testDirectionCCW() throws Exception {
        
        /*
         * 
         * 10                             p2
         *  9
         *  8                          /
         *  7     p3
         *  6      \
         *  5                     /
         *  4        \
         *  3
         *  2           \     /
         *  1              p1
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * directionCCW of change from P1:P2 to P1:P3
         */
        
        float x2 = 10;
        float y2 = 10;
        float x1 = 5;
        float y1 = 1;
        float x3 = 2;
        float y3 = 7;
        
        double direction = LinesAndAngles.directionCCW(x1, y1, x2, y2, x3, y3);
        
        assertTrue(direction > 0);
        
        // swap p1 and p3
        x3 = 10;
        y3 = 10;
        x2 = 2;
        y2 = 7;
        direction = LinesAndAngles.directionCCW(x1, y1, x2, y2, x3, y3);
        
        assertTrue(direction < 0);

    }

    public void testLinesIntersect() {

        /*
                        x2,y2
             x3,y3
         x1,y1      x4,y4
         */
        float x1 = 1;
        float y1 = 1;
        float x2 = 5;
        float y2 = 4;
        float x3 = 2;
        float y3 = 2;
        float x4 = 3;
        float y4 = 1;

        boolean intersect = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertTrue(intersect);

        x3=5;
        y3=5;
        intersect = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertTrue(intersect);

        x3=6;
        y3=5;
        intersect = LinesAndAngles.linesIntersect(x1, y1, x2, y2, x3, y3, x4, y4);
        assertFalse(intersect);

    }
    
    public void testDirection_CCW() throws Exception {
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2                 
         *  1              p2
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * directionCCW of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives negative number
        int direction = LinesAndAngles.directionCW(
            2, 4, 
            5, 1, 
            3, 7);
        assertTrue(direction < 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2                 
         *  1     p2       
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * directionCCW of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives negative number
        direction = LinesAndAngles.directionCW(
            2, 4, 
            2, 1, 
            3, 7);
        assertTrue(direction < 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                      
         *  4     p1     
         *  3
         *  2  p2             
         *  1         
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * directionCCW of change from P1:P2 to P1:P3
         */
        
        // CCW sweep of P1:P2 to P1:P3 gives positive number
        direction = LinesAndAngles.directionCW(
            2, 4, 
            1, 2, 
            3, 7);
        assertTrue(direction > 0);
        
        /*
         * 
         * 10                             
         *  9
         *  8                           
         *  7        p3
         *  6       
         *  5                    p1  
         *  4            
         *  3
         *  2                 
         *  1              p2
         *  0
         *   0  1  2  3  4  5  6  7  8  9  10
         * 
         * directionCCW of change from P1:P2 to P1:P3
         */
        // CCW sweep of P1:P2 to P1:P3 gives positive number
        direction = LinesAndAngles.directionCW(
            7, 5, 
            5, 1, 
            3, 7);
        assertTrue(direction > 0);
       
    }
    
    public void testPointIsInLine() throws Exception {
        
        boolean isInLine = LinesAndAngles.pointIsInLine(0, 5, 0, 0, 0, 10);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(1, 5, 0, 0, 0, 10);
        assertFalse(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(0.f, 5.f, 0.f, 0.f, 0.f, 10.f);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(2, 2, 0, 0, 4, 4);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(3, 0, 2, 0, 4, 0);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(3, 1, 2, 0, 4, 0);
        assertFalse(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(2, 0, 2, 0, 4, 0);
        assertTrue(isInLine);
        
        isInLine = LinesAndAngles.pointIsInLine(4, 4, 0, 0, 4, 4);
        assertTrue(isInLine);
    }
    
    public void testCalcClockwiseAngle() {
        
        float x1, y1, x2, y2, x3, y3;
        double angle, expected; 
        double eps = 0.001;
        
        /*
           1     2
              3
        */
        x1 = 0;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 20;
        y2 = 10;
        expected = Math.PI/2;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3   2
        */
        x1 = 0;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 20;
        y2 = 0;
        expected = 3*Math.PI/4;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
                 2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 20;
        y2 = 0;
        expected = Math.PI;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);

        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
              2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 10;
        y2 = 0;
        expected = 5.*Math.PI/4.;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
           1     
              3 
           2
        */
        x1 = 0;
        y1 = 20;
        x3 = 10;
        y3 = 10;
        x2 = 0;
        y2 = 0;
        expected = 6.*Math.PI/4.;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
              1     
           2  3 
        */
        x1 = 10;
        y1 = 10;
        x3 = 10;
        y3 = 0;
        x2 = 0;
        y2 = 0;
        expected = 6.*Math.PI/4.;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);
        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);
        
        /*
              3  2
              1
        */
        x1 = 0;
        y1 = 0;
        x3 = 0;
        y3 = 10;
        x2 = 10;
        y2 = 10;
        expected = 6.*Math.PI/4.;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1, 
            x2, y2, x3, y3);

        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);

        /*
              3
                 2
              1
        */
        x1 = 0;
        y1 = 0;
        x3 = 0;
        y3 = 1;
        x2 = 0 + (float)(Math.sqrt(2)/2.);
        y2 = 1 - (float)(Math.sqrt(2)/2.);
        expected = 7.*Math.PI/4.;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1,
                x2, y2, x3, y3);

        //System.out.println("a=" + angle + " expected=" + expected);
        assertTrue(Math.abs(angle - expected) < eps);

        x1=11.000000f;
        y1=2.000000f;
        x2=9.000000f;
        y2=2.000000f;
        x3=10.000000f;
        y3=2.000000f;
        expected = Math.PI;
        angle = LinesAndAngles.calcClockwiseAngle(x1, y1,
                x2, y2, x3, y3);
        assertTrue(Math.abs(angle - expected) < eps);
    }
    
    public void testCalcPolarPerp() {
        
        int[] ep;
        int w, h, t;
        
        int radius = 4;
        
        w = 5; h = 5; t = 0;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        
        assertEquals(radius, ep[0]);
        assertEquals(0, ep[1]);
        assertEquals(radius, ep[2]);
        assertEquals(h - 1, ep[3]);
        
        w = 5; h = 5; t = 180;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        
        assertEquals(-radius, ep[0]);
        assertEquals(0, ep[1]);
        assertEquals(-radius, ep[2]);
        assertEquals(h - 1, ep[3]);
        
        w = 5; h = 5; t = 270;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        
        assertEquals(0, ep[0]);
        assertEquals(-radius, ep[1]);
        assertEquals(w - 1, ep[2]);
        assertEquals(-radius, ep[3]);
        
        w = 5; h = 5; t = 90;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        
        assertEquals(0, ep[0]);
        assertEquals(radius, ep[1]);
        assertEquals(w - 1, ep[2]);
        assertEquals(radius, ep[3]);
        
        // ----- 30 degrees ------
        w = 5; h = 10; t = 30;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(8, ep[1]);
        assertEquals(5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 5; h = 5; t = 30;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(2, ep[0]);
        assertEquals(4, ep[1]);
        assertEquals(5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 4; h = 5; t = 30;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(2, ep[0]);
        assertEquals(4, ep[1]);
        assertEquals(3, ep[2]);
        assertEquals(3, ep[3]);
       
        w = 4; h = 10; t = 30;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(8, ep[1]);
        assertEquals(3, ep[2]);
        assertEquals(3, ep[3]);
      
        // ----- 60 degrees
        w = 10; h = 5; t = 60;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(5, ep[1]);
        assertEquals(8, ep[2]);
        assertEquals(0, ep[3]);
        
        // --- 150 degrees ---
        w = 5; h = 10; t = 150;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(8, ep[1]);
        assertEquals(-5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 5; h = 5; t = 150;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(-2, ep[0]);
        assertEquals(4, ep[1]);
        assertEquals(-5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 4; h = 5; t = 150;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(-2, ep[0]);
        assertEquals(4, ep[1]);
        assertEquals(-3, ep[2]);
        assertEquals(3, ep[3]);
       
        w = 4; h = 10; t = 150;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(8, ep[1]);
        assertEquals(-3, ep[2]);
        assertEquals(3, ep[3]);
        
        // --- 210 degrees ---
        w = 5; h = 10; t = 210;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(-8, ep[1]);
        assertEquals(-5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 5; h = 5; t = 210;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(-2, ep[0]);
        assertEquals(-4, ep[1]);
        assertEquals(-5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 4; h = 5; t = 210;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(-2, ep[0]);
        assertEquals(-4, ep[1]);
        assertEquals(-3, ep[2]);
        assertEquals(-3, ep[3]);
       
        w = 4; h = 10; t = 210;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(-8, ep[1]);
        assertEquals(-3, ep[2]);
        assertEquals(-3, ep[3]);
        
        // --- 330 degrees ---
        w = 5; h = 10; t = 330;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(-8, ep[1]);
        assertEquals(5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 5; h = 5; t = 330;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(2, ep[0]);
        assertEquals(-4, ep[1]);
        assertEquals(5, ep[2]);
        assertEquals(0, ep[3]);
        
        w = 4; h = 5; t = 330;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(2, ep[0]);
        assertEquals(-4, ep[1]);
        assertEquals(3, ep[2]);
        assertEquals(-3, ep[3]);
       
        w = 4; h = 10; t = 330;
        ep = LinesAndAngles.calcPolarLineEndPoints(t, radius, w, h);
        //System.out.println(Arrays.toString(ep));
        assertEquals(0, ep[0]);
        assertEquals(-8, ep[1]);
        assertEquals(3, ep[2]);
        assertEquals(-3, ep[3]);
        
    }

    public void testAreaPolygon() {
        double[] x = new double[]{
            2,5,7,4,4
        };
        double[] y = new double[]{
            4,5,3,1,3
        };
        double expArea = 17./2;
        double area = LinesAndAngles.areaPolygon(x, y);
        assertTrue(Math.abs(expArea - area) < 1E-7);
    }

    public void testAreaTriangle() {
        // areaTriangle
        //double x1, double y1, double x2, double y2, double x3, double y3)
        /*

         */

        double x2 = 2;
        double y2 = 2;

        double x1 = x2+2;
        double y1 = y2+2;

        double x3 = x2 + 6;
        double y3 = y2;

        double expArea = 2 + ((4*2)/2);

        double[] x = new double[]{x1, x2, x3};
        double[] y = new double[]{y1, y2, y3};
        double area1 = LinesAndAngles.areaPolygon(x, y);

        double area2 = LinesAndAngles.areaTriangle(x1, y1, x2, y2, x3, y3);

        assertTrue(Math.abs(expArea - area1) < 1E-7);
        assertTrue(Math.abs(expArea - area2) < 1E-7);
        double expDist = 2;
        double minDist = LinesAndAngles.minDistPointToLine(x1, y1,
                x2,y2,x3,y3);
        assertTrue(Math.abs(expDist - minDist) < 1E-7);
    }

    public void testCalcAngle() {
        long seed = System.nanoTime();
        System.out.printf("seed = %d\n", seed);
        Random rand = new Random(seed);

        /*
        P3:P1  sweeping clockwise to P3:P2
              P1  P2
                P3
        */

        final double x3 = 0.f;
        final double y3 = 0.f;
        final double x2 = 1.f;
        final double y2 = 0.f;
        int offsetX, offsetY;
        double _x1, _y1, _x2, _y2, _x3, _y3, _angleR, angleR;

        for (double angle = 0; angle < 360; angle+=1) {
            offsetX = 0;//rand.nextInt(100);
            offsetY = 0;//rand.nextInt(100);
            _x3 = x3 + offsetX;
            _y3 = y3 + offsetY;
            _x2 = x2 + offsetX;
            _y2 = y2 + offsetY;

            _x1 = Math.cos(angle * Math.PI/180.) + offsetX;
            _y1 = Math.sin(angle * Math.PI/180.) + offsetY;

            _angleR = LinesAndAngles.calcAngle(_x1, _y1, _x2, _y2, _x3, _y3);
            angleR = angle * Math.PI/180.;

            if (angle > 180) {
                _angleR += 2.*Math.PI;
            }

            assertTrue(Math.abs(angleR - _angleR) < 1E-6);
        }
    }
}
