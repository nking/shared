package algorithms.signalProcessing;

import junit.framework.TestCase;

public class ResamplerTest extends TestCase {

    public void testResample() {
        int n1 = 3;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] expected = new double[][]{{0., 0.4, 0.8, 1.2, 1.6, 2.0}, {0., 0.4, 0.8, 1.2, 1.6, 2.0}};
        double[][] xyOut = Resampler.resample(xy, 2*n1);
        assertEquals(expected.length, xyOut.length);
        assertEquals(expected[0].length, xyOut[0].length);
        for (int i = 0; i < xyOut.length; ++i) {
            for (int j = 0; j < xyOut[0].length; ++j) {
                assertTrue(Math.abs(xyOut[i][j] - expected[i][j]) < 1E-11);
            }
        }
    }

    public void testUpDownSample() {
        int n1 = 12;
        int n2 = 16;
        //int lcm = 48;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] xyOut = Resampler.upDownSample(xy, n2);
        assertEquals(2, xyOut.length);
        assertEquals(n2, xyOut[0].length);

        assertTrue(Math.abs(xyOut[0][0] - xy[0][0]) < 1E-11);
        assertTrue(Math.abs(xyOut[1][0] - xy[1][0]) < 1E-11);

        assertTrue(xyOut[0][n2-1] <= xy[0][n1-1]);
        assertTrue(xyOut[1][n2-1] <= xy[1][n1-1]);

        // increasing
        double prev = -1;
        for (int i = 0; i < n2; ++i) {
            assertTrue(xyOut[0][i] > prev);
            assertTrue(xyOut[1][i] > prev);
            prev = xyOut[0][i];
        }

        double[][] xyOut2 = Resampler.upDownSample(xyOut, n1);
        assertEquals(2, xyOut2.length);
        assertEquals(n1, xyOut2[0].length);

        for (int i = 0; i < n1; ++i) {
            assertTrue(Math.abs(xy[0][i] - xyOut2[0][i]) < 1E-9);
            assertTrue(Math.abs(xy[1][i] - xyOut2[1][i]) < 1E-9);
        }

    }

}
