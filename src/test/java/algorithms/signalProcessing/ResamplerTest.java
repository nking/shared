package algorithms.signalProcessing;

import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

import java.io.IOException;

public class ResamplerTest extends TestCase {

    public void testUpsample1() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);

        double[][] xy2 = CurveResampler.linearUpsample(xy, xy[0].length);
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xy2[0], xy2[1], null, null, "upsample1");
        plot.writeFile("upsample1");

        //TODO: assert they're the same.
    }

    public void testDownsample1() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);

        double[][] xy2 = CurveResampler.linearDownsample(xy, xy[0].length-5);
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xy2[0], xy2[1], null, null, "downsample1");
        plot.writeFile("downsample1");
    }

    public void testUpSample2() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);
        int n1 = xy[0].length;
        int n2 = 30;
        double[][] xyR = CurveResampler.linearUpsample(xy, n2);
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xyR[0], xyR[1], null, null, "resampled2");
        plot.writeFile("tri_resampled");
        plot = new PolygonAndPointPlotter();
        plot.addPlot(xy[0], xy[1], null, null, "triangle2");
        plot.writeFile("tri");

        double[][] xyUpDown = CurveResampler.upDownSample(xy, n2);
        plot = new PolygonAndPointPlotter();
        plot.addPlot(xyUpDown[0], xyUpDown[1], null, null, "up_down2");
        plot.writeFile("up_down2");

    }

    public void _testResample() {
        int n1 = 3;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] expected = new double[][]{{0., 0.4, 0.8, 1.2, 1.6, 2.0}, {0., 0.4, 0.8, 1.2, 1.6, 2.0}};
        double[][] xyOut = CurveResampler.linearUpsample(xy, 2*n1);
        assertEquals(expected.length, xyOut.length);
        assertEquals(expected[0].length, xyOut[0].length);
        for (int i = 0; i < xyOut.length; ++i) {
            for (int j = 0; j < xyOut[0].length; ++j) {
                assertTrue(Math.abs(xyOut[i][j] - expected[i][j]) < 1E-11);
            }
        }
    }

    public void _testUpDownSample() {
        int n1 = 12;
        int n2 = 16;
        //int lcm = 48;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] xyOut = CurveResampler.upDownSample(xy, n2);
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

        double[][] xyOut2 = CurveResampler.upDownSample(xyOut, n1);
        assertEquals(2, xyOut2.length);
        assertEquals(n1, xyOut2[0].length);

        for (int i = 0; i < n1; ++i) {
            assertTrue(Math.abs(xy[0][i] - xyOut2[0][i]) < 1E-9);
            assertTrue(Math.abs(xy[1][i] - xyOut2[1][i]) < 1E-9);
        }

    }

    protected double[][] getTriangle(int topY, int baseY, int leftX) {
        /*
        7                    *
        6                 *     *
        5              *           *
        4           *                 *
        3        *                       *
        2     *  *  *  *  *  *  *  *  *  *  *
        1
        0
           0  1  2  3  4  5  6  7  8  9 10 11
        */
        PairIntArray p = new PairIntArray();
        int i = leftX;
        int j = baseY;
        while (j <= topY) {
            p.add(i, j);
            ++i;
            ++j;
        }
        j -= 2;
        while (j >= baseY) {
            p.add(i, j);
            ++i;
            --j;
        }
        i -= 2;
        ++j;
        while (i > leftX) {
            p.add(i, j);
            --i;
        }

        double[][] xy = new double[2][p.getN()];
        for (i = 0; i < p.getN(); ++i) {
            xy[0][i] = p.getX(i);
            xy[1][i] = p.getY(i);
        }

        return xy;
    }
}
