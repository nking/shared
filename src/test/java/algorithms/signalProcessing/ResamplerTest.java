package algorithms.signalProcessing;

import algorithms.matrix.MatrixUtil;
import algorithms.util.PairIntArray;
import algorithms.util.PolygonAndPointPlotter;
import junit.framework.TestCase;

import java.io.IOException;

public class ResamplerTest extends TestCase {

    public void testResample1() throws IOException {
        int n1 = 120;
        int n2 = 174;
        float[][] pxy = new float[2][120];
        pxy[0] = new float[]{215.000f,215.000f,215.000f,216.000f,218.000f,220.000f,222.000f,224.000f,225.000f,225.000f,225.000f,227.000f,226.000f,226.000f,225.000f,225.000f,224.000f,223.000f,222.000f,220.000f,218.000f,216.000f,215.000f,216.000f,217.000f,219.000f,221.000f,223.000f,225.000f,227.000f,229.000f,231.000f,233.000f,235.000f,237.000f,239.000f,241.000f,243.000f,244.000f,246.000f,248.000f,250.000f,252.000f,254.000f,256.000f,258.000f,260.000f,260.000f,258.000f,256.000f,254.000f,252.000f,250.000f,249.000f,248.000f,247.000f,247.000f,247.000f,247.000f,248.000f,249.000f,250.000f,251.000f,252.000f,254.000f,256.000f,258.000f,260.000f,261.000f,261.000f,261.000f,260.000f,258.000f,256.000f,254.000f,252.000f,250.000f,248.000f,248.000f,248.000f,249.000f,250.000f,250.000f,251.000f,253.000f,255.000f,257.000f,259.000f,257.000f,255.000f,253.000f,251.000f,249.000f,247.000f,245.000f,243.000f,241.000f,239.000f,237.000f,235.000f,233.000f,231.000f,229.000f,227.000f,225.000f,223.000f,221.000f,223.000f,225.000f,226.000f,226.000f,226.000f,227.000f,227.000f,225.000f,223.000f,221.000f,219.000f,217.000f,216.000f};
        pxy[1] = new float[]{147.000f,149.000f,151.000f,153.000f,155.000f,156.000f,157.000f,159.000f,161.000f,163.000f,165.000f,166.000f,168.000f,170.000f,172.000f,174.000f,176.000f,178.000f,180.000f,181.000f,181.000f,183.000f,185.000f,187.000f,189.000f,189.000f,189.000f,188.000f,187.000f,187.000f,186.000f,184.000f,183.000f,183.000f,185.000f,186.000f,188.000f,190.000f,191.000f,192.000f,193.000f,193.000f,193.000f,193.000f,192.000f,191.000f,189.000f,187.000f,185.000f,184.000f,182.000f,181.000f,180.000f,178.000f,176.000f,174.000f,172.000f,170.000f,168.000f,166.000f,164.000f,162.000f,161.000f,159.000f,157.000f,156.000f,155.000f,154.000f,152.000f,150.000f,148.000f,146.000f,144.000f,143.000f,143.000f,143.000f,143.000f,141.000f,139.000f,137.000f,135.000f,133.000f,131.000f,129.000f,127.000f,127.000f,127.000f,127.000f,127.000f,127.000f,127.000f,127.000f,127.000f,127.000f,125.000f,124.000f,123.000f,122.000f,121.000f,121.000f,122.000f,123.000f,125.000f,127.000f,128.000f,128.000f,128.000f,129.000f,131.000f,133.000f,135.000f,137.000f,139.000f,141.000f,143.000f,143.000f,143.000f,143.000f,144.000f,146.000f};
        float[][] xy = CurveResampler.resample(pxy, n2);
        int t = 2;
    }
    public void testUpsample1() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);

        double[][] xy2 = CurveResampler.resample(xy, xy[0].length);
        assertEquals(xy[0].length, xy2[0].length);
        for (int i = 0; i < xy[0].length; ++i) {
            assertTrue(Math.abs(xy[0][i] - xy2[0][i]) < 1E-11);
            assertTrue(Math.abs(xy[1][i] - xy2[1][i]) < 1E-11);
        }

        float[][] xy3 = MatrixUtil.convertToFloat(xy);
        float[][] xy4 = CurveResampler.resample(xy3, xy3[0].length);
        assertEquals(xy3[0].length, xy3[0].length);
        for (int i = 0; i < xy3[0].length; ++i) {
            assertTrue(Math.abs(xy3[0][i] - xy4[0][i]) < 1E-9);
            assertTrue(Math.abs(xy3[1][i] - xy4[1][i]) < 1E-9);
        }
    }

    public void testDownsample1() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);

        double[][] xy2 = CurveResampler.resample(xy, xy[0].length-5);
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xy2[0], xy2[1], null, null, "downsample1");
        plot.writeFile("downsample1");

    }

    public void testUpSample2() throws IOException {
        double[][] xy = getTriangle(7, 2, 1);
        int n1 = xy[0].length;
        int n2 = 30;
        double[][] xyR = CurveResampler.resample(xy, n2);
        PolygonAndPointPlotter plot = new PolygonAndPointPlotter();
        plot.addPlot(xyR[0], xyR[1], null, null, "resampled2");
        plot.writeFile("tri_resampled");
        plot = new PolygonAndPointPlotter();
        plot.addPlot(xy[0], xy[1], null, null, "triangle2");
        plot.writeFile("tri");

        double[][] xyUpDown = CurveResampler.resample(xyR, n1);
        plot = new PolygonAndPointPlotter();
        plot.addPlot(xyUpDown[0], xyUpDown[1], null, null, "up_down2");
        plot.writeFile("up_down2");

        float[][] xyUpDown2 = CurveResampler.resample(MatrixUtil.convertToFloat(xyR), n1);
        plot = new PolygonAndPointPlotter();
        double[][] xyUpDown3 = MatrixUtil.convertToDouble(xyUpDown2);
        plot.addPlot(xyUpDown3[0], xyUpDown3[1], null, null, "up_down2 float");
        plot.writeFile("up_down2 float");
    }

    public void testResample3() {
        int n1 = 3;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] expected = new double[][]{{0., 0.4, 0.8, 1.2, 1.6, 2.0}, {0., 0.4, 0.8, 1.2, 1.6, 2.0}};
        double[][] xyOut = CurveResampler.resample(xy, 2*n1);
        assertEquals(expected.length, xyOut.length);
        assertEquals(expected[0].length, xyOut[0].length);
        for (int i = 0; i < xyOut.length; ++i) {
            for (int j = 0; j < xyOut[0].length; ++j) {
                assertTrue(Math.abs(xyOut[i][j] - expected[i][j]) < 1E-11);
            }
        }
    }

    public void testUpDownSample4() {
        int n1 = 12;
        int n2 = 16;
        //int lcm = 48;
        double[][] xy = new double[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        double[][] up = CurveResampler.resample(xy, n2);
        assertEquals(2, up.length);
        assertEquals(n2, up[0].length);
        double[][] upDown = CurveResampler.resample(up, n1);
        assertEquals(2, upDown.length);
        assertEquals(n1, upDown[0].length);

        int[][] xyInt = new int[2][n1];
        for (int i = 0; i < n1; ++i) {
            xy[0][i] = i;
            xy[1][i] = i;
        }
        int[][] upInt = CurveResampler.resample(xyInt, n2);
        int[][] upDownInt = CurveResampler.resample(upInt, n1);
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
