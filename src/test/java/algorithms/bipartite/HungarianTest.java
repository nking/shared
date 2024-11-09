package algorithms.bipartite;

import algorithms.matrix.MatrixUtil;
import junit.framework.TestCase;
import no.uib.cipr.matrix.Matrix;

import java.util.Arrays;
import java.util.Random;

public class HungarianTest extends TestCase {

    public void test0() {
        double[][] w = makeGraph1();
        double expSum = 37;
        int[][] expM = getExpectedMGraph1();

        Hungarian h = new Hungarian(w);
        int[][] m = h.solve();

        assertEquals(expM.length, m.length);

        for (int i = 0; i < m.length; ++i) {
            assertEquals(expM[i][0], m[i][0]);
            assertEquals(expM[i][1], m[i][1]);
        }
    }

    public void test1() {

        long seed = System.nanoTime();
        //seed = 28747942308221L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        // if generate a random solution and want to check best answer by brute force,
        // need to keep nL, nR <= 9 for n! <= 362880

        int nTests = 10;
        for (int nTest = 0; nTest < nTests; ++nTest) {
            // for the brute force I'm using, need nL==nR
            int nL = 4 + rand.nextInt(6);
            int nR = nL;//4 + rand.nextInt(6);

            double[][] w = makeGraph0(nL, nR, rand);
            int[][] expM = getExpectedMGraph0(w);
            double expSum = sum(w, expM);

            Hungarian h = new Hungarian(w);
            int[][] m = h.solve();
            double s = sum(w, m);

            //System.out.printf("expSum=%f, res=%f\n", expSum, s);

            assertEquals(expM.length, m.length);

            assertTrue(Math.abs(expSum - s) < 1E-11);

            //TODO: fix getExpectedMGraph0 to return all maximal matchings
            for (int i = 0; i < m.length; ++i) {
                assertEquals(expM[i][0], m[i][0]);
                assertEquals(expM[i][1], m[i][1]);
            }
        }
    }

    protected static class Res {
        double sum = Double.NEGATIVE_INFINITY;
        int[] m;
        public void setIfLarger(int[] m, double[][] w) {
            double s = sum(w, m);
            if (s > this.sum) {
                this.sum = s;
                this.m = Arrays.copyOf(m, m.length);
            }
        }
    }

    protected int[][] getExpectedMGraph0(double[][] w) {
        if (w.length != w[0].length) {
            throw new IllegalArgumentException("cannot use this method unless w.length=w[0].length");
        }
        int[] m = new int[w.length];
        for (int i = 0; i < w.length; ++i) {
            m[i] = i;
        }
        Res res = new Res();
        recursion(w.length, res, m, w);

        int[][] outM = new int[w.length][];
        for (int i = 0; i < w.length; ++i) {
            outM[i] = new int[]{i, res.m[i]};
        }
        return outM;
    }
    protected static double sum(double[][] w, int[] m) {
        double s = 0;
        for (int i = 0; i < m.length; ++i) {
            s += w[i][m[i]];
        }
        return s;
    }
    protected static double sum(double[][] w, int[][] m) {
        double s = 0;
        for (int i = 0; i < m.length; ++i) {
            s += w[m[i][0]][m[i][1]];
        }
        return s;
    }
    protected void recursion(int k, Res res, int[] m, final double[][] w) {
        if (k == 0) {
            res.setIfLarger(m, w);
            return;
        }
        for (int i = 0; i < k; ++i) {

            recursion(k-1, res, m, w);

            if (i < k - 1) {
                if ((k&1)!=1) {// even
                    if (i != k-1) {
                        // swap
                        m[i] ^= m[k-1];
                        m[k-1] ^= m[i];
                        m[i] ^= m[k-1];
                    }
                } else {//odd
                    if (0 != k-1) {
                        // swap
                        m[0] ^= m[k-1];
                        m[k-1] ^= m[0];
                        m[0] ^= m[k-1];
                    }
                }
            }
        }
    }

    protected double[][] makeGraph0(int nL, int nR, Random rand) {
        boolean transpose = (nR > nL);
        if (transpose) {
            nL ^= nR;
            nR ^= nL;
            nL ^= nR;
        }
        double[][] w = new double[nL][nR];
        for (int vL = 0; vL < nL; ++vL) {
            for (int vR = 0; vR < nR; ++vR) {
                w[vL][vR] = rand.nextInt(1000) * rand.nextDouble();
            }
        }

        if (transpose) {
            return MatrixUtil.transpose(w);
        }
        return w;
    }

    protected int[][] getExpectedMGraph1() {
        int[][] expM = new int[4][];
        expM[0] = new int[]{0, 0};
        expM[1] = new int[]{1, 2};
        expM[2] = new int[]{2, 1};
        expM[3] = new int[]{3, 3};
        return expM;
    }

    protected double[][] makeGraph1() {
        /*
        0 : 0 = 11
          : 1 = 12
        1 : 2 = 9
        2:1 = 10
        2:2 = 11
        3:3 = 7
        expSum = 37
         */
        double[][] w = new double[4][4];
        w[0][0] = 11;
        w[0][1] = 12;
        w[1][2] = 9;
        w[2][1] = 10;
        w[2][2] = 11;
        w[3][3] = 7;
        return w;

    }
}
