package algorithms.range;

import algorithms.misc.MiscMath0;

import junit.framework.TestCase;

import java.util.*;

public class MosAlgorithmTest extends TestCase {

    public void testSums() {

        long seed = System.nanoTime();
        //seed = 268028218057513L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 10;
        int nMax = 64;
        int nMin = 3;

        int nQueriesBound;
        for (int nTest = 0; nTest < nTests;++nTest) {

            int n = nMin + rand.nextInt(nMax - nMin);
            int[] a = new int[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();
            }

            long[] prefixA = prefix(a);

            long nQNCK = MiscMath0.computeNDivKTimesNMinusK(n, 2);

            int nQ = 1 + rand.nextInt((int)nQNCK);
            int[][] qs = new int[nQ][];
            int i0, i1;
            for (int j = 0; j < nQ; ++j) {
                i0 = rand.nextInt(n); // 9; 10-9=1-1
                i1 = i0 + rand.nextInt(n - i0);
                qs[j] = new int[]{i0, i1};
            }
            long[] sums = MosAlgorithm.querySums(a, qs);
            for (int j = 0; j < nQ; ++j) {
                assertEquals(prefixA[qs[j][1] + 1] - prefixA[qs[j][0]], sums[j]);
            }
        }
    }

    protected long[] prefix(int[] a) {
        long[] out = new long[a.length + 1];
        for (int i = 0; i < a.length; ++i) {
            out[i+1] = a[i];
        }
        for (int i = 1; i < out.length; ++i) {
            out[i] += out[i-1];
        }
        return out;
    }

    public void testFrequencies() {

        long seed = System.nanoTime();
        //seed = 268028218057513L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int nTests = 10;
        int nMax = 64;
        int nMin = 3;

        for (int nTest = 0; nTest < nTests;++nTest) {

            int n = nMin + rand.nextInt(nMax - nMin);
            int[] a = new int[n];
            for (int i = 0; i < n; ++i) {
                a[i] = rand.nextInt();
            }

            long[] prefixA = prefix(a);

            long nQNCK = MiscMath0.computeNDivKTimesNMinusK(n, 2);

            int nQ = 1 + rand.nextInt((int)nQNCK);
            int[][] qs = new int[nQ][];
            int i0, i1;
            for (int j = 0; j < nQ; ++j) {
                i0 = rand.nextInt(n); // 9; 10-9=1-1
                i1 = i0 + rand.nextInt(n - i0);
                qs[j] = new int[]{i0, i1};
            }

            List<Map<Integer, Integer>> frequencies = MosAlgorithm.queryFrequencies(a, qs);

            for (int j = 0; j < nQ; ++j) {
                // brute force check
                Map<Integer, Integer> expected = new HashMap<>();
                for (int idx = qs[j][0]; idx <= qs[j][1]; ++idx){
                    expected.put(a[idx], expected.getOrDefault(idx, 0) + 1);
                }
                Map<Integer, Integer> f = frequencies.get(j);
                assertEquals(expected.size(), f.size());
                for (int key : expected.keySet()) {
                    assertEquals(expected.get(key), f.get(key));
                    f.remove(key);
                }
                assertTrue(f.isEmpty());
            }
        }
    }
}
