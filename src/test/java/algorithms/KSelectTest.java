package algorithms;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class KSelectTest extends TestCase {

    public void testTopK_int() {

        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int n = 1024;
        int[] a = new int[n];
        for (int i = 0; i <n; ++i) {
            a[i] = rand.nextInt();
        }

        int k = 43;

        int[] sorted = Arrays.copyOf(a, a.length);
        Arrays.sort(sorted);
        Map<Integer, Integer> expected = new HashMap<>();
        for (int i = n - 1; i > n - k - 1; --i) {
            if (expected.containsKey(sorted[i])) {
                expected.put(sorted[i], expected.get(sorted[i]) + 1);
            } else {
                expected.put(sorted[i], 1);
            }
        }

        int[] ans = KSelect.topK(a, k);
        for (int i = 0; i < k; ++i) {
            assert(expected.containsKey(ans[i]));
            if (expected.get(ans[i]) == 1) {
                expected.remove(ans[i]);
            } else {
                expected.put(ans[i], expected.get(ans[i]) - 1);
            }
        }
        assertTrue(expected.isEmpty());

    }
    public void testTopK_double() {
        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int n = 1024;
        double[] a = new double[n];
        for (int i = 0; i <n; ++i) {
            a[i] = rand.nextDouble();
        }

        int k = 43;

        double[] sorted = Arrays.copyOf(a, a.length);
        Arrays.sort(sorted);
        Map<Double, Integer> expected = new HashMap<>();
        for (int i = n - 1; i > n - k - 1; --i) {
            if (expected.containsKey(sorted[i])) {
                expected.put(sorted[i], expected.get(sorted[i]) + 1);
            } else {
                expected.put(sorted[i], 1);
            }
        }

        double[] ans = KSelect.topK(a, k);
        for (int i = 0; i < k; ++i) {
            assert(expected.containsKey(ans[i]));
            if (expected.get(ans[i]) == 1) {
                expected.remove(ans[i]);
            } else {
                expected.put(ans[i], expected.get(ans[i]) - 1);
            }
        }
        assertTrue(expected.isEmpty());
    }

    public void testBottomK_int() {

        long seed = System.nanoTime();
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);
        int n = 1024;
        int[] a = new int[n];
        for (int i = 0; i <n; ++i) {
            a[i] = rand.nextInt();
        }

        int k = 43;

        int[] sorted = Arrays.copyOf(a, a.length);
        Arrays.sort(sorted);
        Map<Integer, Integer> expected = new HashMap<>();
        for (int i = 0; i < k; ++i) {
            if (expected.containsKey(sorted[i])) {
                expected.put(sorted[i], expected.get(sorted[i]) + 1);
            } else {
                expected.put(sorted[i], 1);
            }
        }

        int[] ans = KSelect.bottomK(a, k);
        for (int i = 0; i < k; ++i) {
            assert(expected.containsKey(ans[i]));
            if (expected.get(ans[i]) == 1) {
                expected.remove(ans[i]);
            } else {
                expected.put(ans[i], expected.get(ans[i]) - 1);
            }
        }
        assertTrue(expected.isEmpty());

    }
    public void testBottomK_double() {

    }
}
