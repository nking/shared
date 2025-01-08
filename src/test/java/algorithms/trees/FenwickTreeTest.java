package algorithms.trees;

import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

public class FenwickTreeTest extends TestCase {

    public void testPositive0() {

        long seed = System.nanoTime();
        Random rand = new Random(seed);

        int n = 10;
        int upper = 100;

        long[] a = new long[n];
        long[] b = new long[n];
        int i;
        for (i = 0; i < n; ++i) {
            a[i] = rand.nextInt(upper);
            if (rand.nextBoolean()) {
                b[i] = rand.nextInt(upper);
            }
        }

        run(a, b, true);
    }

    public void _testPositive1() {

        long seed = System.nanoTime();
        Random rand = new Random(seed);

        int n = 10;
        int upper = 100;

        long[] a = new long[n+1];
        long[] b = new long[n+1];
        int i;
        for (i = 1; i < a.length; ++i) {
            a[i] = rand.nextInt(upper);
            if (rand.nextBoolean()) {
                b[i] = rand.nextInt(upper);
            }
        }

        run(a, b, false);
    }

    public void _test0() {

        long seed = System.nanoTime();
        Random rand = new Random(seed);

        int n = 10;

        long[] a = new long[n];
        long[] b = new long[n];
        int i;
        for (i = 0; i < n; ++i) {
            a[i] = rand.nextInt();
            if (rand.nextBoolean()) {
                b[i] = rand.nextInt();
            }
        }

        run(a, b, true);
    }

    public void _test1() {

        long seed = System.nanoTime();
        //seed = 1234567L;
        System.out.printf("seed=%d\n", seed);
        Random rand = new Random(seed);

        int n = 10;

        long[] a = new long[n+1];
        long[] b = new long[n+1];
        int i;
        for (i = 1; i < a.length; ++i) {
            a[i] = rand.nextInt();
            if (rand.nextBoolean()) {
                b[i] = rand.nextInt();
            }
        }

        run(a, b, false);
    }

    private void run(long[] a, long[] b, boolean use0Based) {
        int i;
        int n = a.length;
        int start = use0Based ? 0 : 1;

        FenwickTreeLong tree = new FenwickTreeLong(n, use0Based);
        // test set,get
        for (i = start; i < n; ++i) {
            tree.set(i, a[i]);
        }
        for (i = start; i < n; ++i) {
            assertEquals(a[i], tree.get(i));
        }

        // test add,get
        tree = new FenwickTreeLong(n, use0Based);
        for (i = start; i < n; ++i) {
            tree.add(i, a[i]);
        }
        for (i = start; i < n; ++i) {
            assertEquals(a[i], tree.get(i));
        }

        // test constructor
        tree = new FenwickTreeLong(a, use0Based);
        for (i = start; i < n; ++i) {
            assertEquals(a[i], tree.get(i));
        }

        // test construct, add, get
        tree = new FenwickTreeLong(a, use0Based);
        for (i = start; i < n; ++i) {
            tree.add(i, b[i]);
        }
        for (i = start; i < n; ++i) {
            assertEquals(a[i] + b[i], tree.get(i));
        }

        // test construct, set, get
        tree = new FenwickTreeLong(a, use0Based);
        for (i = start; i < n; ++i) {
            tree.set(i, b[i]);
        }
        for (i = start; i < n; ++i) {
            assertEquals(b[i], tree.get(i));
        }

    }
}
