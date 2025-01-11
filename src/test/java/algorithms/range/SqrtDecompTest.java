package algorithms.range;

import junit.framework.TestCase;

public class SqrtDecompTest extends TestCase {

    public void test0() {

        boolean use0BasedIndexes  = true;
        long[] a = new long[]{5,8,6,3,2,7,2,6};

        SqrtDecomp sd = new SqrtDecomp(a);
        long s = sd.sum(2, 7, use0BasedIndexes);
        System.out.printf("s=%d, expect=%d\n", s, 26);
        assertEquals(26, s);

        int[][] updates = new int[][]{
                {2,4,3},
                {1,5,4},
                {5,7,1},
                {3,3,3}
        };
        sd.updateAdd(updates, use0BasedIndexes);

        s = sd.sum(0, 7, use0BasedIndexes);
        //System.out.printf("s=%d, expect=%d\n", s, 74);
        assertEquals(74, s);

        s = sd.sum(4, 6, use0BasedIndexes);
        //System.out.printf("s=%d, expect=%d\n", s, 24);
        assertEquals(24, s);

    }
}
