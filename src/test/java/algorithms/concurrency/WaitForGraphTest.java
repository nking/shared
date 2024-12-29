package algorithms.concurrency;

import junit.framework.TestCase;

public class WaitForGraphTest extends TestCase {

    public void test0() {

        WaitForGraph wg = new WaitForGraph();

        for (int i = 0; i < 3; ++i) {
            assertTrue(wg.cannotDeadlock(i, i));
            wg.addAlloc(i, i);
        }
        assertTrue(wg.addRequestIfCannotDeadlock(1, 2));

        assertTrue(wg.addRequestIfCannotDeadlock(2, 0));

        assertFalse(wg.addRequestIfCannotDeadlock(0, 1));

    }
}
