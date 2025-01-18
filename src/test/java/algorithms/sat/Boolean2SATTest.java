package algorithms.sat;

import junit.framework.TestCase;

import java.util.Map;

public class Boolean2SATTest extends TestCase {

    public void test1() {

        int[][] l1Clauses = new int[][] {
                {2, -1}, {-1, -2}, {1, 3}, {-2, -3}, {1, 4}
        };
        int m = 4;

        Boolean2SAT sat = new Boolean2SAT();
        Map<Integer, Boolean> soln = sat.solve(l1Clauses, m);
        assertNotNull(soln);
        assertEquals(false, soln.get(1).booleanValue());
        assertEquals(false, soln.get(2).booleanValue());
        assertEquals(true, soln.get(3).booleanValue());
        assertEquals(true, soln.get(4).booleanValue());

    }

    public void test2() {

        int[][] l2Clauses = new int[][] {
                {1,2},{1,-2},{-1,3},{-1,-3}
        };
        int m = 3;

        Boolean2SAT sat = new Boolean2SAT();
        Map<Integer, Boolean> soln = sat.solve(l2Clauses, m);
        assertNull(soln);

    }
}
