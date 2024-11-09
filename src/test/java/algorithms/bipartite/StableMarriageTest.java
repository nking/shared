package algorithms.bipartite;

import junit.framework.TestCase;

import java.util.*;

public class StableMarriageTest extends TestCase {

    @SuppressWarnings({"unchecked", "rawtypes"})
    public void test1() {
        /*
        for CLRS chap 25.2
        test:
        a = Wanda, Emma, Lacey, Karen
             0      1      2     3
        b = Oscar, Davis, Brent, Hank
             0       1     2     3
        */
        LinkedList<Integer>[] aPrefs = new LinkedList[4];
        Map<Integer, Integer>[] bPrefs = new HashMap[4];
        for (int i = 0; i < 4; ++i) {
            aPrefs[i] = new LinkedList<Integer>();
            bPrefs[i] = new HashMap<Integer, Integer>();
        }
        aPrefs[0].add(2); aPrefs[0].add(3); aPrefs[0].add(0); aPrefs[0].add(1);
        aPrefs[1].add(1); aPrefs[1].add(3); aPrefs[1].add(0); aPrefs[1].add(2);
        aPrefs[2].add(2); aPrefs[2].add(1); aPrefs[2].add(3); aPrefs[2].add(0);
        aPrefs[3].add(2); aPrefs[3].add(3); aPrefs[3].add(1); aPrefs[3].add(0);

        bPrefs[0].put(0, 0); bPrefs[0].put(3, 1); bPrefs[0].put(2,2); bPrefs[0].put(1,3);
        bPrefs[1].put(0, 0); bPrefs[1].put(2, 1); bPrefs[1].put(3,2); bPrefs[1].put(1,3);
        bPrefs[2].put(2, 0); bPrefs[2].put(3,1); bPrefs[2].put(0,2); bPrefs[2].put(1,3);
        bPrefs[3].put(2, 0); bPrefs[3].put(0,1); bPrefs[3].put(1,2); bPrefs[3].put(3,3);

        int[] m0 = new int[]{3, 0, 2, 1};

        int[] m = StableMarriage.match(aPrefs, bPrefs);

        assertTrue(Arrays.equals(m0, m));
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    public void test2() {
        /*
        for CLRS chap 25.2
        test:
        a = Monica, Phoebe, Rachel
             0        1      2
        b = Chandler, Joey, Ross
             0         1     2
        */
        LinkedList<Integer>[] aPrefs = new LinkedList[3];
        Map<Integer, Integer>[] bPrefs = new HashMap[3];
        for (int i = 0; i < 3; ++i) {
            aPrefs[i] = new LinkedList<Integer>();
            bPrefs[i] = new HashMap<Integer, Integer>();
        }
        aPrefs[0].add(0); aPrefs[0].add(1); aPrefs[0].add(2);
        aPrefs[1].add(1); aPrefs[1].add(2); aPrefs[1].add(0);
        aPrefs[2].add(2); aPrefs[2].add(0); aPrefs[2].add(1);

        bPrefs[0].put(1, 0); bPrefs[0].put(2, 1); bPrefs[0].put(0, 2);
        bPrefs[1].put(2, 0); bPrefs[1].put(0, 1); bPrefs[1].put(1, 2);
        bPrefs[2].put(0, 0); bPrefs[2].put(1, 1); bPrefs[2].put(2, 2);

        int[] m0 = new int[]{0, 1, 2};
        int[] m1 = new int[]{2, 0, 1};
        int[] m2 = new int[]{1, 2, 0};

        int[] m = StableMarriage.match(aPrefs, bPrefs);

        assertTrue(Arrays.equals(m0, m) || Arrays.equals(m1, m) || Arrays.equals(m2, m));
    }
}
