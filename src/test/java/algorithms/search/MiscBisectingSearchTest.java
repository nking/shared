package algorithms.search;

import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

public class MiscBisectingSearchTest extends TestCase {

    public void testRightIncr() {

        int[] a;
        int eRes, res, srch;

        a = new int[]{2, 2};
        eRes = 0;
        srch = 1;
        res = MiscBisectingSearch.bisectRightForIncreasingList(toList(a), srch);
        assertEquals(eRes, res);

        eRes = 2;
        srch = 2;
        res = MiscBisectingSearch.bisectRightForIncreasingList(toList(a), srch);
        assertEquals(eRes, res);
    }

    public void testLeftIncr() {
        int[] a;
        int eRes, res, srch;

        a = new int[]{2, 2};
        eRes = 0;
        srch = 1;
        res = MiscBisectingSearch.bisectLeftForIncreasingList(toList(a), srch);
        assertEquals(eRes, res);

        eRes = 0;
        srch = 2;
        res = MiscBisectingSearch.bisectLeftForIncreasingList(toList(a), srch);
        assertEquals(eRes, res);
    }

    protected List<Integer> toList(int[] a) {
        List<Integer> list = new ArrayList<>();
        for (int i : a) {
            list.add(i);
        }
        return list;
    }
}
