package algorithms.search;

import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MiscBisectingSearchTest extends TestCase {

    public void testFindInsertIndex() {
        int[] a;
        int ans, srch, expAns;
        //System.out.println("findInsertIndex1");

        a = new int[]{0,1,2,2,3};
        srch = 2;
        expAns = 2;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        //System.out.printf("srch=%d, a=%s\n  ret=%d\n", srch, Arrays.toString(a), ans);
        assertEquals(expAns, ans);

        a = new int[]{0,1,2,2,2,3};
        srch = 2;
        expAns = 2;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2,2,2,3,3,5,7,8};
        srch = 2;
        expAns = 0;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{-2,-2,0,0,0,1,1,2,2,2};
        srch = 2;
        expAns = 7;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2, 2};
        srch = 2;
        expAns = 0;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2, 2};
        srch = 1;
        expAns = 0;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2, 2};
        srch = 5;
        expAns = 2;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{0,3,4,4,6};
        srch = 5;
        expAns = 4;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2,3,4,4,6};
        srch = 0;
        expAns = 0;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2,3,4,4,6};
        srch = 7;
        expAns = 5;
        ans = MiscBisectingSearch.findInsertIndex(a, srch);
        assertEquals(expAns, ans);

    }
    public void testFloor() {
        // all indexes <= returned index hold values LEQ srch key

        int[] a;
        int expAns, ans, srch;

        a = new int[]{0,1,2,2,3};
        srch = 2;
        expAns = 2;
        ans = MiscBisectingSearch.floor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.floor(toList(a), srch);
        assertEquals(expAns, ans);

        a = new int[]{0,3,4,4,6};
        srch = 5;
        expAns = 3;
        ans = MiscBisectingSearch.floor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.floor(toList(a), srch);
        assertEquals(expAns, ans);

        a = new int[]{2,3,4,4,6};
        srch = 0;
        expAns = -1;
        ans = MiscBisectingSearch.floor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.floor(toList(a), srch);
        assertEquals(expAns, ans);

        a = new int[]{2,3,4,4,6};
        srch = 7;
        expAns = 4;
        ans = MiscBisectingSearch.floor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.floor(toList(a), srch);
        assertEquals(expAns, ans);

    }

    public void testSuccessorForIncreasingList() {

        // all indexes < returned index hold values LT srch key

        int[] a;
        int expAns, ans, srch;

        a = new int[]{2, 2};
        srch = 1;
        expAns = 0;

        ans = MiscBisectingSearch.successor(toList(a), srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.successor(a, srch);
        assertEquals(expAns, ans);

        srch = 2;
        expAns = 2; // beyond array indexes
        ans = MiscBisectingSearch.successor(toList(a), srch);
        assertEquals(expAns, ans);

        srch = 5;
        expAns = 2; // beyond array indexes
        ans = MiscBisectingSearch.successor(toList(a), srch);
        assertEquals(expAns, ans);

        a = new int[]{0,1,2,2,3};
        srch = 2;
        expAns = 4;
        ans = MiscBisectingSearch.successor(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{0,3,4,4,6};
        srch = 5;
        expAns = 4;
        ans = MiscBisectingSearch.successor(a, srch);
        assertEquals(expAns, ans);

        a = new int[]{2,3,4,4,6};
        srch = 0;
        expAns = 0;
        ans = MiscBisectingSearch.successor(a, srch);
        assertEquals(expAns, ans);
    }

    public void testCeiling() {

        // all indexes < returned index hold values LT srch key

        int[] a;
        int expAns, ans, srch;

        a = new int[]{1, 2, 3, 3, 4};
        srch = 3;
        expAns = 3;
        ans = MiscBisectingSearch.ceiling(toList(a), srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.ceiling(a, srch);
        assertEquals(expAns, ans);

        srch = 4;
        expAns = 4;
        ans = MiscBisectingSearch.ceiling(toList(a), srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.ceiling(a, srch);
        assertEquals(expAns, ans);

        srch = 10;
        expAns = 5;
        ans = MiscBisectingSearch.ceiling(toList(a), srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.ceiling(a, srch);
        assertEquals(expAns, ans);

        srch = 0;
        expAns = 0;
        ans = MiscBisectingSearch.ceiling(toList(a), srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.ceiling(a, srch);
        assertEquals(expAns, ans);
    }

    public void testPredecessor() {
        int[] a;
        int expAns, ans, srch;

        a = new int[]{1, 2, 3, 3, 4};
        srch = 3;
        expAns = 1;
        ans = MiscBisectingSearch.predecessor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.predecessor(toList(a), srch);
        assertEquals(expAns, ans);

        srch = 4;
        expAns = 3;
        ans = MiscBisectingSearch.predecessor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.predecessor(toList(a), srch);
        assertEquals(expAns, ans);

        srch = 10;
        expAns = 4;
        ans = MiscBisectingSearch.predecessor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.predecessor(toList(a), srch);
        assertEquals(expAns, ans);

        srch = 0;
        expAns = -1;
        ans = MiscBisectingSearch.predecessor(a, srch);
        assertEquals(expAns, ans);
        ans = MiscBisectingSearch.predecessor(toList(a), srch);
        assertEquals(expAns, ans);
    }

    protected List<Integer> toList(int[] a) {
        List<Integer> list = new ArrayList<>();
        for (int i : a) {
            list.add(i);
        }
        return list;
    }
}
