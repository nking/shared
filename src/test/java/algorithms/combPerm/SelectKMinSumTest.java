package algorithms.combPerm;

import junit.framework.TestCase;

public class SelectKMinSumTest extends TestCase {

    public void test0() {
        int[][] kByN = new int[][] {
                {2,    3,  4,  5,  1,  3},// 3
                {2+2,3+2,4+2,5+2,1+2,3+2},// 5
                {2+4,3+4,4+4,5+4,1+4,3+4}, //6
                {2+5,3+5,4+5,5+5,1+5,3+5}  //6
        };
        long expAns = 20;
        long ans = SelectKMinSum.selectKMinSum(kByN);
        assertEquals(expAns, ans);
    }
}
