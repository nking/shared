package algorithms.optimization;

import junit.framework.TestCase;

import java.util.*;

public class LongestIncreasingSubsequenceTest extends TestCase {

    public void testFindAll() {
        int[] a;
        List<int[]> ans;

        a = new int[]{3,2,1,7,6,9,8};
        ans = LongestIncreasingSubsequence.findAllStrictlyIncreasing(a);

        /*
        Comparator<int[]> comp = new Comparator<int[]>() {
            public int compare(int[] o1, int[] o2) {
                if (Arrays.equals(o1, o2)) {
                    return 0;
                }
                return Arrays.compare(o1, o2);
            }
        };*/

        /*
         a = [3,2,1,7,6,9,8]
         we would have piles:
         pile 0   pile 1    pile 2
         3        7         9
         2        6         8
         1

         The length of the longest increasing subsequence is 3.

         Every LIS is (3,7,9), (3,7,8), (3,6,9), (3,6,8)
         (2,7,9), (2,7,8), (2,6,9), (2,6,8)
         (1,7,9), (1,7,8), (1,6,9), (1,6,8)

         There are 12 of LIS (3*2*2=12)

        by value:
        (3,7,9), (3,7,8), (3,6,9), (3,6,8)
         (2,7,9), (2,7,8), (2,6,9), (2,6,8)
         (1,7,9), (1,7,8), (1,6,9), (1,6,8)
         by index:
         (0,3,5), (0,3,6), (0,4,5), (0,4,6)
         (1,3,5), (1,3,6), (1,4,5), (1,4,6)
         (2,3,5), (2,3,6), (2,4,5), (2,4,6)
         */
        // cannot build with a comparator
        Set<int[]> expected = new HashSet<>();
        expected.add(new int[]{0,3,5});
        expected.add(new int[]{0,3,6});
        expected.add(new int[]{0,4,5});
        expected.add(new int[]{0,4,6});
        expected.add(new int[]{1,3,5});
        expected.add(new int[]{1,3,6});
        expected.add(new int[]{1,4,5});
        expected.add(new int[]{1,4,6});
        expected.add(new int[]{2,3,5});
        expected.add(new int[]{2,3,6});
        expected.add(new int[]{2,4,5});
        expected.add(new int[]{2,4,6});

        for (int[] ai : ans) {
            int[] m = null;
            for (int[] s : expected) {
                if (Arrays.equals(ai, s)) {
                    m = s;
                    break;
                }
            }
            assertNotNull(m);
            assertTrue(expected.remove(m));
        }
        assertTrue(expected.isEmpty());
    }

    public void testFindAll2() {
        int[] a;
        List<int[]> ans;

        a = new int[]{2,2,1,7,6};
        ans = LongestIncreasingSubsequence.findAllStrictlyIncreasing(a);

        /*
             2      7
             2
             1      6
             by value produces LIS:
             [2,7], [2,6], [2,7], [2,6], [1,7], [1,6]
             by index:
             [0,3], [0,4], [1,3], [1,4], [2,3], [2,4]
         */
        // cannot build with a comparator
        Set<int[]> expected = new HashSet<>();
        expected.add(new int[]{0, 3});
        expected.add(new int[]{0, 4});
        expected.add(new int[]{1, 3});
        expected.add(new int[]{1, 4});
        expected.add(new int[]{2, 3});
        expected.add(new int[]{2, 4});

        for (int[] ai : ans) {
            int[] m = null;
            for (int[] s : expected) {
                if (Arrays.equals(ai, s)) {
                    m = s;
                    break;
                }
            }
            assertNotNull(m);
            assertTrue(expected.remove(m));
        }
        assertTrue(expected.isEmpty());
    }

    public void testFindAny() {
        int[] a;
        List<Integer> ans;

        a = new int[]{3, 2, 1, 7, 6, 9, 8};
        ans = LongestIncreasingSubsequence.findAnyStrictlyIncreasing(a);

        Set<int[]> expected = new HashSet<>();
        expected.add(new int[]{0,3,5});
        expected.add(new int[]{0,3,6});
        expected.add(new int[]{0,4,5});
        expected.add(new int[]{0,4,6});
        expected.add(new int[]{1,3,5});
        expected.add(new int[]{1,3,6});
        expected.add(new int[]{1,4,5});
        expected.add(new int[]{1,4,6});
        expected.add(new int[]{2,3,5});
        expected.add(new int[]{2,3,6});
        expected.add(new int[]{2,4,5});
        expected.add(new int[]{2,4,6});

        //int[] ansI = ans.stream().collect(
        //       Collectors.toList()).stream().mapToInt(Integer::intValue).toArray();

        int[] ansI = ans.stream().mapToInt(Integer::intValue).toArray();

        int[] m = null;
        for (int[] s : expected) {
            if (Arrays.equals(ansI, s)) {
                m = s;
                break;
            }
        }
        assertNotNull(m);
        assertTrue(expected.remove(m));

        int[] ans2 = LongestIncreasingSubsequence.findSizeAndNumberStrictlyIncreasing(a);
        assertEquals(ans2[0], 3);
        assertEquals(ans2[1], 12);

    }

    public void testFindAllStrictlyIncreasing() {
        List<int[]> ans;
        int[][] a;

        a = new int[4][];
        a[0] = new int[]{4,1};
        a[1] = new int[]{2,3};
        a[2] = new int[]{2,2};
        a[3] = new int[]{6,3};

        ans = LongestIncreasingSubsequence.findAllStrictlyIncreasing(a);

        /*
        [[2,3], [6,3]], = [1,3]
        [[2,2], [6,3]], = [2,3]
        [[4,1], [6,3]]. = [0,3]
        */
        assertEquals(3, ans.size());
        Set<int[]> expected = new HashSet<>();
        expected.add(new int[]{1,3});
        expected.add(new int[]{2,3});
        expected.add(new int[]{0,3});

        for (int[] ai : ans) {
            int[] m = null;
            for (int[] s : expected) {
                if (Arrays.equals(ai, s)) {
                    m = s;
                    break;
                }
            }
            assertNotNull(m);
            assertTrue(expected.remove(m));
        }
        assertTrue(expected.isEmpty());

        List<Integer> ans2 = LongestIncreasingSubsequence.findAnyStrictlyIncreasing(a);
        assertEquals(2, ans2.size());
        assertEquals(0, ans2.get(0).intValue());
        assertEquals(3, ans2.get(1).intValue());
        int[] ans3 = LongestIncreasingSubsequence.findSizeAndNumberStrictlyIncreasing(a);
        assertEquals(2, ans3[0]);
        assertEquals(3, ans3[1]);

    }

}
