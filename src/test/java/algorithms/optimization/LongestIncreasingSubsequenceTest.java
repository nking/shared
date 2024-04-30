package algorithms.optimization;

import junit.framework.TestCase;

import java.util.*;
import java.util.stream.Collectors;

public class LongestIncreasingSubsequenceTest extends TestCase {

    public void testFindGTAll() {
        int[] a;
        List<int[]> ans;

        a = new int[]{3,2,1,7,6,9,8};
        ans = LongestIncreasingSubsequence.findGTAll(a);

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

    public void testFindGTAll2() {
        int[] a;
        List<int[]> ans;

        a = new int[]{2,2,1,7,6};
        ans = LongestIncreasingSubsequence.findGTAll(a);

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

    public void testFindStrictlyAll2() {
        int[] a;
        List<int[]> ans;

        a = new int[]{2,2,1,7,6};
        ans = LongestIncreasingSubsequence.findStrictlyIncrAll(a);

        /*
             patience piles by value:
             pile 0   pile 1   pile 2
             2          2        7
                        1        6

         resulting sequences by value:
         [2,2,7], [2,2,6], [2,1,7], [2,1,6]

         resulting sequences by index:
         [0,1,3], [0,1,4], [0,2,3], [0,2,4]
         */
        // cannot build with a comparator
        Set<int[]> expected = new HashSet<>();
        expected.add(new int[]{0, 1, 3});
        expected.add(new int[]{0, 1, 4});
        expected.add(new int[]{0, 2, 3});
        expected.add(new int[]{0, 2, 4});

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

    public void testFindGTAny() {
        int[] a;
        List<Integer> ans;

        a = new int[]{3, 2, 1, 7, 6, 9, 8};
        ans = LongestIncreasingSubsequence.findGTAny(a);

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

        int size = LongestIncreasingSubsequence.findGTSize(a);
        assertEquals(3, size);

        int[] ans2 = LongestIncreasingSubsequence.findGTSizeAndNumber(a);
        assertEquals(ans2[0], size);
        assertEquals(ans2[1], 12);

    }

}
