package algorithms.search;

import java.util.List;

/**
 * useful methods for bisecting search, a.k.a. binary search
 */
public class MiscBisectingSearch {

    /**
     * bisecting search tailored to find the position in list 'a',
     * such that all elements at indexes < the returned position are LEQ srch.
     *
     * if srch is larger than every element in ;ost 'a', this method will return a.size()
     * which is beyond the indexes of list 'a'.
     * if this method returns a value < a.size(), it is where this item can
     * be inserted or can replace an item.
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @return the position in list 'a',
     * such that all elements at indexes < the returned position are LEQ srch.
     */
    public static int bisectRightForIncreasingList(List<Integer> a, int srch) {
        int lo = 0;
        int hi = a.size() - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a.get(mid) <= srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    }

    public static int bisectLeftForIncreasingList(List<Integer> a, int srch) {
        return bisectLeftForIncreasingList(a, srch, 0, a.size() - 1);
    }
    public static int bisectLeftForIncreasingList(List<Integer> a, int srch, int lo, int hi) {
        int mid;
        while (lo <= hi) {
            mid = (lo + hi)/2;
            if (a.get(mid) < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return Math.max(lo, hi);
    }
}
