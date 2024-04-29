package algorithms.search;

import java.util.Arrays;
import java.util.List;

/**
 * useful methods for bisecting search, a.k.a. binary search
 */
public class MiscBisectingSearch {

    /**
     * finds where srch would be inserted such that all keys at indexes < resulting index
     * have smaller values than srch.   if srch is larger than all array elements,
     * this method returns an index larger than last array index.
     <pre>

     a = [2, 2]
     srch = 1,
     returns index 0

     a = [2, 2]
     srch = 2,
     returns index 2

     a = [2, 2]
     srch = 2,
     returns index 2
     </pre>
     * @param a
     * @param srch
     * @return smallest index i such that a[i] is GT srch.
     * range of returned value is is [0, n] where n is a.length.
     */
    public static int successor(int[] a, int srch) {
        int lo = 0;
        int hi = a.length - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a[mid] <= srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    }

    /**
     * find the least index i whose value a]i] EQ srch, else finds the least index i
     * whose value a[i] is GT srch.
     * The method returns indexes in the range [-1, a.length], inclusive.
     * The return valye of a.length means that srch is larger than all values in the array.
     *
     * Note that one could instead use the java.util.Arrays.binarySearch(...) and if the index is
     * negative, you can convert it with mult by -1 and then subtract 1 (
     * i.e. idx = -(insertion point) - 1, so insertion point = -idx - 1).
     * @param a
     * @param srch
     * @return
     * range is [0, n] inclusive where n = a.length.
     */
    public static int findInsertIndex(int[] a, int srch) {
        int lo = 0;
        int hi = a.length - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            // NOTE: the successor method uses <= here:
            if (a[mid] < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    }

    /**
     * find the least index i whose value a]i] EQ srch, else finds the least index i
     * whose value a[i] is GT srch.
     * The method returns indexes in the range [-1, a.length], inclusive.
     * The return valye of a.length means that srch is larger than all values in the array.
     *
     * Note that one could instead use the java.util.Arrays.binarySearch(...) and if the index is
     * negative, you can convert it with mult by -1 and then subtract 1 (
     * i.e. idx = -(insertion point) - 1, so insertion point = -idx - 1).
     * @param a
     * @param srch
     * @return
     * range is [0, n] inclusive where n = a.length.
     */
    public static int findInsertIndex(List<Integer> a, int srch) {
        int lo = 0;
        int hi = a.size() - 1;
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            // NOTE: the successor method uses <= here:
            if (a.get(mid) < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return lo;
    }

    /**
     * finds where srch would be inserted such that all keys at indexes < resulting index
     * have smaller values than srch.   if srch is larger than all array elements,
     * this method returns an index larger than last array index.
     <pre>

     a = [2, 2]
     srch = 1,
     returns index 0

     a = [2, 2]
     srch = 2,
     returns index 2

     a = [2, 2]
     srch = 2,
     returns index 2
     </pre>
     * @param a
     * @param srch
     * @return smallest index i such that a[i] is GT srch.
     * range of returned value is is [0, n] where n is a.length.
     */
    public static int successor(List<Integer> a, int srch) {
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

    /**
     Finds the smallest index for which value == srch, or if doesn't exist,
     finds the largest index for which the value is LT srch.

     * This is a floor function.
     *
     * If the method returns -1, the srch is smaller than all elements in the
     * array.
     * if srch is larger than every element in list 'a', this method will return
     * the last index of a.
     <pre>
     e.g.

     a = [0,1,2,2,3]
     srch = 2
     returns  2

     a = [0,3,4,4,6};
     srch = 5;
     returns 3

     a = [2,3,4,4,6};
     srch = 0;
     returns -1

     a = [2,3,4,4,6};
     srch = 7;
     returs 4

     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @return the floor index, that is,
     * finds the smallest index for which value == srch, or if doesn't exist,
     *      finds the largest index for which the value is LT srch.
     *      range of values from floor is [-1, a.length-1] inclusive.
     */
    public static int floor(int[] a, int srch) {
        return floor(a, srch, 0, a.length - 1);
    }

    /**
     Finds the smallest index for which value == srch, or if doesn't exist,
     finds the largest index for which the value is LT srch.

     * This is a floor function.
     *
     * If the method returns -1, the srch is smaller than all elements in the
     * array.
     * if srch is larger than every element in list 'a', this method will return
     * the last index of a.
     <pre>
     e.g.

     a = [0,1,2,2,3]
     srch = 2
     returns  2

     a = [0,3,4,4,6};
     srch = 5;
     returns 3

     a = [2,3,4,4,6};
     srch = 0;
     returns -1

     a = [2,3,4,4,6};
     srch = 7;
     returs 4

     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @param lo smallest index for search range
     * @param hi largets index for search range, inclusive.
     * @return the floor index, that is,
     * finds the smallest index for which value == srch, or if doesn't exist,
     *      finds the largest index for which the value is LT srch.
     *      range of values from floor is [-1, a.length-1] inclusive.
     */
    public static int floor(int[] a, int srch, int lo, int hi) {
        if (srch < a[lo]) {
            return -1;
        }
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a[mid] == srch) {
                //we want the smallest index whose value == srch, so search further in lower half of array
                if (mid == lo) {
                    return mid;
                }
                int idx = floor(a, srch, lo, mid - 1);
                if (idx == -1 || a[idx] < srch) {
                    return mid;
                }
                return idx;
            } else if (a[mid] < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return Math.min(lo, hi);
    }

    /**
     Finds the smallest index for which value == srch, or if doesn't exist,
     finds the largest index for which the value is LT srch.

     * This is a floor function.
     *
     * If the method returns -1, the srch is smaller than all elements in the
     * array.
     * if srch is larger than every element in list 'a', this method will return
     * the last index of a.
     <pre>
     e.g.

     a = [0,1,2,2,3]
     srch = 2
     returns  2

     a = [0,3,4,4,6};
     srch = 5;
     returns 3

     a = [2,3,4,4,6};
     srch = 0;
     returns -1

     a = [2,3,4,4,6};
     srch = 7;
     returs 4

     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @return the floor index, that is,
     * finds the smallest index for which value == srch, or if doesn't exist,
     *      finds the largest index for which the value is LT srch.
     */
    public static int floor(List<Integer> a, int srch) {
        return floor(a, srch, 0, a.size() - 1);
    }
    public static int floor(List<Integer> a, int srch, int lo, int hi) {
        if (srch < a.get(lo)) {
            return -1;
        }
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a.get(mid) == srch) {
                //we want the smallest index whose value == srch, so search further in lower half of array
                int idx = floor(a, srch, lo, mid - 1);
                if (idx == -1 || a.get(idx) < srch) return mid;
                return idx;
            } else if (a.get(mid) < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return Math.min(lo, hi);
    }

    /**
     * find the largest index i whose value a[i] EQ srch,
     * else if no value exists this method returns the smallest index whose
     * value is GT srch.  if all values in the array are less than srch, the
     * method returns index == array size which is beyond bounds of array.
     <pre>
     for example:
     a       = [1, 2, 3, 3, 4]
     indexes =  0  1  2  3  4

     ceil of srch=3 is index 3
     ceil of srch=4 is index 4
     ceil of srch=9 is index 5
     ceil of srch=0 is index 0
     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @return
    the greatest index i for which a[i] == srch, else if no such value exists,
    the index for the smallest value GT srch.
    range is [0, n] inclusive where n = a.length.
     */
    public static int ceiling(List<Integer> a, int srch) {
        return ceiling(a, srch, 0, a.size() - 1);
    }

    /**
     * find the largest index i whose value a[i] EQ srch,
     * else if no value exists this method returns the smallest index whose
     * value is GT srch.  if all values in the array are less than srch, the
     * method returns index == array size which is beyond bounds of array.
     <pre>
     for example:
         a       = [1, 2, 3, 3, 4]
         indexes =  0  1  2  3  4

         ceil of srch=3 is index 3
         ceil of srch=4 is index 4
         ceil of srch=9 is index 5
         ceil of srch=0 is index 0
     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @param lo smallest index of search range
     * @param hi largest index of search range
     * @return
    the greatest index i for which a[i] == srch, else if no such value exists,
    the index for the smallest value GT srch.
    range is [0, n] inclusive where n = a.length.
     */
    public static int ceiling(List<Integer> a, int srch, int lo, int hi) {
        //System.out.printf("a=%s, srch=%d, lo=%d, hi=%d\n",
        //        Arrays.toString(a.stream().mapToInt(Integer::intValue).toArray()),
        //        srch, lo, hi);
        if (srch > a.get(hi)) {
            return a.size();
        } else if (srch < a.get(lo)) {
            return lo;
        }
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a.get(mid) == srch) {
                if (mid + 1 > hi){
                    return mid;
                }
                //we want the largest index whose value == srch, so search further in higher half of array
                /*int idx = ceiling(a, srch, mid+1, hi);
                if (idx == -1 || idx == hi+1 || a.get(idx) > srch) return mid;
                return idx;*/
                lo = mid + 1;
            } else if (a.get(mid) < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return Math.min(lo, hi);
    }

    /**
     find the largest index i whose value a[i] EQ srch,
     * else if no value exists this method returns the smallest index whose
     * value is GT srch.  if all values in the array are less than srch, the
     * method returns index == array size which is beyond bounds of array.
     <pre>
     for example:
     a       = [1, 2, 3, 3, 4]
     indexes =  0  1  2  3  4

     ceil of srch=3 is index 3
     ceil of srch=4 is index 4
     ceil of srch=9 is index 5
     ceil of srch=0 is index 0
     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @param lo smallest index of search range
     * @param hi largest index of search range
     * @return
     the greatest index i for which a[i] == srch, else if no such value exists,
     the index for the smallest value GT srch.
     range is [0, n] inclusive where n = a.length.
     */
    public static int ceiling(int[] a, int srch) {
        return ceiling(a, srch, 0, a.length - 1);
    }

    /**
     * find the largest index i whose value a[i] EQ srch,
     * else if no value exists this method returns the smallest index whose
     * value is GT srch.  if all values in the array are less than srch, the
     * method returns index == array size which is beyond bounds of array.
     <pre>
     for example:
     a       = [1, 2, 3, 3, 4]
     indexes =  0  1  2  3  4

     ceil of srch=3 is index 3
     ceil of srch=4 is index 4
     ceil of srch=9 is index 5
     ceil of srch=0 is index 0
     </pre>
     * @param a an ascending ordered list (non-decreasing, adjacent values can be ==).
     * @param srch the value to search for in a
     * @param lo smallest index of search range
     * @param hi largest index of search range
     * @return the greatest index i for which a[i] == srch, else if no such value exists,
     * the index for the smallest value GT srch.
     * range is [0, n] inclusive where n = a.length.
     */
    public static int ceiling(int[] a, int srch, int lo, int hi) {

        if (srch > a[hi]) {
            return a.length;
        } else if (srch < a[lo]) {
            return lo;
        }
        while (lo <= hi) {
            int mid = lo + (hi - lo) / 2;
            if (a[mid] == srch) {
                if (mid + 1 > hi) {
                    return mid;
                }
                //we want the largest index whose value == srch, so search further in higher half of array
                /*int idx = ceiling(a, srch, mid+1, hi);
                if (idx == -1 || idx == hi+1 || a[idx] > srch) {
                    return mid;
                }
                return idx;*/
                lo = mid + 1;
            } else if (a[mid] < srch) {
                lo = mid + 1;
            } else {
                hi = mid - 1;
            }
        }
        return Math.min(lo, hi);
    }

    /**
     * find the greatest index i whose value a[i] is LT srch.
     * @param a list sorted by ascending values
     * @param srch value to srch for predecessor of.
     * @return the greatest index i whose value a[i] is LT srch.
     * range is [-1, n-1] where n = a.length.
     */
    public static int predecessor(int[] a, int srch) {

        // range of values from floor is [-1, a.length-1]
        int floorIdx = MiscBisectingSearch.floor(a, srch);

        if (floorIdx == -1) {
            return floorIdx;
        }
        if (a[floorIdx] == srch) {
            return floorIdx - 1;
        }
        return floorIdx;
    }

    /**
     * find the greatest index i whose value a[i] is LT srch.
     * @param a list sorted by ascending values
     * @param srch value to srch for predecessor of.
     * @return the greatest index i whose value a[i] is LT srch.
     * range is [-1, n-1] where n = a.length.
     */
    public static int predecessor(List<Integer> a, int srch) {

        // range of values from floor is [-1, a.length-1]
        int floorIdx = MiscBisectingSearch.floor(a, srch);

        if (floorIdx == -1) return floorIdx;
        if (a.get(floorIdx) == srch) return floorIdx - 1;
        return floorIdx;
    }

}
