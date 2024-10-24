package algorithms;

import java.util.Arrays;

public class MedianOfMedians {

    /**
     find the median of array a with runtime complexity O(n) where n = a.length.
     <pre>
     references :
     CLRS Select, section 9.3

     </pre>
     * @param a
     * @return the median of array a
     */
    public static double medianOfMedians(double[] a) {
        return select(a, 0, a.length-1, a.length/2);
    }

    /**
     *
     * @param a
     * @param idxLo
     * @param idxHi
     * @param i the rank of the item to select
     * @return
     */
    static double select(double[] a, int idxLo, int idxHi, int i) {
        // TODO: transform to 0-based indexing
        while ((idxHi - idxLo + 1) %5 != 0) {
            for (int j = idxLo + 1; j <= idxHi; ++j) {
                if (a[idxLo] > a[j] && idxLo != j) {
                    double tmp = a[idxLo];
                    a[idxLo] = a[j];
                    a[j] = tmp;
                }
            }
            if (i == 1) {
                return a[idxLo];
            }
            ++idxLo;
            --i;
        }
        int g = (idxHi - idxLo + 1)/5;
        assert(idxLo + 5*g - 1 == idxHi);
        for (int j = idxLo; j <= (idxLo + g - 1); ++j) {
            quickSort5(a, j, j + g*4, g);
            // sort each group (A[j], A[j+g], A[j+2g], A[j+3g, A[j+4g] in place
            // length is 5.
            // sort in place algoriths: quicksort, insertsort, bubblesort, selectionsort
            // all but the first are n^2.  quicksort is O(n*log(n)) on average.
        }
        //debugPrint(a, idxLo, idxHi, g);
        // all group medians now lie in the middle fifth of A[idxLo:idcHi]
        // find the pivot x recursively as the median iof the group medians.
        // this range holds the medians:
        double x = select(a, idxLo + 2*g, idxLo + 3*g - 1, (int)Math.ceil(g/2.));

        // a is index of pivot x
        int q = partitionAround(a, idxLo, idxHi, x);
        //System.out.printf("pivotIdx=%d, pivot=%f\n", q, x);
        //debugPrint(a, idxLo, idxHi, g);

        // the result is just like lines 3-9 of randomized select
        int k = q - idxLo + 1;
        if (i == k) {
            return a[q];
        } else if (i < k) {
            return select(a, idxLo, q - 1, i);
        } else {
            return select(a, q + 1, idxHi, i-k);
        }
    }

    private static void debugPrint(double[] a, int idxLo, int idxHi, int g) {
        int idx = idxLo;
        System.out.printf("idxLo=%d, idxHi=%d, g=%d\n", idxLo, idxHi, g);
        for (int i = idxLo; i <= idxHi; i += g) {
            for (int j = 0; j < g; ++j) {
                System.out.printf("%6.1f ", a[idx++]);
            }
            System.out.println();
        }
        System.out.flush();;
    }

    protected  static void quickSort5(double[] a, int idxLo, int idxHi, int delta) {
        // sort each 5 elements at idxLo + i*delta where i=0 to 4, inclusive
        if (idxLo < idxHi) {
            int idxMid = partition5(a, idxLo, idxHi, delta);
            quickSort5(a, idxLo, idxMid - delta, delta);
            quickSort5(a, idxMid + delta, idxHi, delta);
        }
    }
    private static int partition5(double[] a, int idxLo, int idxHi, int delta) {
        final double x = a[idxHi];
        int i = idxLo - delta;
        double tmp;
        for (int k = idxLo; k < idxHi; k += delta) {
            if (a[k] <= x) {
                i += delta;
                if (i != k) {
                    tmp = a[i];
                    a[i] = a[k];
                    a[k] = tmp;
                }
            }
        }
        i += delta;
        if (i != idxHi) {
            tmp = a[i];
            a[i] = a[idxHi];
            a[idxHi] = tmp;
        }
        return i;
    }

    protected static int partitionAround(double[] a, int idxLo, int idxHi, double pivot) {
        // all values in a[idxLo, idxHi] arranged to < pivot and > pivot
        int i = idxLo - 1;
        double swap;
        for (int j = idxLo; j <= idxHi ; j++ ) {
            if (a[j] <= pivot) {
                i++;
                swap = a[i];
                a[i] = a[j];
                a[j] = swap;
            }
        }

        // find the pivotIdx as last matching value in a[idxLo, idxHi]
        int pivotIdx = -1;
        for (int j = idxHi; j >= 0 ; --j) {
            if (a[j] == pivot) {
                pivotIdx = j;
                break;
            }
        }
        if (pivotIdx == -1) {
            throw new IllegalArgumentException("did not find x in array a in range [idxLo, idxHi");
        }
        if (pivotIdx < i) {
            swap = a[i];
            a[i] = a[pivotIdx];
            a[pivotIdx] = swap;
            pivotIdx = i;
        }
        return pivotIdx;
    }
}
