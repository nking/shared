package algorithms;

import algorithms.util.FormatArray;

public class MedianOfMedians {

    /**
     * NOTE: this implementation is not correct.

     find the median of array a with runtime complexity O(n) where n = a.length.
     <pre>
     references :
     CLRS Select, section 9.3
     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a
     * @return the median of array a
     */
    public static double medianOfMediansCLRS(double[] a) {
        return selectCLRS(a, 0, a.length-1, a.length/2);
    }

    /**
     NOTE: this implementation is not correct.

     find the median of array a with runtime complexity O(n) where n = a.length.
     <pre>
     references :
     CLRS Select, section 9.3
     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a
     * @param idxLo smallest index to of range
     * @param idxHi largest index of range, inclusive
     * @param i the rank of the item to select in 0-based numbering.
     * @return the median
     */
    static double selectCLRS(double[] a, int idxLo, int idxHi, int i) {
        double result = select(a, idxLo, idxHi, i);
        // by this time, the array is largely sorted.
        System.out.printf("Do these agree? %f, and a[i]=%f?\n", result, a[i]);
        return result;
        //return a[i];
    }

    /**
     NOTE: this one is not correct.

     find the median of array a with runtime complexity O(n) where n = a.length.
     <pre>
     references :
     CLRS Select, section 9.3
     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a
     * @param idxLo smallest index to of range
     * @param idxHi largest index of range, inclusive
     * @param i the rank of the item to select in 0-based numbering.
     * @return the median
     */
    private static double select(double[] a, int idxLo, int idxHi, int i) {

        System.out.printf("*select idxLo=%d, idxHi=%d, i=%d\n", idxLo, idxHi, i);

        final int _idxLo = idxLo;
        int _idxHi = idxHi;
        final int _i = i;
        while ((idxHi - idxLo + 1) %5 != 0) {
            for (int j = idxLo + 1; j <= idxHi; ++j) {
                if (a[idxLo] > a[j] && idxLo != j) {
                    double tmp = a[idxLo];
                    a[idxLo] = a[j];
                    a[j] = tmp;
                }
            }
            if (i == 0 && (idxHi - idxLo) < 5) {
                System.out.printf("return 0\n");
                return a[idxLo];
            }
            ++idxLo;
            --i;
        }
        int g = (idxHi - idxLo + 1)/5;

        System.out.printf(" => idxLo=%d, idxHi=%d, i=%d, g=%d\n", idxLo, idxHi, i, g);

        if (idxHi < idxLo) {
            // idxHi + i or idxLo - 1 + i?
            System.out.printf("   ERROR idxHi<idxLo a=%s\n", FormatArray.toString(a, "%.0f"));
            //return a[idxLo + i];
        }
        if (g==0) {
            System.out.printf("   ** a=%s\n", FormatArray.toString(a, "%.0f"));
            //return a[idxLo + i];
        }

        assert(idxLo + 5*g - 1 == idxHi);
        for (int j = idxLo; j <= (idxLo + g - 1); ++j) {
            quickSort5(a, j, j + g*4, g);
            // sort each group (A[j], A[j+g], A[j+2g], A[j+3g, A[j+4g] in place
            // length is 5.
            // sort in place algoriths: quicksort, insertsort, bubblesort, selectionsort
            // all but the first are n^2.  quicksort is O(n*log(n)) on average.
        }
        debugPrint(a, idxLo, idxHi, g);

        //Blum et al. 1973 1.b, pick recursively if n/5 > 1
        if (g < 2 && (idxHi - idxLo) < 5) {
            System.out.printf("   **** i=%d, idxLo=%d, idxHi=%d\n    a=%s\n", i, idxLo, idxHi, FormatArray.toString(a, "%.0f"));
            return a[idxLo + i];
        }

        System.out.printf("NEXT select 1\n");
        System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        // all group medians now lie in the middle fifth of A[idxLo:idcHi]
        // find the pivot x recursively as the median iof the group medians.
        // this range holds the medians:
        // but if idxLo was offset above, this possibly needs to be lowered up to 2 "rows"

        int nextI = (int)Math.ceil(g/2);//((idxHi - idxLo) + 4) / 5;
        double x = select(a, idxLo + 2*g, idxLo + 3*g - 1, nextI);

        // q is index of pivot x, 0-based
        int q = partitionAround(a, idxLo, idxHi, x);

        int k = q - idxLo + 1;

        System.out.printf("pivotIdx=q=%d, pivot=%.0f, k=%d\n", q, x, k);
        System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        if (k==i) {
            System.out.printf("NEXT select 2\n");
            return a[q];
        } else if (k>i) {
            System.out.printf("NEXT select 3\n");
            return select(a, idxLo, q - 1, i);
        } else {
            System.out.printf("NEXT select 4\n");
            //return select(a, q + 1, idxHi, i-k);
            return select(a, q + 1, idxHi, i-k);
        }
    }

    private static void debugPrint(double[] a, int idxLo, int idxHi, int g) {
        int idx = idxLo;
        System.out.printf("idxLo=%d, idxHi=%d, g=%d\n", idxLo, idxHi, g);
        for (int i = idxLo; i <= idxHi; i += g) {
            for (int j = 0; j < g; ++j) {
                System.out.printf("%6.1f ", a[idx++]);
                if (idxHi == idx) break;
            }
            System.out.println();
            if (g == 0 || idxHi == idx) break;
        }
        System.out.flush();
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
