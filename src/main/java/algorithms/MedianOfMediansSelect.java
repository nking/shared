package algorithms;

import java.util.Arrays;

public class MedianOfMediansSelect {

    /**
     find the median of array a with runtime complexity O(n) where n = a.length.
     The worst case runtime complexity is O(n).

     <pre>
     references :
     Select algorithm, section 9.3
     "Introduction to Algorithms" by
     Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein

     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a and unsorted array
     * @return the median of array a
     */
    static double medianOfMediansCLRS(double[] a) {
        return selectCLRS(a, 0, a.length-1, a.length/2);
    }

    /**
     find the value with rank i in array a with runtime complexity O(n) where n = a.length where i as a rank is 0-based.
     The worst case runtime complexity is O(n).

     <pre>
     references :
     Select, section 9.3
     "Introduction to Algorithms" by
     Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein

     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a
     * @param idxLo smallest index to of range
     * @param idxHi largest index of range, inclusive
     * @param i the rank of the item to select in 0-based numbering.
     * @return the value of a's rank i item where i is 0-based
     */
    static double selectCLRS(double[] a, int idxLo, int idxHi, int i) {
        if (idxLo < 0 || idxLo >= a.length) {
            throw new IllegalArgumentException("idxLo is out of bounds");
        }
        if (idxHi < 0 || idxHi >= a.length) {
            throw new IllegalArgumentException("idxHi is out of bounds");
        }
        if (i < 0 || i > idxHi) {
            throw new IllegalArgumentException("i is out of bounds");
        }

        int n = (idxHi - idxLo + 1);
        if (n <= 5) {
            Arrays.sort(a, idxLo, idxHi + 1);
            return a[idxLo + i];
        }
        double result = select(a, idxLo, idxHi, i);
        //System.out.printf("Do these agree? %f, and a[i]=%f?\n", result, a[i]);
        //return result;
        return a[i];
    }

    /**
     find the value with rank i in array a with runtime complexity O(n) where n = a.length where i as a rank is 0-based.
     The worst case runtime complexity is O(n).

     <pre>
     references :
     Select, section 9.3
     "Introduction to Algorithms" by
     Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein

     Blum et al. 1983, "Time Bounds for Selection", STAN-CS-73-349
     </pre>
     * @param a
     * @param idxLo smallest index to of range
     * @param idxHi largest index of range, inclusive
     * @param i the rank of the item to select in 0-based numbering.
     * @return the median
     */
    private static double select(double[] a, int idxLo, int idxHi, int i) {

        final int _idxLo = idxLo;
        int _idxHi = idxHi;
        final int _i = i;

        //System.out.printf("* select idxLo=%d; idxHi=%d; i=%d\n", idxLo, idxHi, i);

        while ((idxHi - idxLo + 1) %5 != 0) {
            for (int j = idxLo + 1; j <= idxHi; ++j) {
                if (a[idxLo] > a[j] && idxLo != j) {
                    double tmp = a[idxLo];
                    a[idxLo] = a[j];
                    a[j] = tmp;
                }
            }
            if (i == 0 && (idxHi - idxLo) < 5) {
                //System.out.printf("return 0\n");
                return a[idxLo];
            }
            ++idxLo;
            --i;
        }

        int n = (idxHi - idxLo + 1);

        int g = n/5;
        int nRem = n - g*5;

        //System.out.printf("idxLo=%d; idxHi=%d; i=%d; g=%d; nRem=%d\n", idxLo, idxHi, i, g, nRem);

        if (idxHi < idxLo) {
            throw new IllegalArgumentException("   ERROR idxHi<idxLo\n");
        }

        for (int j = idxLo; j <= (idxLo + g - 1); ++j) {
            quickSort5(a, j, j + g*4, g);
            // sort each group (A[j], A[j+g], A[j+2g], A[j+3g, A[j+4g] in place
            // length is 5.
            // sort in place algoriths: quicksort, insertsort, bubblesort, selectionsort
            // all but the first are n^2.  quicksort is O(n*log(n)) on average.
        }
        //debugPrint(a, idxLo, idxHi, g);

        //Blum et al. 1973 1.b, pick recursively if n/5 > 1
        if (nRem==0 && g < 2 && (idxHi - idxLo) < 5) {
            //System.out.printf("   **** i=%d; idxLo=%d; idxHi=%d; g=%d, nRem=%d\n    a=%s\n",
            //        i, idxLo, idxHi, g, nRem,
            //        FormatArray.toString(a, "%.0f"));
            return a[idxLo + i];
        }

        //System.out.printf("NEXT select 1\n");
        //System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        // all group medians now lie in the middle fifth of A[idxLo:idcHi]

        int nextI = g - 1;//(int)Math.ceil(g/2);//((idxHi - idxLo) + 4) / 5;
        double x = select(a, idxLo + 2*g, idxLo + 3*g - 1, nextI);

        //if nAux == even number, we should consider both central numbers.  the other is at index (nAux/2) - 1.
        // or consider whether there is a way to append another number to aux (making the array 'odd' in length)
        // in a manner that finds the true ith rank number.

        //System.out.printf("aux pivot=%.0f\n",x);
        //System.out.printf("i=%d; idxLo=%d; idxHi=%d; g=%d; nRem=%d\n    a=%s\n", i, idxLo, idxHi, g, nRem,
        //        FormatArray.toString(a, "%.0f"));

        int q = partitionAround(a, idxLo, idxHi, x);
        // q is pivotIndex w.r.t 0
        // k is its rank
        int k = q - idxLo;

        //System.out.printf("pivotIdx = q = %d; pivot = %.0f; k = %d\n", q, x, k);
        //System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        double result;
        if (k==i) {
            //System.out.printf("NEXT select 2 (==q)\n");
            result = a[q];
        } else if (k>i) {
            //System.out.printf("NEXT select 3 (lower)\n");
            result = select(a, idxLo, q - 1, i);
        } else {
            //System.out.printf("NEXT select 4 (higher)\n");
            //return select(a, q + 1, idxHi, i-k);
            result = select(a, q + 1, idxHi, i - k - 1);
        }

        return result;
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
