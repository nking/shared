package algorithms;

import java.util.Arrays;
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
        // by this time, the array is largely sorted.
        System.out.printf("Do these agree? %f, and a[i]=%f?\n", result, a[i]);
        //return result;
        return a[i];
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

        int n = (idxHi - idxLo + 1);

        final int _idxLo = idxLo;
        int _idxHi = idxHi;
        final int _i = i;

        /*
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
        }*/

        int g = n/5;
        int nRem = n - g*5;

        System.out.printf("* select idxLo=%d; idxHi=%d; i=%d; g=%d; nRem=%d\n", idxLo, idxHi, i, g, nRem);

        if (idxHi < idxLo) {
            // idxHi + i or idxLo - 1 + i?
            System.out.printf("   ERROR idxHi<idxLo a=%s\n", FormatArray.toString(a, "%.0f"));
            //return a[idxLo + i];
        }

        //assert(idxLo + 5*g - 1 == idxHi);
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
            System.out.printf("   **** i=%d; idxLo=%d; idxHi=%d; g=%d, nRem=%d\n    a=%s\n",
                    i, idxLo, idxHi, g, nRem,
                    FormatArray.toString(a, "%.0f"));
            return a[idxLo + i];
        }

        System.out.printf("NEXT select 1\n");
        System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        // all group medians now lie in the middle fifth of A[idxLo:idcHi]
        // find the pivot x recursively as the median of the group medians.
        // this range holds the medians:
        // but if idxLo was offset above, this possibly needs to be lowered up to 2 "rows"

        // take median of x and the remainder.
        // how to include xRem into calcs for x without creating auxiliary arrays?
        // will use auxilliary array method first then when method is correct, can improve it
        // considering just passing in an an empty auxiliary array large enough to be reused for max g possible + nRem

        int nAux = g;
        double[] aux;
        if (nRem == 0) {
            aux = new double[nAux];
        } else {
            Arrays.sort(a, idxHi - nRem + 1, idxHi + 1);
            ++nAux;
            aux = new double[nAux];
            aux[nAux - 1] = a[ (idxHi-nRem+1) + (nRem/2)];
        }
        int _ii = 0;
        for (int ii = idxLo + 2*g; ii <= idxLo + 3*g - 1; ++ii) {
            aux[_ii++] = a[ii];
        }

        System.out.printf("check aux=%s\n", FormatArray.toString(aux, "%.0f"));

        //int nextI = (int)Math.ceil(g/2);//((idxHi - idxLo) + 4) / 5;
        //double x = select(a, idxLo + 2*g, idxLo + 3*g - 1, nextI);
        double x;
        if (nAux == 1) {
            x = aux[0];
        } else {
            x = select(aux, 0, nAux - 1, nAux/2);
        }

        //TODO: if nAux == even number, we should consider both central numbers.  the other is (nAux/2) - 1.
        // or consider whether there is a way to append another number (making the array odd in length)
        // in a manner that finds the true ith rank number.

        System.out.printf("aux pivot=%.0f\n",x);
        System.out.printf("i=%d; idxLo=%d; idxHi=%d; g=%d; nRem=%d\n    a=%s\n", i, idxLo, idxHi, g, nRem,
                FormatArray.toString(a, "%.0f"));

        // q is index of pivot x, 0-based
        int q = partitionAround(a, idxLo, idxHi, x);

        int k = q - idxLo + 1;
        System.out.printf("pivotIdx = q = %d; pivot = %.0f; k = %d\n", q, x, k);
        System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

        double result;
        if (k==i) {
            System.out.printf("NEXT select 2 (==q)\n");
            result = a[q];
        } else if (k>i) {
            System.out.printf("NEXT select 3 (lower)\n");
            result = select(a, idxLo, q - 1, i);
        } else {
            System.out.printf("NEXT select 4 (higher)\n");
            //return select(a, q + 1, idxHi, i-k);
            result = select(a, q + 1, idxHi, i-k);
        }

        if ((nAux & 1) == 0 && nAux > 1) {
            //Looks like this result is correct in context of my pivot algorithm

            System.out.printf("trying the other median of the even-sized aux array\n");
            double x2 = aux[(nAux/2) - 1];
            int q2 = partitionAround(a, idxLo, idxHi, x2);

            int k2 = q2 - idxLo + 1;

            System.out.printf("pivotIdx2 = q2 = %d; pivot2=%.0f; k2=\n", q2, x2, k2);
            System.out.printf("a=%s\n", FormatArray.toString(a, "%.0f"));

            double result2;
            if (k2==i) {
                System.out.printf("*NEXT select 2 (==q2)\n");
                result2 = a[q2];
            } else if (k2>i) {
                System.out.printf("*NEXT select 3 (lower)\n");
                result2 = select(a, idxLo, q2 - 1, i);
            } else {
                System.out.printf("*NEXT select 4 (higher)\n");
                //return select(a, q + 1, idxHi, i-k);
                result2 = select(a, q2 + 1, idxHi, i - k2);
            }
            System.out.printf("COMPARE result=%.0f; result2=%.0f\n", result, result2);
            result = result2;
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
