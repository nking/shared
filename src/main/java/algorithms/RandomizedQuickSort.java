package algorithms;

import java.util.Random;

public class RandomizedQuickSort {

    public static int partition(double[] a, int idxLo, int idxHi, Random rand) {
        int i = idxLo + rand.nextInt(idxHi - idxLo + 1);
        if (idxHi != i) {
            double swap = a[i];
            a[i] = a[idxHi];
            a[idxHi] = swap;
        }
        return partition(a, idxLo, idxHi);
    }

    public static void sort(double[] a, int idxLo, int idxHi, Random rand) {
        if (idxLo < idxHi) {
            int idxMid = partition(a, idxLo, idxHi, rand);
            sort(a, idxLo, idxMid - 1, rand);
            sort(a, idxMid + 1, idxHi, rand);
        }
    }

    private static int partition(double[] a, int idxLo, int idxHi) {
        double x = a[idxHi];
        int i = idxLo - 1;
        double swap;
        for (int k = idxLo; k < idxHi; ++k) {
            if (a[k] <= x) {
                ++i;
                if (i != k) {
                    swap = a[i];
                    a[i] = a[k];
                    a[k] = swap;
                }
            }
        }
        ++i;
        if (i != idxHi) {
            swap = a[idxHi];
            a[idxHi] = a[i];
            a[i] = swap;
        }
        return i;
    }
}
