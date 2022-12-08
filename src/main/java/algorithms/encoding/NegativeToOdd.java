package algorithms.encoding;

/**
 *
 * @author nichole
 */
public class NegativeToOdd {
    
    /**
     *
     @param a
     */
    public static void negativeToOddEncode(int[] a) {
        for (int i = 0; i < a.length; ++i) {
            a[i] *= 2;
            if (a[i] < 0) {
                a[i] *= -1;
                a[i] += 1;
            }
        }
    }

    /**
     *
     @param a
     */
    public static void negativeToOddDecode(int[] a) {
        for (int i = 0; i < a.length; ++i) {
            if ((a[i] & 1) != 0) {
                a[i] -= 1;
                a[i] *= -1;
            }
            a[i] /= 2;
        }
    }
}
