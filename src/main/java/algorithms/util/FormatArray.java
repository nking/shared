package algorithms.util;

/**
 *
 * @author nichole
 */
public class FormatArray {
    
    public static String toString(double[] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            sb.append(String.format(decimalFormat, a[i]));
            if (i < (a.length-1)) {
                sb.append(", ");
            }
        }
        return sb.toString();
    }
}
