package algorithms.util;

import no.uib.cipr.matrix.DenseMatrix;

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
    
    public static String toString(double[][] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(String.format(decimalFormat, a[i][j]));
                if (j < (a.length-1)) {
                    sb.append(", ");
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    public static String toString(DenseMatrix a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.numRows(); ++i) {
            for (int j = 0; j < a.numColumns(); ++j) {
                sb.append(String.format(decimalFormat, a.get(i, j)));
                if (j < (a.numColumns()-1)) {
                    sb.append(", ");
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
}
