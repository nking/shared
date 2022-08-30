package algorithms.util;

import algorithms.misc.Complex;
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
            if (i < (a.length - 1)) {
                sb.append(",");
            }
            sb.append(" ");
        }
        return sb.toString();
    }
    
    public static String toString(int[] a, String format) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            sb.append(String.format(format, a[i]));
            if (i < (a.length - 1)) {
                sb.append(",");
            }
            sb.append(" ");
        }
        return sb.toString();
    }
    
    public static String toString(boolean[] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            sb.append(String.format(decimalFormat, a[i]));
            if (i < (a.length - 1)) {
                sb.append(",");
            }
            sb.append(" ");
        }
        return sb.toString();
    }
    
    public static String toString(double[][] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(String.format(decimalFormat, a[i][j]));
                if (j < (a[i].length - 1)) {
                    sb.append(",");
                }
                sb.append(" ");
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
                if (j < (a.numColumns() - 1)) {
                    sb.append(",");
                }
                sb.append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    public static String toString(Complex[] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            sb.append("(").append(String.format(decimalFormat, a[i].re()))
            .append(", ").append(String.format(decimalFormat, a[i].im()))
            .append("j, abs=").append(String.format(decimalFormat, a[i].abs()))
            .append(", p=").append(String.format(decimalFormat, a[i].phase()))
            .append(")");
            if (i < a.length - 1) {
                sb.append(", ");
            }
        }
        return sb.toString();
    }

    public static String toString(int[][] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(String.format(decimalFormat, a[i][j]));
                if (j < (a[i].length - 1)) {
                    sb.append(",");
                }
                sb.append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
    
    public static String toString(float[][] a, String decimalFormat) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                sb.append(String.format(decimalFormat, a[i][j]));
                if (j < (a[i].length - 1)) {
                    sb.append(",");
                }
                sb.append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }
}
