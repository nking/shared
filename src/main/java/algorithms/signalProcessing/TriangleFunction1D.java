package algorithms.signalProcessing;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class TriangleFunction1D {
    

    /**
     * <pre>
     * An interpolation function called the triangle function
     * that uses a base 2 spacing. The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     *
     * The runtime complexity is O(N_points).
     *
     * c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     *
     * @param input
     * @param j level associated with input image. The output is calculated
     * using 2^j as spacing for interpolation points.
     * @return a sampling of input, interpolated over spacings 2^j.
     */
    public float[] calculateNextLevel(float[] input, int j) {

        return addOrSubtract(input, j, true);
    }
    
    /**
     * <pre>
     * An interpolation function called the triangle function
     * that uses a base 2 spacing to subtract to transformed levels. 
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     *
     * The runtime complexity is O(N_points).
     *
     * w_(j+1,k) = c_(j,k) − c_(j+1,k)
     *           = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     *
     * @param input
     * @param j level associated with input image. The output is calculated
     * using 2^j as spacing for interpolation points.
     * @return a sampling of input, interpolated over spacings 2^j.
     */
    public float[] subtractLevels(float[] input, int j) {

        return addOrSubtract(input, j, false);
    }

    private float[] addOrSubtract(float[] input, int j, boolean add) {

        int len = input.length;

        int s = 1 << j;

        float[] output = Arrays.copyOf(input, len);

        // use separability, that is 1D operation on columns, then rows
        
        for (int col = 0; col < len; ++col) {
            int x0 = col - s;
            int x2 = col + s;

            // choosing "continuity" for boundary corrections
            if (x0 < 0) {
                x0 = col;
            }
            if (x2 > (len - 1)) {
                x2 = col;
            }

            // add:
            //    c_(j + 1,k) = (1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) + (1/4)*c_(j,k+(2^j))
            // subtract:
            //    w_(j+1,k) = (-1/4)*c_(j,k-(2^j)) + (1/2)*c_(j,k) - (1/4)*c_(j,k+(2^j))
            float v0 = 0.25f * input[x0];
            float v1 = 0.5f * input[col];
            float v2 = 0.25f * input[x2];
            float vSum;
            if (add) {
                vSum = v0 + v1 + v2;
            } else {
                vSum = -1*v0 + v1 - v2;
            }

            output[col] = vSum;
        }

        return output;
    }

}
