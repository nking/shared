package algorithms.signalProcessing;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class B3SplineFunction1D {
    
    /**
     * <pre>
     * An interpolation function for B-Spline, 3rd order.
     * The implementation follows pseudocode in
     * http://www.multiresolution.com/svbook.pdf
     * "Handbook of Astronomical Data Analysis" by 
     * Jean-Luc Starck and Fionn Murtagh
     * 
     * reduced to 1-D.
     * 
     * The runtime complexity is O(N_points).
     *
     * Handling boundaries:
     * "mirror" :      c(k + N) = c(N −k)
     * "periodicity" : (c(k + N) = c(N))
     * "continuity"  : (c(k + N) = c(k))
     * </pre>
     * @param input
     * @return
    */
    public float[] calculate(float[] input) {

        int len = input.length;

        float[] output = Arrays.copyOf(input, len);

        for (int col = 0; col < len; ++col) {

            output[col] = interpolate1DX(col, input);
        }

        return output;
    }

    /**
     * interpolate values around (idx) using a B3 spline.
     *
     * @param idx
     * @param data
     * @return
     */
    public float interpolate1DX(int idx, float[] data) {

        int len = data.length;

        /*
        (1/12)*(|x−2|^3 − 4*|x−1|^3 + 6*|x|^3 − 4*|x+1|^3 + |x+2|^3)
        1/16, 1/4, 3/8, 1/4, 1/16
        1/16*(1, 4, 6, 4, 1)
        */
        float vSum = 0;
        for (int dx = -2; dx <= 2; ++dx) {

            int xi = idx + dx;
            if ((xi < 0) || (xi > (len - 1))) {
                xi = idx;
            }

            float v = data[xi];

            switch(dx) {
                // -2 and +2
                case -1:
                case 1:
                    v *= 4.;
                    break;
                case 0:
                    v *= 6.;
                    break;
                // case -2 and +2 are factor 1
                default:
                    break;
            }

            vSum += v;
        }

        vSum /= 16.;

        return vSum;
    }

}

