package algorithms.signalProcessing;

/**
 *
 * @author nichole
 */
public class Interp {
    
    public float[] bin(float[] input, int binFactor) {
        
        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = input.length;

        int w1 = w0/binFactor;

        float[] out = new float[w1];

        for (int j = 0; j < w1; j++) {

            float vSum = 0;
            int count = 0;

            for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {

                if ((jj < 0) || (jj > (w0 - 1))) {
                    continue;
                }

                float v = input[jj];

                vSum += v;
                count++;
            }

            if (count > 0) {
                float v = vSum/(float)count;
                vSum = Math.round(v);
            }

            out[j] = vSum;
        }

        return out;
    }

    public float[] unbin(float[] input, int binFactor) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        int w0 = input.length;

        float[] out = new float[binFactor* w0];

        int w1 = out.length;

        for (int j = 0; j < w0; j++) {

            float v = input[j];

            for (int jj = (j*binFactor); jj < ((j + 1)*binFactor); jj++) {
                out[jj] = v;
            }
            for (int jj = ((j + 1)*binFactor); jj < w1; jj++) {
                out[jj] = v;
            }
        }

        return out;
    }
  
    public float[] linearInterp(float[] input,
        int outLength, float minValue, float maxValue) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        final int w0 = input.length;
        final int w2 = outLength;

        if (w2 < w0) {
            throw new IllegalArgumentException("output dimensions cannot be"
                + " less than input dimensions for upsample");
        }
        
        float[] output = new float[w2];

        final float yFactor = (float)w2/(float)(w0 - 1);

        for (int i = 0; i < w2; ++i) {

            float i0 = i / yFactor;
            int i0_0 = (int)i0;
            int i0_1 = (int)Math.ceil(i0);
            
            if (i0_1 == i0_0) {
                output[i] = input[i0_0];
                continue;
            }
            
            float va = input[i0_0];
            float vb = input[i0_1];
            
            float v = ((i0_1 - i0)/(i0_1 - i0_0)) * va +
                ((i0 - i0_0)/(i0_1 - i0_0)) * vb;

            output[i] = v;
        }

        return output;
    }
}
