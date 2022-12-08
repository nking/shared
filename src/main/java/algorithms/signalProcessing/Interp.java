package algorithms.signalProcessing;

/**
 *
 * @author nichole
 */
public class Interp {
    
    /**
     *
     @param input
     @param binFactor
     @return
     */
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

                vSum += input[jj];

                count++;
            }

            if (count > 0) {
                vSum /= (float)count;
            }

            out[j] = vSum;
        }

        return out;
    }

    /**
     *
     @param input
     @param binFactor
     @return
     */
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
  
    /**
     *
     @param input
     @param outLength
     @param minValue
     @param maxValue
     @return
     */
    public float[] linearInterp(float[] input,
        int outLength, float minValue, float maxValue) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }

        final int w0 = input.length;
        final int w2 = outLength;

        float[] output = new float[w2];

        final float yFactor = (float)(w0 - 1)/(float)(w2 - 1);

        //System.out.println("w0=" + w0 + " w1=" + w2 + " f=" + yFactor);
        
        for (int i = 0; i < w2; ++i) {

            float i0 = (float)i * yFactor;
            int i0_0 = (int)i0;
            int i0_1 = (int)Math.ceil(i0);
            
            if (i0_1 >= input.length) {
                i0_1 = input.length - 1;
            }
             
            if (i0_1 == i0_0) {
                output[i] = input[i0_0];
                continue;
            }
                       
            float va = input[i0_0];
            float vb = input[i0_1];
            
            float d = (float)(i0_1 - i0_0);
            float fa = ((float)i0_1 - i0)/d;
            float fb = (i0 - (float)i0_0)/d;
                        
            float v = va * fa + vb * fb;
            
            output[i] = v;
        }

        return output;
    }
}
