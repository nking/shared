package algorithms.signalProcessing;

import algorithms.util.PairIntArray;

public class Bilinear {

    public static PairIntArray upsample(PairIntArray data, int upsamplingFactor) {

        PairIntArray out = new PairIntArray();
        float factor = (float)upsamplingFactor;
        out.add(data.getX(0), data.getY(0));

        for (int i = 0; i < data.getN(); i++) {
            int x0 = data.getX(i);
            int y0 = data.getY(i);
            int x1, y1;
            if (i + 1 == data.getN()) {
                x1 = data.getX(0);
                y1 = data.getY(0);
            } else {
                x1 = data.getX(i + 1);
                y1 = data.getY(i + 1);
            }
            for (int j = 1; j < upsamplingFactor; j++) {
                float t = (float) j / factor;
                float x = x0 + t * (x1 - x0);
                float y = y0 + t * (y1 - y0);
                out.add(Math.round(x), Math.round(y));
            }
            if (i < data.getN() - 1) {
                out.add(x1, y1);
            }
        }
        return out;
    }
}
