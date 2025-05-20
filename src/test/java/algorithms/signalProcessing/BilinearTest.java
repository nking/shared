package algorithms.signalProcessing;

import algorithms.util.PairIntArray;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class BilinearTest extends TestCase {

    public BilinearTest(String testName) {
        super(testName);
    }

    public void testUpsample() {
        
        //    public static PairIntArray upsample(PairIntArray data, int upsamplingFactor) {
        PairIntArray data = new PairIntArray();
        data.add(0, 0);
        data.add(0, 6);
        data.add(6, 6);
        data.add(6, 0);

        PairIntArray expected = new PairIntArray();
        expected.add(0, 0);
        expected.add(0, 2);
        expected.add(0, 4);

        expected.add(0, 6);
        expected.add(2, 6);
        expected.add(4, 6);

        expected.add(6, 6);
        expected.add(6, 4);
        expected.add(6, 2);

        expected.add(6, 0);
        expected.add(4, 0);
        expected.add(2, 0);

        PairIntArray upsampled = Bilinear.upsample(data, 3);
        assertEquals(expected.getN(), upsampled.getN());

        for (int i = 0; i < upsampled.getN(); ++i) {
            assertEquals(expected.getX(i), upsampled.getX(i));
            assertEquals(expected.getY(i), upsampled.getY(i));
        }
    }
}
