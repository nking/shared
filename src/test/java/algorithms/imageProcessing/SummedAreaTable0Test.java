package algorithms.imageProcessing;

import java.security.SecureRandom;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class SummedAreaTable0Test extends TestCase {
    
    public SummedAreaTable0Test() {
    }
    
    public void testRandomExtract() throws Exception {
        
        SecureRandom random = SecureRandom.getInstance("SHA1PRNG");
        
        int w = 16;
        int h = 25;
        int nT = 100;//1<<29;
        
        double[][] img = new double[w][h];
        
        for (int i = 0; i < w; ++i) {
            for (int j = 0; j < h; ++j) {
                int r = random.nextInt(256);
                img[i][j] = r;
            }
        }
        
        double[] outSumAndN = new double[2];
        
        double diff;
        double eps = 1.e-15;
        
        SummedAreaTable0 sat = new SummedAreaTable0();
        double[][] sImg = sat.createAbsoluteSummedAreaTable(img);
        for (int i = 0; i < nT; ++i) {
            
            int x2 = random.nextInt(w - 1);
            x2++;
            int x1 = random.nextInt(x2);
            
            int y2 = random.nextInt(h - 1);
            y2++;
            int y1 = random.nextInt(y2);
            
            sat.extractWindowFromSummedAreaTable(sImg, x1, x2, y1, y2, outSumAndN);
            
            int sum = 0;
            for (int ii = x1; ii <= x2; ++ii) {
                for (int jj = y1; jj <= y2; ++jj) {
                    sum += img[ii][jj];
                }
            }
            
            diff = sum - outSumAndN[0];
            assertTrue(Math.abs(diff) < eps);
        }
    }
}
