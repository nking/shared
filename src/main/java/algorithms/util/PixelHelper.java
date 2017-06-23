package algorithms.util;

/**
 *
 * @author nichole
 */
public class PixelHelper {
    
    public int toPixelIndex(PairInt p, int width) {
        return (p.getY() * width) + p.getX();
    }
    
    public int toPixelIndex(int x, int y, int width) {
        return (y * width) + x;
    }
    
    public void toPixelCoords(int pixIdx, int width, 
        int[] outputXY) {
        outputXY[1] = pixIdx/width;
        outputXY[0] = pixIdx - (outputXY[1] * width);
    }
}
