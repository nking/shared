package algorithms.util;

import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import gnu.trove.iterator.TIntIterator;
import java.util.Set;
import java.util.HashSet;

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

    public TIntSet convert(Set<PairInt> points, int width) {
        TIntSet set = new TIntHashSet();
        for (PairInt p : points) {
            set.add(toPixelIndex(p, width));
        }
        return set;
    }

    public Set<PairInt> convert(TIntSet pixIdxs, int width) {
        HashSet<PairInt> set = new HashSet<PairInt>();
        TIntIterator iter = pixIdxs.iterator();
        int[] xy = new int[2];
        while (iter.hasNext()) {
            int pixIdx = iter.next();
            toPixelCoords(pixIdx, width, xy);
            set.add(new PairInt(xy[0], xy[1]));
        }
        return set;
    }
}
