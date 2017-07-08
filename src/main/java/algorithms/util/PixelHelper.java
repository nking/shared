package algorithms.util;

import gnu.trove.set.TIntSet;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.iterator.TLongIterator;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.Set;
import java.util.HashSet;

/**
 *
 * @author nichole
 */
public class PixelHelper {
    
    public long toPixelIndex(PairInt p, int width) {
        return ((long)width * p.getY()) + p.getX();
    }
    
    public long toPixelIndex(int x, int y, int width) {
        return ((long)width * y) + x;
    }
    
    public void toPixelCoords(long pixIdx, int width, 
        int[] outputXY) {
        outputXY[1] = (int)(pixIdx/width);
        outputXY[0] = (int)(pixIdx - (outputXY[1] * width));
    }

    public TLongSet convert(Set<PairInt> points, int width) {
        TLongSet set = new TLongHashSet();
        for (PairInt p : points) {
            long pixIdx = toPixelIndex(p, width);
            set.add(pixIdx);
        }
        return set;
    }

    @SuppressWarnings({"unchecked"})
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
    
    @SuppressWarnings({"unchecked"})
    public Set<PairInt> convert(TLongSet pixIdxs, int width) {
        HashSet<PairInt> set = new HashSet<PairInt>();
        TLongIterator iter = pixIdxs.iterator();
        int[] xy = new int[2];
        while (iter.hasNext()) {
            long pixIdx = iter.next();
            toPixelCoords(pixIdx, width, xy);
            set.add(new PairInt(xy[0], xy[1]));
        }
        return set;
    }
}
