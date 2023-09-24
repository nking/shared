package algorithms.util;

/**
 * a class with identity (x,y) that carries additional variable index.
 * only x, and y are used in equals and hash methods.
 *
 * @author nichole
 */
public class PairIntWithIndex extends PairInt {

    int pixIdx;

    public PairIntWithIndex(int xPoint, int yPoint, int thePixIndex) {
        super(xPoint, yPoint);
        pixIdx = thePixIndex;
    }
    
    public int getPixIndex() {
        return pixIdx;
    }

    @Override
    public boolean equals(Object obj) {

        if (!(obj instanceof PairInt)) {
            return false;
        }

        PairInt other = (PairInt) obj;

        return (getX() == other.getX()) && (getY() == other.getY());
    }

    @Override
    public int hashCode() {
        return FNVHash.hash32a(new int[]{getX(), getY()});
    }

    @Override
    public PairInt copy() {
         return new PairIntWithIndex(getX(), getY(), pixIdx);
    }

    @Override
    public String toString() {

        StringBuilder sb = new StringBuilder(super.toString());
        sb.append(" pixIdx=").append(Integer.toString(pixIdx));

        return sb.toString();
    }

}
