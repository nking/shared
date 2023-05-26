package thirdparty.ods;

import algorithms.util.FNVHash;
import algorithms.util.ObjectSpaceEstimator;

import java.math.BigInteger;

/**
 * @author nichole
 @param <T>
 @param <T>
 */
public class XFastTrieNodeLong<T> extends BinaryTrieNode<T> {
    
    long prefix;

    private static final FNVHash fnv = new FNVHash();

    @SuppressWarnings({"unchecked"})
    @Override
    public boolean equals(Object u) {
        
        if (!(u instanceof XFastTrieNodeLong<?>)) {
            return false;
        }
        
        XFastTrieNodeLong<T> uu = (XFastTrieNodeLong<T>)u;
        
        if (this.prefix == uu.prefix) {
            return true;
        }
        
        return false;
    }

    @Override
    public int hashCode() {
        return fnv.hash(new long[]{prefix});
    }

    /**
     *
     @return
     */
    public static long estimateSizeOnHeap() {
        
        long total = BinaryTrieNode.estimateSizeOnHeap();
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNLongFields(1);
        
        total += est.estimateSizeOnHeap();
                   
        return total;
    }
}
