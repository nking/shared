package thirdparty.ods;

import algorithms.util.FNVHash;
import algorithms.util.ObjectSpaceEstimator;

/**
 * @author nichole
 @param <T> parameter type of node data
 */
public class XFastTrieNodeLong<T> extends BinaryTrieNode<T> {
    
    long prefix;

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
        return FNVHash.hash64a(new long[]{prefix});
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
