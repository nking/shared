package thirdparty.ods;

import algorithms.util.ObjectSpaceEstimator;

/**
 * @author nichole
 @param <T>
 */
public class XFastTrieNode<T> extends BinaryTrieNode<T> {
    
    int prefix;

    @SuppressWarnings({"unchecked"})
    @Override
    public boolean equals(Object u) {
        
        if (!(u instanceof XFastTrieNode<?>)) {
            return false;
        }
        
        XFastTrieNode<T> uu = (XFastTrieNode<T>)u;
        
        if (this.prefix == uu.prefix) {
            return true;
        }
        
        return false;
    }

    @Override
    public int hashCode() {
        return prefix;
    }
   
}
