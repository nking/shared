package thirdparty.ods;

import algorithms.util.ObjectSpaceEstimator;

/**
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
        return fnvHashCode(prefix);
    }
     
    protected final static int fnv321aInit = 0x811c9dc5;
    protected final static int fnv32Prime = 0x01000193;

    protected int fnvHashCode(long p) {

        /*
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c
         */

        int hash = 0;

        int sum = fnv321aInit;

        int mask = Integer.MAX_VALUE;
        int shift = 31;
        int i0 = (int)(p & mask);
        int i1 = (int)((p >> shift) & mask);
        
        // xor the bottom with the current octet.
        sum ^= i0;

        // multiply by the 32 bit FNV magic prime mod 2^32
        sum *= fnv32Prime;
        
        sum ^= i1;
        
        sum *= fnv32Prime;
        
        hash = sum;

        return hash;
    }

    
    public static long estimateSizeOnHeap() {
        
        long total = BinaryTrieNode.estimateSizeOnHeap();
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNLongFields(1);
        
        total += est.estimateSizeOnHeap();
                   
        return total;
    }
}
