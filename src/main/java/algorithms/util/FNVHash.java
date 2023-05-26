package algorithms.util;

import java.math.BigInteger;

/**
 * class for using FNV algorithm in making hashcodes for object identity.
 *
 <pre>
 references:
 http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
 http://www.isthe.com/chongo/src/fnv/fnv-5.0.3.tar.gz

 The FNV algorithm is under the public domain:

 CC0 - Public Domain
 FNV hash algorithms and source code been been put into the public domain via the following Creative Commons license:

 CC0 1.0 Universal (CC0 1.0) Public Domain Dedication
 No Copyright - CC0 - No Rights Reserved

 CC0 Public Domain

 The CC0 license means that the hash algorithms and source code has been dedicated to the public domain by waiving all of our rights to the work worldwide under copyright law, including all related and neighboring rights, to the extent allowed by law.

 You can copy, modify, distribute and perform the work, even for commercial purposes, all without asking permission.

 See the Creative Commons CC0 page for more details.
 </pre>
 */
public class FNVHash {

    protected final static BigInteger fnv321aInit = new BigInteger("2166136261");//0x811c9dc5

    protected final static BigInteger fnv32Prime = new BigInteger("16777619");//0x01000193

    protected final static BigInteger fnv64Init = new BigInteger("14695981039346656037");//0xcbf29ce484222325

    protected final static BigInteger fnv64Prime = new BigInteger("1099511628211");//0x100000001b3

    private int oldHash(int[] params) {
        int fnv321aInit = 0x811c9dc5;
        int fnv32Prime = 0x01000193;
        int hash = fnv321aInit;
        for (int p : params) {
            hash = hash ^ p;
            hash = hash * fnv32Prime;
        }
        return hash;
    }

    /**
     * create a hashcode for use in Object.hashCode() from the data in the object.
     * The method uses the FNV-1a algorithm which is in the public domain.
     * The algorithm is edited for use with java's signed integers produce a 31-bit number here that won't overflow an int.
     <pre>
     References:
     http://www.isthe.com/chongo/src/fnv/hash_32a.c
     http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
     Changing the FNV hash size - xor-folding
     </pre>
     * @param params parameters to use in making the hashcode
     * @return a non-cryptographic hash for use in object identity
     */
    public int hash(int[] params) {

        if (false) {
            // TODO:  finish testing.  reenable old hash meanwhile
            return oldHash(params);
        }

        /*
        FNV-1a alternate algorithm
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c

         and the fix for 31 bit for signed integers
         see http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
         see section Changing the FNV hash size - xor-folding
         */

        // BigInteger is needed as this variable overflows a signed integer
        BigInteger hash = fnv321aInit;

        for (int p : params) {
            hash = hash.xor(BigInteger.valueOf(p));
            hash = hash.multiply(fnv32Prime);
        }

        return to31Bit(hash);
    }

    /**
     * create a hashcode for use in Object.hashCode() from the data in the object.
     * The method uses the FNV-1a algorithm which is in the public domain.
     * The algorithm is edited for use with java's signed integers produce a 31-bit number here that won't overflow an int.
     <pre>
     References:
     http://www.isthe.com/chongo/src/fnv/hash_32a.c
     http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
     Changing the FNV hash size - xor-folding
     </pre>
     * @param params parameters to use in making the hashcode
     * @return a non-cryptographic hash for use in object identity
     */
    public int hash(float[] params) {

        /*
        FNV-1a alternate algorithm
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c

         and the fix for 31 bit for signed integers
         see http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
         see section Changing the FNV hash size - xor-folding
         */

        // BigInteger is needed as this variable overflows a signed integer
        BigInteger hash = fnv321aInit;

        for (float p : params) {
            hash = hash.xor(BigInteger.valueOf(Float.floatToIntBits(p)));
            hash = hash.multiply(fnv32Prime);
        }

        return to31Bit(hash);
    }

    /**
     * create a hashcode for use in Object.hashCode() from the data in the object.
     * The method uses the FNV-1a algorithm which is in the public domain.
     * The algorithm is edited for use with java's signed integers produce a 31-bit number here that won't overflow an int.
     <pre>
     References:
     http://www.isthe.com/chongo/src/fnv/hash_64a.c
     http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
     Changing the FNV hash size - xor-folding
     </pre>
     * @param params parameters to use in making the hashcode
     * @return a non-cryptographic hash for use in object identity
     */
    public int hash(long[] params) {

        /*
        FNV-1a alternate algorithm
         * hash = offset_basis
         * for each octet_of_data to be hashed
         *     hash = hash xor octet_of_data
         *     hash = hash * FNV_prime
         * return hash
         *
         * Public domain:  http://www.isthe.com/chongo/src/fnv/hash_32a.c

         and the fix for 31 bit for signed integers
         see http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
         see section Changing the FNV hash size - xor-folding
         */

        BigInteger hash = fnv64Init;

        for (long p : params) {
            hash = hash.xor(BigInteger.valueOf(p));
            hash = hash.multiply(fnv64Prime);
        }

        return to31Bit(hash);
    }

    protected int to31Bit(BigInteger hash) {
        // make it 31 bit: hash = (((hash>>31) ^ hash) & ((1<<31)-1))
        // make it n bit: hash = (((hash>>nBits) ^ hash) & ((1<<nBits)-1))
        BigInteger upper = hash.shiftRight(31);
        upper = upper.xor(hash);

        BigInteger mask = BigInteger.ONE.shiftLeft(31);
        mask = mask.subtract(BigInteger.ONE);

        hash = upper.and(mask);
        String s = hash.toString();
        return Integer.valueOf(s);
    }

    // for testing
    protected short _toShort(BigInteger hash) {

        // make it n bit: hash = (((hash>>nBits) ^ hash) & ((1<<nBits)-1))
        BigInteger upper = hash.shiftRight(15);
        upper = upper.xor(hash);

        BigInteger mask = BigInteger.ONE.shiftLeft(15);
        mask = mask.subtract(BigInteger.ONE);

        hash = upper.and(mask);
        String s = hash.toString();
        return Short.valueOf(s);
    }

    /**
     * for test use
     * @param params
     * @return
     */
    protected short _hash(short[] params) {

        // BigInteger is needed as this variable overflows a signed integer
        BigInteger hash = fnv321aInit;

        for (short p : params) {
            hash = hash.xor(BigInteger.valueOf(p));
            hash = hash.multiply(fnv32Prime);
        }

        return _toShort(hash);
    }
}
