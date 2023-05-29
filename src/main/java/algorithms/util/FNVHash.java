package algorithms.util;

import java.math.BigInteger;

/**
  The FNV algorithm is a non-cyrptographic hash useful for object identity.
  FNV hashes are designed to be fast while maintaining a low collision rate.
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

 NOTE: a list of other non-cryptographic hashes is in https://en.wikipedia.org/wiki/List_of_hash_functions#Non-cryptographic_hash_functions

 */
public class FNVHash {

    /**
     * hash FNV32-a.
     * @param params
     * @return
     */
    public static int hash32a(int[] params) {
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
     * @return a non-cryptographic hashcode for use in object identity
     */
    public static int hash31a(int[] params) {

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
        BigInteger fnv321aInit = new BigInteger("2166136261");//0x811c9dc5

        BigInteger fnv32Prime = new BigInteger("16777619");//0x01000193

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
     * @return a non-cryptographic hashcode for use in object identity
     */
    public static int hash31a(float[] params) {

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
        BigInteger fnv321aInit = new BigInteger("2166136261");//0x811c9dc5

        BigInteger fnv32Prime = new BigInteger("16777619");//0x01000193

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
     <pre>
     References:
     http://www.isthe.com/chongo/src/fnv/hash_32a.c
     http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
     Changing the FNV hash size - xor-folding
     </pre>
     * @param params parameters to use in making the hashcode
     * @return a non-cryptographic hashcode for use in object identity
     */
    public static int hash32a(float[] params) {

        int fnv321aInit = 0x811c9dc5;
        int fnv32Prime = 0x01000193;
        int hash = fnv321aInit;
        for (float p : params) {
            hash = hash ^ Float.floatToIntBits(p);
            hash = hash * fnv32Prime;
        }
        return hash;
    }

    /**
     * create a hashcode for use in Object.hashCode() from the data in the object.
     * The method uses the FNV-1a algorithm which is in the public domain.
     <pre>
     References:
     http://www.isthe.com/chongo/src/fnv/hash_32a.c
     http://www.isthe.com/chongo/tech/comp/fnv/index.html#xor-fold
     Changing the FNV hash size - xor-folding
     </pre>
     * @param params parameters to use in making the hashcode
     * @return a non-cryptographic hashcode for use in object identity
     */
    public static int hash63a(long[] params) {

        BigInteger fnv64Init = new BigInteger("14695981039346656037");//0xcbf29ce484222325

        BigInteger fnv64Prime = new BigInteger("1099511628211");//0x100000001b3

        BigInteger hash = fnv64Init;

        for (long p : params) {
            hash = hash.xor(BigInteger.valueOf(p));
            hash = hash.multiply(fnv64Prime);
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
     * @return a non-cryptographic hashcode for use in object identity
     */
    public static int hash64a(long[] params) {

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

        long fnv64Init = 0xcbf29ce484222325L;

        long fnv64Prime = 0x100000001b3L;

        long hash = fnv64Init;

        for (long p : params) {
            hash = hash ^ p;
            hash = hash * fnv64Prime;
        }

        long upper = hash >> 31;
        upper = upper ^ hash;

        long mask = (1 << 31) - 1;

        hash = upper & mask;
        return (int)hash;
    }

    protected static int to31Bit(BigInteger hash) {
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
    protected static short _toShort(BigInteger hash) {

        // make it n bit: hash = (((hash>>nBits) ^ hash) & ((1<<nBits)-1))
        BigInteger upper = hash.shiftRight(15);
        upper = upper.xor(hash);

        BigInteger mask = BigInteger.ONE.shiftLeft(15);
        mask = mask.subtract(BigInteger.ONE);

        hash = upper.and(mask);
        String s = hash.toString();
        return Short.valueOf(s);
    }

    protected static short _oldhash(short[] params) {
        // use old hash then truncate to short
        int fnv321aInit = 0x811c9dc5;
        int fnv32Prime = 0x01000193;
        int hash = fnv321aInit;
        for (int p : params) {
            hash = hash ^ p;
            hash = hash * fnv32Prime;
        }
        // make it n bit: hash = (((hash>>nBits) ^ hash) & ((1<<nBits)-1))
        int h = (((hash>>15) ^ hash) & ((1<<15)-1));
        return (short)h;
    }

    /**
     * for test use
     * @param params
     * @return
     */
    protected static short _hash(short[] params) {
        BigInteger fnv321aInit = new BigInteger("2166136261");//0x811c9dc5

        BigInteger fnv32Prime = new BigInteger("16777619");//0x01000193

        // BigInteger is needed as this variable overflows a signed integer
        BigInteger hash = fnv321aInit;

        for (short p : params) {
            hash = hash.xor(BigInteger.valueOf(p));
            hash = hash.multiply(fnv32Prime);
        }

        return _toShort(hash);
    }
}
