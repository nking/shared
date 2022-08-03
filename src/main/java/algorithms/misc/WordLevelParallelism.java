package algorithms.misc;

/**
 * word-level parallel operations.
 * useful for classes such as the FusionTree.
 * <pre>
 *  methods are implemented following lecture notes in
 *  http://web.stanford.edu/class/cs166/
 *  and
 *  ____ add other reference ---
 *  The MSB methods are ports of the c code at
 *  http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
 *  refactored here to use a variable tile size.
 * </pre>
 */
public class WordLevelParallelism {

    /**
     *
     * given an array of bitstringLength values, concatenate them and insert 0's on the high
     * end of each value.
     * e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *     tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     *     each tileBitLength bitstring, making a bitstring of length 16.
     * @param values array of bitstrings, each of length bitstringLength
     * @param bitstringLength the bitlength of each value in values
     * @return the bitarray holding the bitarray of replicated values with '0' separators.
     */
    public static long createTiledBitstring0(int[] values, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 63");
        }
        int n = values.length;
        if (n == 0) {
            return 0;
        }
        int i0 = 0;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc
        long kMult = 0;
        for (int i = 0; i < n; ++i) {
            kMult |= (values[i]*(1L<<i0));
            i0 += d;
        }
        // clear the gap bits
        i0 = bitstringLength;
        //kMask=(1<<15)|(1<<7) etc
        for (int i = 0; i < n; ++i) {
            kMult &= ~(1L << i0);
            i0 += d;
        }

        return kMult;
    }

    /**
     * create a bitmask array of set bits at the location of separators in the concatenation of
     * nTiles of length bitstringLength.
     * e.g. for nTiles=2 and bitstringLength=7, the resulting bitmask is 0b1000000010000000
     * which is 16 bits.
     * @param nTiles the number of tiles of bitstringLength for which this mask will be calculated.
     * @param bitstringLength the bit-length of each tile
     * @return the bitarray holding the bitarray of '1' separators for nTiles of length bitstringLength.
     * e.g. for nTiles=2 and bitstringLength=7, the resulting bitmask is 0b1000000010000000
     * which is 16 bits.
     */
    static long createTiledBitMask1(int nTiles, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 63");
        }
        if (nTiles == 0) {
            return 0;
        }
        int i1 = bitstringLength;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMask=(1<<15)|(1<<7) etc
        long kMask = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMask |= (1 << i1);
            i1 += d;
        }
        return kMask;
    }

    /**
     * given a bitstring called value which is much smaller than a machine word,
     * create a bitarray (word) with nTiles number of copies of value,
     * concatenated, with 1's in between them and on the high end.
     * e.g. for 7-bit value 0b1100111 and nTiles=2, the returned bitarray would be 0b1110011111100111.
     * @param value bitstring of length .lte. bitstringLength
     * @param nTiles the number of copies of value to set in the returned bitarray
     * @param bitstringLength the length of tiling before the 1's are concatenated as separators.
     * e.g. for a bitstringLength of 5 and nTiles=10, the resulting bitarray is length 10*(5+1)=60 bits
     * @return the bitarray holding the bitarray of replicated values with '1' separators.
     */
    public static long createTiledBitstring1(int value, int nTiles, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 63");
        }
        if (nTiles < 1 || nTiles > 63) {
            throw new IllegalArgumentException("nTiles must be greater than 0 and less than 63");
        }
        int i0 = 0;
        int i1 = bitstringLength;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc and kMask=(1<<15)|(1<<7) etc
        long kMult = 0;
        long kMask = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMult |= (1<<i0);
            kMask |= (1<<i1);
            i0 += d;
            i1 += d;
        }

        long kTiled = (value * kMult) | kMask;
        System.out.printf("value=%s, tiled=%s\n", Integer.toBinaryString(value), Long.toBinaryString(kTiled));
        return kTiled;
    }

    /**
    parallel compare of tiled1 to tiled2 and return a masked bit array whose set bits indicate which
    tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     </pre>
    @param tiled1 a bit array holding numbers of length tileBitLength (called tiles) separated by 0's.
    e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
    tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
    each tileBitLength bitstring, making a bitstring of length 16.
    @param tiled2 a bit array holding numbers of length tileBitLength separated by 0's.
     @param nTiles the number of tiles in the bitarray tiled1 or tiled2 (which should be the same number of tiles).
    @param tileBitLength the length of each tile in the bit arrays.
    @return a bit array of same size as tiled1 and tiled2 in which the bit of each
    tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
    in tiled2.
    */
    public static long parallelCompare00(long tiled1, long tiled2, int nTiles, int tileBitLength) {

        long mask1 = createTiledBitMask1(nTiles, tileBitLength);
        tiled1 |= mask1;

        return parallelCompare10(tiled1, tiled2, nTiles, tileBitLength, mask1);
    }

    /**
     parallel compare of tiled1 to tiled2 and return a masked bit array whose set bits indicate which
     tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     </pre>
     @param tiled1 a bit array holding numbers of length tileBitLength (called tiles) separated by 1's.
     e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     tiled1 is 0b1010010011100111, where 1's have been concatenated onto the high end of
     each tileBitLength bitstring, making a bitstring of length 16.
     @param tiled2 a bit array holding numbers of length tileBitLength separated by 0's.
     @param tileBitLength the length of each tile in the bit arrays.
     @param mask1 the 1's mask (same used in setting the gap bits in tiled1)
     @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     in tiled2.
     */
    public static long parallelCompare10(long tiled1, long tiled2, int nTiles,
                                         int tileBitLength, long mask1) {

        //3. Compute X – Y. The bit preceding xi – yi is 1 if xi ≥ yi and 0 otherwise.
        long diff = tiled1 - tiled2;

        long comparison = diff & mask1;

        System.out.printf("tiled1=%s\ntiled2=%s\ndiff=%s\ncomp=%s\n", Long.toBinaryString(tiled1),
                Long.toBinaryString(tiled2), Long.toBinaryString(diff), Long.toBinaryString(comparison));

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp

        //TODO: optimize this code later to reuse kMult
        int i0 = 0;
        final int d = tileBitLength + 1;
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc
        long kMult = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMult |= (1<<i0);
            i0 += d;
        }


        /*
        caveat: the first example is for a 64-bit string, and in java we're limited to unsigned 63 bit length.
        sumOf for 8 bit tiling of 7 bit bitstrings:
                   6         5         4         3         2         1
                3210987654321098765432109876543210987654321098765432109876543210
        value=0b1000000010000000100000000000000000000000000000000000000000000000#<== 8 tiles comparison
        kMult=0b0000000000000001000000010000000100000001000000010000000100000001#<== 7 set bits
        kMask=0b0000001110000000000000000000000000000000000000000000000000000000#<== 3 bit mask
        kShift = 64 - 8 - 1
        sum=(((value * kMult) & kMask) >> kShift) + (value >> 63)

        sumOf for 6 bit tiling of 5 bit bitstrings, 10 tiles (total was 64 bits, is now 60 bits):
                   6         5         4         3         2         1
                3210987654321098765432109876543210987654321098765432109876543210
        value=0b0000100000100000100000100000100000100000100000100000100000100000 #<== 10 tiles comparison
        kMult=0b0000000000000001000001000001000001000001000001000001000001000001 #<== 9 set bits
        kMask=0b0000000111100000000000000000000000000000000000000000000000000000 #<== 4 bit mask
                 # the number of tiles fits in 4 bits, so set the 4 bits at start of 2nd block
        kShift = 60 - 6 - 1
        sum=(((value * kMult) & kMask) >> kShift) + (value >> 59)
*/

        throw new UnsupportedOperationException("not yet implemente");
    }

}
