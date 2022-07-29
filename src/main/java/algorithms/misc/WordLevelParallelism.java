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
     * given an array of bitstringLength bitstrings, concatenate them and insert 0's on the high
     * end of each value.
     * e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *     tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     *     each tileBitLength bitstring, making a bitstring of length 16.
     * @param values array of bitstrings, each of length bitstringLength
     * @param bitstringLength the bitlength of each value in values
     * @return the bitarray holding the bitarray of replicated values with '0' separators.
     */
    public static long createTiledBitstring(int[] values, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 62");
        }
        int n = values.length;
        if (n == 0) {
            return 0;
        }
        int i0 = 0;
        int i1 = bitstringLength;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc and kMask=kmask=(1<<15)|(1<<7) etc
        long kMult = 0;
        long kMask = 0;
        for (int i = 0; i < n; ++i) {
            kMult |= (values[i]*(1<<i0));
            kMask |= (1<<i1);
            i0 += d;
            i1 += d;
        }

        long kTiled = kMult | kMask;
        {
            StringBuilder sb = new StringBuilder();
            for (int i = n-1; i >= 0; --i) {
                sb.append("0").append(Integer.toBinaryString(values[i]));
            }
            System.out.printf("expected=%s, tiled=%s\n", sb.toString(), Long.toBinaryString(kTiled));
        }
    }

    /**
     * given a bitstring called value, create a bitarray with nTiles number of copies of value,
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
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 62");
        }
        if (nTiles < 1 || nTiles > 62) {
            throw new IllegalArgumentException("nTiles must be greater than 0 and less than 62");
        }
        int i0 = 0;
        int i1 = bitstringLength;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc and kMask=kmask=(1<<15)|(1<<7) etc
        long kMult = 0;
        long kMask = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMult |= (1<<i0);
            kMask |= (1<<i1);
            i0 += d;
            i1 += d;
        }

        long kTiled=(value * kMult) | kMask;
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
    @param tileBitLength the length of each tile in the bit arrays.
    @return a bit array of same size as tiled1 and tiled2 in which the bit of each
    tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
    in tiled2.
    */
    public static long parallelCompare(long tiled1, long tiled2, int tileBitLength) {

        /*
        1. Use a bitwise OR to place 1s between the xi’s.
        2. Use a bitwise AND to place 0s between the yi’s.
        3. Compute X – Y. The bit preceding xi – yi is 1 if xi ≥ yi and 0 otherwise.
         */
    }
}
