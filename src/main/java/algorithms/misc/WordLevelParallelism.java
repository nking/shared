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

    //TODO: optimize code

    /**
     * given a bitarray packed full of tiles separated by flags (= blocks),
     * returns the index of the highest nonzero block.
     *
     * <pre>
     *     reference http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     *     then edited here for variable block size and number of tiles packed into tiled.
     * </pre>
     *
     * @param tiled         a bitarray of concatenated bitstrings of length tileBitLength separated by flag bits.
     *                      the portion of tiled read is the first nTiles * (tileBitLength + 1) bits.
     * @param nTiles        the number of tiles packed into the bitarray tiled.
     * @param tileBitLength the size of a tile before a gap is appended to it.  the block size is tileBitlength + 1.
     * @return the index of the highest nonzero block in tiled
     */
    public static long highestOneBitIn(long tiled, int nTiles, int tileBitLength) {

        int bSz = tileBitLength + 1;

        /* Step 1: Identify the index of the highest block with a 1 bit in it. */
        long highBlockIndex = highestBlockSetIn(tiled, nTiles, tileBitLength);

        /* Step 2: Identify the highest bit within that block. To do so, we're going
         * to shift that block down to the proper position and mask out the other
         * bits.
         */
        //long highBlock = tiled >> (highBlockIndex * 8);
        long highBlock = tiled >> (highBlockIndex * bSz);

        return highBlockIndex * bSz + highestBitSetIn(highBlock, bSz);
    }

    /**
     * Given a 64-bit integer, returns the index of the block within that integer
     * that contains a 1 bit, where the index is zero-based numbering.
     * @param tiled the bitarray of tiled integers.
     * @param nTiles the number of tiles embedded in tiled.
     * @param tileBitLength the length of each tile in tiled, not counting the surrounding single flag bits.
     *                      the block size is tileBitLength + 1.
     * @return
     */
    static long highestBlockSetIn(long tiled, int nTiles, int tileBitLength) {

        int bSz = tileBitLength + 1;

        long usedBlocksIn = usedBlocksIn(tiled, bSz);

        // the block number in bits, e.g. 6th block is 0b1000000
        long sketch = sketch(usedBlocksIn, nTiles, tileBitLength);

        return highestBitSetIn(sketch, tileBitLength);
    }

    /**
     * given an array of bitstringLength values, concatenate them and insert 0's on the high
     * end of each value.
     * e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     * tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     * each tileBitLength bitstring, making a bitstring of length 16.
     * <p>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (values.length * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent values.length is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from values array in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param values          array of bitstrings, each of length bitstringLength
     * @param bitstringLength the bitlength of each value in values.  the tile for each will be bitstringLength + 1
     *                        bits long.  the total tiled result will be values.length * (bitstringLength + 1) bits.
     * @return the bitarray holding the bitarray of replicated values with '0' separators.
     */
    public static long createTiledBitstring0(int[] values, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 62");
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
            kMult |= (values[i] * (1L << i0));
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
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param nTiles          the number of tiles of bitstringLength for which this mask will be calculated.
     * @param bitstringLength the bit-length of each tile
     * @return the bitarray holding the bitarray of '1' separators for nTiles of length bitstringLength.
     * e.g. for nTiles=2 and bitstringLength=7, the resulting bitmask is 0b1000000010000000
     * which is 16 bits.
     */
    static long createTiledBitMask1(int nTiles, int bitstringLength) {
        if (bitstringLength < 1 || bitstringLength > 61) {
            throw new IllegalArgumentException("bitstringLength must be greater than 0 and less than 62");
        }
        if (nTiles == 0) {
            return 0;
        }
        int i1 = bitstringLength;
        final int d = bitstringLength + 1;
        // e.g. for bitstringLength=7, kMask=(1<<15)|(1<<7) etc
        long kMask = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMask |= (1L << i1);
            i1 += d;
        }
        return kMask;
    }

    /**
     * given a bitstring called value which is much smaller than a machine word,
     * create a bitarray (word) with nTiles number of copies of value,
     * concatenated, with 1's in between them and on the high end.
     * e.g. for 7-bit value 0b1100111 and nTiles=2, the returned bitarray would be 0b1110011111100111.
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bitstringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value           bitstring of length .lte. bitstringLength
     * @param nTiles          the number of copies of value to set in the returned bitarray
     * @param bitstringLength the length of tiling before the 1's are concatenated as separators.
     *                        e.g. for a bitstringLength of 5 and nTiles=10, the resulting bitarray is length 10*(5+1)=60 bits
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
        // e.g. for bitstringLength=7, kMult=(1<<8)|(1<<0) etc and kMask=(1<<15)|(1<<7) etc
        long kMult = 0;
        long kMask = 0;
        for (int i = 0; i < nTiles; ++i) {
            kMult |= (1L << i0);
            kMask |= (1L << i1);
            i0 += d;
            i1 += d;
        }

        return (value * kMult) | kMask;
    }

    /**
     * parallel compare of tiled1 to tiled2 and return a masked bit array whose set bits indicate which
     * tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bitstringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param tiled1        a bit array holding numbers of length tileBitLength (called tiles) separated by 0's.
     *                      e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *                      tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     *                      each tileBitLength bitstring, making a bitstring of length 16.
     * @param tiled2        a bit array holding numbers of length tileBitLength separated by 0's.
     * @param nTiles        the number of tiles in the bitarray tiled1 or tiled2 (which should be the same number of tiles).
     * @param tileBitLength the length of each tile in the bit arrays.
     * @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     * tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     * in tiled2.
     */
    public static long parallelCompare00(long tiled1, long tiled2, int nTiles, int tileBitLength) {

        long mask1 = createTiledBitMask1(nTiles, tileBitLength);
        tiled1 |= mask1;

        return parallelCompare10(tiled1, tiled2, nTiles, tileBitLength, mask1);
    }

    /**
     * parallel compare of tiled1 to tiled2 which both have block sizes of 8 and have flag bits
     * of '0' separating the embedded 7-bit tiles.
     * returns a masked bit array whose set bits indicate which
     * tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bitstringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param tiled1        a bit array holding numbers of length tileBitLength (called tiles) separated by 0's.
     *                      e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *                      tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     *                      each tileBitLength bitstring, making a bitstring of length 16.
     * @param tiled2        a bit array holding numbers of length tileBitLength separated by 0's.
     * @param nTiles        the number of tiles in the bitarray tiled1 or tiled2 (which should be the same number of tiles).
     * @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     * tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     * in tiled2.
     */
    public static long parallelCompare008(long tiled1, long tiled2, int nTiles) {
        //                6         5         4         3         2         1
        //              210987654321098765432109876543210987654321098765432109876543210
        long kMask1 = 0b000000010000000100000001000000010000000100000001000000010000000L;
        // e.g. for bitstringLength=7, block size=8, kMask=(1<<15)|(1<<7) etc

        tiled1 |= kMask1;

        return parallelCompare10(tiled1, tiled2, nTiles, 7, kMask1);
    }

    /**
     * parallel compare of tiled1 to tiled2 and return a masked bit array whose set bits indicate which
     * tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param tiled1        a bit array holding numbers of length tileBitLength (called tiles) separated by 1's.
     *                      e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *                      tiled1 is 0b1010010011100111, where 1's have been concatenated onto the high end of
     *                      each tileBitLength bitstring, making a bitstring of length 16.
     * @param tiled2        a bit array holding numbers of length tileBitLength separated by 0's.
     * @param tileBitLength the length of each tile in the bit arrays.  the block size is tileBitLength + 1
     *                      because it includes the gap bit between tiles.
     * @param mask1         the 1's mask (same used in setting the gap bits in tiled1)
     * @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     * tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     * in tiled2.
     */
    public static long parallelCompare10(long tiled1, long tiled2, int nTiles, int tileBitLength, long mask1) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        if (nTiles < 1) {
            throw new IllegalArgumentException("nTiles must be > 0");
        }

        final int bSz = tileBitLength + 1;

        final int nMaskBits = (int) Math.ceil(Math.log(nTiles) / Math.log(2));

        // by default the number of blacks used to hold the number nTiles is 1.
        // if number of bits in nTiles > bSz, nBExtra is the number of blacks to add to the default '1 reserved block'
        // to hold the number of bits in nTiles.
        // e.g. nTiles=7 can be held in a block of size 3 bits. if bSz=2, need 1 extra block to hold nTiles.
        //      nTiles=16 can be held in a block of size 5 bits. if bSz=2, need 2 extra blocks to hold nTiles.
        int nBExtra = 0;

        // assert that there is enough space to hold the bits to represent nTiles
        if (nMaskBits > bSz) {
            // how many blocks needed to store nMaskBits?  then subtract 1 which is already reserved for it.
            nBExtra = ((int) Math.ceil((double) nMaskBits / bSz)) - 1;
            int tiledLength = nTiles * bSz;
            /*if ((tiledLength + nBExtra * bSz) > 62) {
                throw new IllegalArgumentException(String.format("nTiles needs %d blocks of size %d bits above the " +
                                "total tiled bit length =%d.\n  That total %d must fit within 62 bits.",
                        nBExtra, bSz, tiledLength, (tiledLength + nBExtra * bSz)));
            }*/
        }

        //3. Compute X – Y. The bit preceding xi – yi is 1 if xi ≥ yi and 0 otherwise.
        long diff = tiled1 - tiled2;

        long comparison = diff & mask1;

        //System.out.printf("tiled1=%30s\ntiled2=%30s\ndiff=%32s\ncomp=%32s\n", Long.toBinaryString(tiled1),
        //        Long.toBinaryString(tiled2), Long.toBinaryString(diff), Long.toBinaryString(comparison));

        return parallelSum(comparison, nTiles, tileBitLength);
    }

    /**
     * parallel compare of tiled1 to tiled2 both of which have block size 8 and
     * embedded tiles of size 7-bits.
     * this method returns a masked bit array whose set bits indicate which
     * tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param tiled1        a bit array holding numbers of length tileBitLength (called tiles) separated by 1's.
     *                      e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *                      tiled1 is 0b1010010011100111, where 1's have been concatenated onto the high end of
     *                      each tileBitLength bitstring, making a bitstring of length 16.
     * @param tiled2        a bit array holding numbers of length tileBitLength separated by 0's.
     * @param mask1         the 1's mask (same used in setting the gap bits in tiled1)
     * @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     * tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     * in tiled2.
     */
    static long parallelCompare108(long tiled1, long tiled2, int nTiles, long mask1) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        if (nTiles < 1) {
            throw new IllegalArgumentException("nTiles must be > 0");
        }

        final int bSz = 8;

        final int nMaskBits = (int) Math.ceil(Math.log(nTiles) / Math.log(2));

        // by default the number of blacks used to hold the number nTiles is 1.
        // if number of bits in nTiles > bSz, nBExtra is the number of blacks to add to the default '1 reserved block'
        // to hold the number of bits in nTiles.
        // e.g. nTiles=7 can be held in a block of size 3 bits. if bSz=2, need 1 extra block to hold nTiles.
        //      nTiles=16 can be held in a block of size 5 bits. if bSz=2, need 2 extra blocks to hold nTiles.
        int nBExtra = 0;

        // assert that there is enough space to hold the bits to represent nTiles
        if (nMaskBits > bSz) {
            // how many blocks needed to store nMaskBits?  then subtract 1 which is already reserved for it.
            nBExtra = ((int) Math.ceil((double) nMaskBits / bSz)) - 1;
            int tiledLength = nTiles * bSz;
            /*if ((tiledLength + nBExtra * bSz) > 62) {
                throw new IllegalArgumentException(String.format("nTiles needs %d blocks of size %d bits above the " +
                                "total tiled bit length =%d.\n  That total %d must fit within 62 bits.",
                        nBExtra, bSz, tiledLength, (tiledLength + nBExtra * bSz)));
            }*/
        }

        //3. Compute X – Y. The bit preceding xi – yi is 1 if xi ≥ yi and 0 otherwise.
        long diff = tiled1 - tiled2;

        long comparison = diff & mask1;

        //System.out.printf("tiled1=%30s\ntiled2=%30s\ndiff=%32s\ncomp=%32s\n", Long.toBinaryString(tiled1),
        //        Long.toBinaryString(tiled2), Long.toBinaryString(diff), Long.toBinaryString(comparison));

        return parallelSum(comparison, nTiles, 7);
    }

    /**
     * sum the set bits of bitstring comparison.  the flags that may have set bits are the MSB if each block.
     *
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param comparison    a bit array with flags at the MSB of each block.  The flags that are set bits
     *                      are summed in this method.
     * @param tileBitLength the length of each tile in the bit arrays.  the block size is tileBitLength + 1
     *                      because it includes the gap bit between tiles.
     * @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     * tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     * in tiled2.
     */
    static long parallelSum(long comparison, int nTiles, int tileBitLength) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        final int bSz = tileBitLength + 1;

        final int nMaskBits = (int) Math.ceil(Math.log(nTiles) / Math.log(2));

        // by default the number of blacks used to hold the number nTiles is 1.
        // if number of bits in nTiles > bSz, nBExtra is the number of blacks to add to the default '1 reserved block'
        // to hold the number of bits in nTiles.
        // e.g. nTiles=7 can be held in a block of size 3 bits. if bSz=2, need 1 extra block to hold nTiles.
        //      nTiles=16 can be held in a block of size 5 bits. if bSz=2, need 2 extra blocks to hold nTiles.
        int nBExtra = 0;

        // assert that there is enough space to hold the bits to represent nTiles
        if (nMaskBits > bSz) {
            // how many blocks needed to store nMaskBits?  then subtract 1 which is already reserved for it.
            nBExtra = ((int) Math.ceil((double) nMaskBits / bSz)) - 1;
            int tiledLength = nTiles * bSz;
            /*if ((tiledLength + nBExtra * bSz) > 62) {
                throw new IllegalArgumentException(String.format("nTiles needs %d blocks of size %d bits above the " +
                                "total tiled bit length =%d.\n  That total %d must fit within 62 bits.",
                        nBExtra, bSz, tiledLength, (tiledLength + nBExtra * bSz)));
            }*/
        }

        final long kMult;
        final long kMask;
        final long kShift;
        long sumBits;

        switch(bSz) {
            case 8:
                // for nTiles=7, need 3-bit mask
                //         6         5         4         3         2         1
                //       210987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000000001000000010000000100000001000000010000000100000001L;
                kMask = 0b00001110000000000000000000000000000000000000000000000000000000L;
                //              A0000000B0000000C0000000D0000000E0000000F0000000G0000000L
                kShift = 55;
                return ((comparison * kMult) & kMask) >> kShift;
            case 7:
                // for nTiles = 8, need 4 bits in mask
                // kMult=(1<<0) | (1<<7) | (1<<14) | (1<<21) | (1<<28) | (1<<35) | (1<<42) | (1<<49)
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000000010000001000000100000010000001000000100000010000001L;
                kMask = 0b00011110000000000000000000000000000000000000000000000000000000L;
                //              10000001000000100000010000001000000100000010000001000000;
                kShift = 55;
                return ((comparison * kMult) & kMask) >> kShift;
            case 6:
                // for nTiles = 10, need 4 bits in mask
                // kMult=(1<<0) | (1<<6) | (1<<12) | (1<<18) | (1<<24) | (1<<30) | (1<<36) | (1<<42) | (1<<48) | (1<<54)
                //           6         5         4         3         2         1
                //          10987654321098765432109876543210987654321098765432109876543210
                //kMult = 0b00000001000001000001000001000001000001000001000001000001000001L;
                //kMask =    0b11110000000000000000000000000000000000000000000000000000000L;
                //            A00000B00000C00000D00000E00000F00000G00000H00000I00000J00000L;
                // not enough space in 62 bits for the bit mask, so handle nTiles = 9,
                // then set a higher bit if block 9 high bit is set
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000000001000001000001000001000001000001000001000001000001L;
                kMask =      0b111100000000000000000000000000000000000000000000000000000L;
                //                  A00000B00000C00000D00000E00000F00000G00000H00000I00000L;
                kShift = 53;
                sumBits = ((comparison * kMult) & kMask) >> kShift;
                sumBits += (comparison >> 59);
                return sumBits;
            case 5:
                // for nTiles = 12, need 4 bits in mask
                // not enough space for mask, so will calculate for nTiles=11 and set the high bit for last
                // kMult=(1<<0) | (1<<5) | (1<<10) | (1<<15) | (1<<20) | (1<<25) | (1<<30) | (1<<35) | (1<<40) | (1<<45) | (1<<50) | (1<<55)
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000000100001000010000100001000010000100001000010000100001L;
                kMask =     0b1111000000000000000000000000000000000000000000000000000000L;
                //           A0000B0000C0000D0000E0000F0000G0000H0000I0000J0000K0000L0000L;
                kShift = 54;
                sumBits = ((comparison * kMult) & kMask) >> kShift;
                sumBits += (comparison >> 59);
                return sumBits;
            case 4:
                // for nTiles = 15, need 4 bits in mask
                // not enough space for nTiles=15 and 4 bits of mask,
                // so will calculate sum for nTiles=14 and set last high bit
                // kMult=(1<<0) | (1<<4) | (1<<8) | (1<<12) | (1<<16) | (1<<20) | (1<<24) | (1<<28) | (1<<32) | (1<<36)
                // | (1<<40) | (1<<44) | (1<<48) | (1<<52)
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000010001000100010001000100010001000100010001000100010001L;
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMask =    0b11110000000000000000000000000000000000000000000000000000000L;
                //           A0000B0000C0000D0000E0000F0000G0000H0000I0000J0000K0000L0000L;
                kShift = 55;
                sumBits = ((comparison * kMult) & kMask) >> kShift;
                sumBits += (comparison >> 59);
                return sumBits;
            case 3:
                // for nTiles = 20, need 5 bits in mask.
                // not enough space in 62 bits to hold the 5 bit mask for nTiles=20,
                // so will handle nTiles=18 and set high bits for blocks 18 and 19
                // kMult=(1<<0) | (1<<3)  | (1<<6)  | (1<<9)  | (1<<12) | (1<<15) | (1<<18) | (1<<21) | (1<<24) | (1<<27)
                //    | (1<<30) | (1<<33) | (1<<36) | (1<<39) | (1<<42) | (1<<45) | (1<<48) | (1<<51)
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000001001001001001001001001001001001001001001001001001001L;
                kMask =     0b11111000000000000000000000000000000000000000000000000000000L;
                kShift = 54;
                sumBits = ((comparison * kMult) & kMask) >> kShift;
                comparison >>= 56;
                sumBits += (comparison & 0b1) + ((comparison >> 3) & 0b1);
                return sumBits;
            case 2:
                editing here for a more efficient way
                // for nTiles = 31, need 5 bits in mask.
                // not enough space in 62 bits to hold the 5 bit mask for nTiles=31,
                // so will handle nTiles=29
                // then set a higher bit if block 29 high bit is set and same for block 30
                // kMult=(1<<0) | (1<<2)  | (1<<4)  | etc
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000101010101010101010101010101010101010101010101010101010101L;
                //              9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1 0 9 8 7 6 5 4 3 2 1
                //kMask = 0b;
                //            1_2_3_4_5_6_7_8_9_0_1_2_3_4_5_6_7_8_9_0_1_2_3_4_5_6_7_8_9_0_1_L
                //kShift = 38;
                sumBits = ((comparison * kMult) & kMask) >> kShift;
                break;
            case 1:
                // for nTiles = 62, need 6 bits in mask.
                //         6         5         4         3         2         1
                //        10987654321098765432109876543210987654321098765432109876543210
                kMult = 0b00000000000000000000000000000000000000000000000000000000000001L;
                kMask = 0b0;
                editing  should be a more efficient way to impl this one
                break;
            default:
                throw new UnsupportedOperationException("not yet implemented");
        }
    }

    /**
     * sum the set bits of the comparison bitstring which has a block size of 8 and embedded tiles of
     * size 7-bits.  this method sums the set bits at the high end of each block and returns the result.
     *
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     * operations.
     * Let block size = (bistringLength + 1).
     * The unsigned long restricts the total bit length of the tiled result of this method to 62 bits,
     * and so (nTiles * block) must be less than or equal to 62.
     * Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     * representation of the number of tiles into the highest blocks of a mask that is the same size as the
     * total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     * block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     * from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     * and the 62 bit limit.
     *
     * @param comparison    a bit array with flags at the MSB of each block.  The flags that are set bits
     *                      are summed in this method.
     * @return the sum of the set bits of the MSB of each 8-bit block.
     */
    static long parallelSum8(long comparison, int nTiles) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        final int bSz = 8;

        // e.g. for bSz=8, kMult=(1<<8)|(1<<0) etc
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        //long kMult = 0b000000100000001000000010000000100000001000000010000000100000001L;
        long kMult = 0b00000000000001000000010000000100000001000000010000000100000001L;
        long kMask = 0b00001110000000000000000000000000000000000000000000000000000000L;

        int kShift = 64-8-1;
        int kShift2 = 63;

        long s1 = (((comparison * kMult) & kMask) >> kShift);
        long s2 = (comparison >> kShift2);

        long sum = s1 + s2;

        return sum;
    }

    /**
     * given a bitarray packed full of tiles separated by flags, extract and return the flags.
     * e.g. if tiled were A0000000B0000000C0000000D0000000, this method would return ABCD.
     <pre>
     references:
     lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
          and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp

     https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge

     </pre>
     *
     * @param tiled         a bitarray of concatenated bitstrings of length tileBitLength separated by flag bits.
     *                      the portion of tiled read is the first nTiles * (tileBitLength + 1) bits.
     * @param nTiles        the number of tiles packed into the bitarray tiled.
     * @param tileBitLength the size of a tile before a gap is appended to it.  the block size is tileBitlength + 1.
     * @return
     */
    public static long sketch(long tiled, int nTiles, int tileBitLength) {

        switch(tileBitLength + 1) {
            case 8:
                return sketch8(tiled, nTiles);
            case 7:
                return sketch7(tiled, nTiles);
            case 6:
                return sketch6(tiled, nTiles);
            case 5:
                return sketch5(tiled, nTiles);
            case 4:
                return sketch4(tiled, nTiles);
            case 3:
                return sketch3(tiled, nTiles);
            case 2:
                return sketch2(tiled);
            case 1:
                return tiled;
            default:
                throw new UnsupportedOperationException("not yet implemented");
        }
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 8 bits
     * and embedded tile size of 7 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A0000000B0000000C0000000D0000000, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 7 bits
     * and embedded tile size of 6 bits
     * @return
     */
    public static long sketch8(long tiled, int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 7) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 7 for block size 8");
        }

        // e.g. for bitstringLength=7 blockSize=8, kMult=(1<<0) | (1<<7) | (1<<14) etc
        //                6         5         4         3         2         1
        //               10987654321098765432109876543210987654321098765432109876543210
        //                     10000001000000100000010000001000000100000010000001000000
        long kMult   = 0b00000000000000000001000000100000010000001000000100000010000001L;
        long kMask   = 0b11111110000000000000000000000000000000000000000000000000L;
        int kShift   = 49;

        /*System.out.printf("\nkMask=\n%63s\n" +
                "kMult=\n%63s\n" +
                "kShift=%d\n" +
                "(tiled * kMult) & kMask=\n%63s\n" +
                "(((tiled * kMult) & kMask) >> kShift)=\n%63s\n",
                Long.toBinaryString(kMask), Long.toBinaryString(kMult), kShift,
                Long.toBinaryString((tiled * kMult) & kMask),
                Long.toBinaryString(((tiled * kMult) & kMask) >> kShift));*/

        long sketch = ((tiled * kMult) & kMask) >> kShift;

        return sketch;

        /*
                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                                         000000010000000100000001000000010000000100000001000000010000000
                                                              0b6       5       4       3       2       1       0       L
                                                       0b6       5       4       3       2       1       0       L
                                                0b6       5       4       3       2       1       0       L
                                         0b6       5       4       3       2       1       0       L
                                  0b6       5       4       3       2       1       0       L
                           0b6       5       4       3       2       1       0       L
                    0b6       5       4       3       2       1       0       L

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask       =         0b11111110000000000000000000000000000000000000000000000000;
        */
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 7 bits
     * and embedded tile size of 6 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A000000B000000C000000D000000, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     also used      https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 7 bits
     * and embedded tile size of 6 bits
     * @param nTiles the number of tiles of block size 7 embedded in tiled.  maximum nTiles is 8.
     * @return
     */
    public static long sketch7(final long tiled, final int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 8) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 8 for block size 7");
        }

        // kMult=(1<<0) | (1<<6) | (1<<12) | (1<<18) | (1<<24) | (1<<30) | (1<<36)
        //                    6         5         4         3         2         1
        //                   10987654321098765432109876543210987654321098765432109876543210
        final long kMult = 0b00000000000000000000000001000001000001000001000001000001000001L;
        final long kMask = 0b1111111000000000000000000000000000000000000000000L;
        int kShift = 42;
        long sketch = ((tiled * kMult) & kMask) >> kShift;

        if (nTiles < 8) {
            return sketch;
        }

        //TODO: consider more efficient ways for these last few lines

        sketch += (((tiled>>56) & 0b1) << 7);
/*
        // set bit 7 of sketch to the high bit of block 8
        if ((tiled & (1L << 56)) != 0) {
            sketch |= (1L << 7);
        } else {
            sketch &= ~(1L << 7);
        }
*/
        return sketch;

        /*
        for nTiles <=7
                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                                         _000000_0000001000000100000010000001000000100000010000001000000

                                                               0b      6      5      4      3      2      1      0      L
                                                         0b      6      5      4      3      2      1      0      L
                                                   0b      6      5      4      3      2      1      0      L
                                             0b      6      5      4      3      2      1      0      L
                                       0b      6      5      4      3      2      1      0      L
                                 0b      6      5      4      3      2      1      0      L
                           0b      6      5      4      3      2      1      0      L

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =                0b1111111000000000000000000000000000000000000000000;
         */
        /*
        editing for nTiles = 8
                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                                         _00000010000001000000100000010000001000000100000010000001000000
sketch overlaps here:
    so will handle block 7 after sketching the first 7 blocks as above

                                                                7      6      5      4      3      2      1      0      L
                                                          7      6      5      4      3      2      1      0      L
                                                    7      6      5      4      3      2      1      0      L
                                              7      6      5      4      3      2      1      0      L
                                        7      6      5      4      3      2      1      0      L
                                  7      6      5      4      3      2      1      0      L
                            7      6      5      4      3      2      1      0      L
                      7      6      5      4      3      2      1      0      L

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =         0b11111111000000000000000000000000000000000000000000000000;
        */

        /*System.out.printf("\nkMask=\n%63s\n" +
                "kMult=\n%63s\n" +
                "kShift=%d\n" +
                "(tiled * kMult) & kMask=\n%63s\n" +
                "(((tiled * kMult) & kMask) >> kShift)=\n%63s\n",
                Long.toBinaryString(kMask), Long.toBinaryString(kMult), kShift,
                Long.toBinaryString((tiled * kMult) & kMask),
                Long.toBinaryString(((tiled * kMult) & kMask) >> kShift));*/

    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 6 bits
     * and embedded tile size of 5 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A00000B00000C00000D00000, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     also used      https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge

     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 6 bits
     * and embedded tile size of 5 bits
     * @param nTiles the number of tiles of block size 6 embedded in tiled.  the maximum number for the java unsigned
     *               long is 10 tiles of block size 6.
     * @return
     */
    public static long sketch6(long tiled, final int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 10) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 10 for block size 6");
        }

        // kMult=(1<<0) | (1<<5) | (1<<10) | (1<<15) | (1<<20) | (1<<25)
        //               6         5         4         3         2         1
        //              10987654321098765432109876543210987654321098765432109876543210
        long kMult = 0b000000000000000000000000000000000000010000100001000010000100001L;
        long kMask = 0b111111000000000000000000000000000000L;
        long kShift = 30;

        long sketch = ((tiled * kMult) & kMask) >> kShift;

        if (nTiles < 7) {
            return sketch;
        }

        /*
        note that cannot edit the multiplier to use a larger interval of 9 for nTiles=10 because
        the highest bit in the product that needs to be masked to extract the sketch
        would be > 62 bits.

        so will make another sketch for remaining tiles:

        for nTiles 7-10, that is, blocks 6-9, inclusive.
        shift tiled down by 6*6, leaving blocks 6:9.
        */

        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));
        // shift down by the 6 blocks we just sketched:
        tiled >>= 36;
        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));

        // change the shift to reserve space of 6 at the end to merge the 2 sketches:
        kShift -= 6;
        long sketch2 = ((tiled * kMult) & kMask) >> kShift;

        //from https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
        // Merge bits from two values according to a mask
        //   a=0b000000001110  <--- similar to sketch
        //   b=0b101111110000  <--- similar to sketch2
        //mask=0b111111110000
        //r = a ^ ((a ^ b) & mask); bin(r)
        //0b10111110'

        //System.out.printf("nTiles=%d, blockSize=6\n", nTiles);
        //System.out.printf("sketch=%62s\n", Long.toBinaryString(sketch));
        //System.out.printf("sketch=%62s\n", Long.toBinaryString(sketch2));

        sketch = sketch ^ ((sketch ^ sketch2) & 0b111111000000L);

        //System.out.printf("merged=%62s\n", Long.toBinaryString(sketch));

        return sketch;

        /*
        for nTiles <=7
                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                                         000_00000_00000_00000_00000100000100000100000100000100000100000

                                                                                  0b5     4     3     2     1     0     L
                                                                             0b5     4     3     2     1     0     L
                                                                        0b5     4     3     2     1     0     L
                                                                   0b5     4     3     2     1     0     L
                                                              0b5     4     3     2     1     0     L
                                                         0b5     4     3     2     1     0     L

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =                             0b111111000000000000000000000000000000;
         */
         /*
        for nTiles 7 - 10
        shift tiled down by 6*6, leaving blocks 6:9
        mask = 0b111100000000000000000000
        last shift = 20
                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                                                                          _00000100000100000100000100000
                                                                             BLOCK              9     8     7     6
                                                                                           9     8     7     6
                                                                                      9     8     7     6
                                                                                 9     8     7     6

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =                                         0b111100000000000000000000;
         */
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 5 bits
     * and embedded tile size of 4 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A0000B0000C0000D0000, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     also used      https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 5 bits
     * and embedded tile size of 4 bits
     * @param nTiles the number of tiles of block size 5 embedded in tiled.  the maximum number of 12 for nTiles for
     *               block size of 5 is limited by the java unsigned long and the location of the mask bits needed after
     *               the sketch multiplier.
     * @return
     */
    public static long sketch5(long tiled, final int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 12) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 12 for block size 5");
        }

        int kShift = 20;
        long kMult = 0b000000000000000000000000000000000000000000000010001000100010001L;
        long kMask = 0b1111100000000000000000000L;

        long sketch = ((tiled * kMult) & kMask) >> kShift;

        if (nTiles < 6) {
            return sketch;
        }

        //TODO: consider how to implement a word-level parallel shift right to reduce the tiled values to just the flags?
        // involves a division which is expensive?

        // sketch 5 tiles at a time, and merge after each

        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));
        // shift down by the 5 blocks we just sketched:
        tiled >>= 25;
        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));

        // change the shift to reserve space of 5 at the end to merge the 2 sketches:
        kShift -= 5;
        long sketch2 = ((tiled * kMult) & kMask) >> kShift;

        //from https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
        // Merge bits from two values according to a mask
        //   a=0b000000001110  <--- similar to sketch
        //   b=0b101111110000  <--- similar to sketch2
        //mask=0b111111110000
        //r = a ^ ((a ^ b) & mask); bin(r)
        //0b10111110'

        sketch = sketch ^ ((sketch ^ sketch2) & 0b1111100000L);

        if (nTiles < 11) {
            return sketch;
        }

        // one more round of sketch and merge
        tiled >>= 25;
        kShift -= 5;
        sketch2 = ((tiled * kMult) & kMask) >> kShift;

        // sketch is 10 bits, sketch2 is 5 bits
        sketch = sketch ^ ((sketch ^ sketch2) & 0b111110000000000L);

        return sketch;

        /*
        for nTiles <=5
                                                           6         5         4         3         2         1
                                                          10987654321098765432109876543210987654321098765432109876543210
                                                          00_0000_0000_0000_0000_0000_0000_00001000010000100001000010000
                                                                                               4    3    2    1    0    L
                                                                                           4    3    2    1    0    L
                                                                                       4    3    2    1    0    L
                                                                                   4    3    2    1    0    L
                                                                               4    3    2    1    0    L

                                                           6         5         4         3         2         1
                                                          10987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =                                        0b1111100000000000000000000;
         */
        /*
        editing for nTiles <= 11

        overlaps here so need to
        use more than one sketch

                                                           6         5         4         3         2         1
                                                          10987654321098765432109876543210987654321098765432109876543210
                                                          00100001000010000100001000010000100001000010000100001000010000
                                                            B    A    9    8    7    6    5    4    3    2    1    0    L
                                                        B    A    9    8    7    6    5    4    3    2    1    0    L
                                                    B    A    9    8    7    6    5    4    3    2    1    0    L
                                                B    A    9    8    7    6    5    4    3    2    1    0    L
                                            B    A    9    8    7    6    5    4    3    2    1    0    L
                                        B    A    9    8    7    6    5    4    3    2    1    0    L
                                    B    A    9    8    7    6    5    4    3    2    1    0    L
                                B    A    9    8    7    6    5    4    3    2    1    0    L
                            B    A    9    8    7    6    5    4    3    2    1    0    L
                        B    A    9    8    7    6    5    4    3    2    1    0    L
                    B    A    9    8    7    6    5    4    3    2    1    0    L
                B    A    9    8    7    6    5    4    3    2    1    0    L

                                                           6         5         4         3         2         1
                                                         210987654321098765432109876543210987654321098765432109876543210
                                        kMask1      =     0b111111111111000000000000000000000000000000000000000000000000;
         */
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 4 bits
     * and embedded tile size of 3 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A000B000C000D000, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     also used      https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 4 bits
     * and embedded tile size of 3 bits
     * @param nTiles the number of tiles of block size 4 embedded in tiled.  the maximum number of 15 for nTiles for
     *               block size of 4 is limited by the java unsigned long and the location of the mask bits needed after
     *               the sketch multiplier.
     * @return
     */
    public static long sketch4(long tiled, final int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 15) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 15 for block size 4");
        }

        // TODO: more efficient ways to implement this?  intrinsics?

        // kMult=(1<<0) | (1<<3) | (1<<6) | (1<<9)
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        long kMult = 0b00000000000000000000000000000000000000000000000000001001001001L;
        long kMask = 0b1111000000000000L;
        int kShift = 12;
        long sketch = ((tiled * kMult) & kMask) >> kShift;

        if (nTiles < 5) {
            return sketch;
        }

        // 15 tiles, each sketch is 4 tiles, would mean 4 sketches

        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));
        // shift down by the 4 blocks we just sketched:
        tiled >>= 16;
        //System.out.printf("tiled=%62s\n", Long.toBinaryString(tiled));

        // change the shift to reserve space of 4 at the end to merge the 2 sketches:
        kShift -= 4;
        long sketch2 = ((tiled * kMult) & kMask) >> kShift;

        //from https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
        // Merge bits from two values according to a mask
        //   a=0b000000001110  <--- similar to sketch
        //   b=0b101111110000  <--- similar to sketch2
        //mask=0b111111110000
        //r = a ^ ((a ^ b) & mask); bin(r)
        //0b10111110'

        sketch = sketch ^ ((sketch ^ sketch2) & 0b11110000L);

        if (nTiles < 9) {
            return sketch;
        }

        // another round of sketch and merge
        tiled >>= 16;
        kShift -= 4;
        sketch2 = ((tiled * kMult) & kMask) >> kShift;

        // sketch is 8 bits, sketch2 is 4 bits
        sketch = sketch ^ ((sketch ^ sketch2) & 0b111100000000L);

        // one more round of sketch and merge
        tiled >>= 16;
        kShift -= 4;
        sketch2 = ((tiled * kMult) & kMask) >> kShift;

        // sketch is 12 bits, sketch2 is 4 bits
        sketch = sketch ^ ((sketch ^ sketch2) & 0b1111000000000000L);

        return sketch;
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of 3 bits
     * and embedded tile size of 2 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A00B00C00D00, this method would return ABCD.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     also used      https://graphics.stanford.edu/~seander/bithacks.html#MaskedMerge
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 3 bits
     * and embedded tile size of 2 bits
     * @param nTiles the number of tiles of block size 3 embedded in tiled.  the maximum number of 20 for nTiles for
     *               block size of 3 is limited by the java unsigned long and the location of the mask bits needed after
     *               the sketch multiplier.
     * @return
     */
    public static long sketch3(long tiled, final int nTiles) {

        if (nTiles == 0) {
            return 0;
        }
        if (nTiles < 0 || nTiles > 20) {
            throw new UnsupportedOperationException("nTiles must be > 0 and <= 20 for block size 3");
        }

        // same as decoding 4 bit Morton 3D, for z
        // uses binary magic numbers
        long w = tiled >> 2;
        w &= 0x1249249249249249L;
        w = (w ^ (w >> 2))  & 0x30c30c30c30c30c3L;
        w = (w ^ (w >> 4))  & 0xf00f00f00f00f00fL;
        w = (w ^ (w >> 8))  & 0x00ff0000ff0000ffL;
        w = (w ^ (w >> 16)) & 0x00ff00000000ffffL;
        w = (w ^ (w >> 32)) & 0x00000000001fffffL;
        return w;
    }

    /**
     * given a bitarray packed full of tiles separated by flags with a block size of  bits
     * and embedded tile size of 1 bits, extract and return the flags as consecutive bits.
     * e.g. if tiled were A0B0C0D0, this method would return ABCD.
     <pre>
     reference:
     https://stackoverflow.com/questions/30539347/2d-morton-code-encode-decode-64bits
     </pre>
     *
     * @param tiled a bitarray packed full of tiles separated by flags with a block size of 2 bits
     * and embedded tile size of 1 bits
     *
     * @return
     */
    public static long sketch2(final long tiled) {

        // can do this with 5 sketches and 5 merges

        /* else, use a down shift by 1, then 6 binary magic number shifts and masks
        from https://stackoverflow.com/questions/30539347/2d-morton-code-encode-decode-64bits
        x = x & 0x5555555555555555;
        x = (x | (x >> 1))  & 0x3333333333333333;
        x = (x | (x >> 2))  & 0x0F0F0F0F0F0F0F0F;
        x = (x | (x >> 4))  & 0x00FF00FF00FF00FF;
        x = (x | (x >> 8))  & 0x0000FFFF0000FFFF;
        x = (x | (x >> 16)) & 0x00000000FFFFFFFF;
        */

        long sketch = tiled >> 1;
        //                                                              6         5         4         3         2         1
        //                                                             10987654321098765432109876543210987654321098765432109876543210
        sketch = sketch & 0x5555555555555555L;                    //0b101010101010101010101010101010101010101010101010101010101010101
        sketch = (sketch | (sketch >> 1))  & 0x3333333333333333L; // 0b11001100110011001100110011001100110011001100110011001100110011
        sketch = (sketch | (sketch >> 2))  & 0x0F0F0F0F0F0F0F0FL; //   0b111100001111000011110000111100001111000011110000111100001111
        sketch = (sketch | (sketch >> 4))  & 0x00FF00FF00FF00FFL; //       0b11111111000000001111111100000000111111110000000011111111
        sketch = (sketch | (sketch >> 8))  & 0x0000FFFF0000FFFFL; //               0b111111111111111100000000000000001111111111111111
        sketch = (sketch | (sketch >> 16)) & 0x00000000FFFFFFFFL; //                               0b11111111111111111111111111111111

        return sketch;
    }

    /**
     * Given an n-bit value, returns the index of the highest 1 bit within that
     * value.
     *
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value a bitlength number.  note that this is the same as the block size for some contexts.
     * @param bitlength the bitlength of value
     * @return the highest bit set in value.  the bit number is w.r.t. 0.
     * e.g. if bitlength is 8, the return value range is [0,7] inclusive.
     */
    public static long highestBitSetIn(long value, int bitlength) {

        // the method name highestBitSetIn8 is a one-based index, but the bit number
        // returned is w.r.t. a zero-based index.
        // e.g. highestBitSetIn8 is the highest bit in an 8 bit value,
        //      but the returned bit range is [0,7] inclusive
        // switch is based on the block size which s bitlength + 1
        switch (bitlength + 1) {
            case 8: {
                return highestBitSetIn8(value);
            } case 7: {
                return highestBitSetIn7(value);
            } case 6: {
                return highestBitSetIn6(value);
            } case 5: {
                return highestBitSetIn5(value);
            } case 4: {
                return highestBitSetIn4(value);
            } case 3: {
                return highestBitSetIn3(value);
            } case 2: {
                return highestBitSetIn2(value);
            } case 1: {
                return highestBitSetIn1(value);
            } default: {
                throw new UnsupportedOperationException(bitlength + " is not a supported bitlength");
            }
        }
    }

    /**
     * Given an 8-bit value, returns the index of the highest 1 bit within that
     * value.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * This subroutine is where much of the magic happens with regards to the
     * overall algorithm. The idea is that if we can get down to an eight-bit
     * number, we can manually check each power of two that could serve as the
     * most-significant bit. This is actually done using a clever parallel
     * comparison step, describe below.
     * @param value
     * @return
     */
     static long highestBitSetIn8(long value) {

        /* We will again use the parallel comparison technique. To get everything to
         * fit cleanly into a machine word, we'll treat the 8-bit value as actually
         * being 7 bits, since if the top bit is set we immediately know the answer.
         *
         * As a result, our first step is to simply check if the highest bit is set
         * and immediately return the answer if it is.
         */
        if ((value & 0b10000000L) != 0) {
            return 7;
        }

        /* The main observation here is that the MSB of an integer is equal to the
         * number of powers of two less than or equal to it, minus one. As an
         * example, the number 00110110 has six powers of two less than or equal to
         * it (each of the powers of two that are less than or equal to 00100000).
         * Therefore, if we can count up how many powers of two are less than or
         * equal to our number, we can quickly determine the most-significant bit.
         *
         * Since we have a seven-bit number at this point, there are only seven
         * powers of two that we need to test, and we can test them all in parallel
         * using our lovely parallel comparison technique! We'll compare against
         * the numbers 1000000, 100000, 10000, 1000, 100, 10, and 1 all at the same
         * time by performing this subtraction, which is a parallel compare:
         *
         *   1aaaaaaa1aaaaaaa1aaaaaaa1aaaaaaa1aaaaaaa1aaaaaaa1aaaaaaa
         * - 01000000001000000001000000001000000001000000001000000001
         *
         * That bottom string is the concatenation of all the powers of two listed
         * above, padded with zeros betweeh the elements.
         */
                           // 6         5         4         3         2         1
                          // 10987654321098765432109876543210987654321098765432109876543210
        long kComparator = 0b00000001000000001000000001000000001000000001000000001000000001L;

        /* We need to spray out seven copies of the value so that we can compare
         * multiple copies of it in parallel. To do this, we essentially want to make
         * a bunch of shifts and add them all together:
         *
         *                                                     aaaaaaaa
         *                                             aaaaaaaa
         *                                     aaaaaaaa
         *                             aaaaaaaa
         *                     aaaaaaaa
         *             aaaaaaaa
         *     aaaaaaaa
         *  + ---------------------------------------------------------
         *
         * This corresponds to a multiplication by 2^0 + 2^8 + 2^16 + 2^24 + ...
         * It is (1<<0) | (1<<8) | (1<<16) | (1<<24) etc
         */
         long kSpreader   = 0b0000000000000001000000010000000100000001000000010000000100000001L;

        /* As before, to make a parallel comparison, we're going to force a 1 bit at
         * the start of each block. (This is why we special-cased away the top bit of
         * the byte - we need to recycle that bit for other purposes.)
         * It is (1<<7) | (1<<15) | (1<<23) | (1<<31) etc.
         */
         long kMask       = 0b0000000010000000100000001000000010000000100000001000000010000000L;

        /* Perform the parallel comparison:
         *
         *  1. Spray out multiple copies of the value.
         *  2. Put 1 bits at the front of each block.
         *  3. Do the subtraction with the comparator to see which powers of two
         *     we're bigger than.
         *  4. Mask off everything except the flag bits.
         */
         long comparison = (((kSpreader * value) | kMask) - kComparator) & kMask;

        /* We now have a flag integer holding bits indicating which powers of two
         * are smaller than us. Summing up those flags using a summation operation
         * gives us the number we want!
         */
        int nTiles = 7;
        int tileBitLength = 7;
        return parallelSum8(comparison, nTiles) - 1;
    }

    /**
     * Given a 7-bit value, returns the index of the highest 1 bit within that
     * value.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value
     * @return
     */
     static long highestBitSetIn7(long value) {

        if ((value & 0b1000000L) != 0) {
            return 6;
        }

        /*
         * Since we have a 6-bit number at this point, there are only 6
         * powers of two that we need to test, and we can test them all in parallel
         * using our lovely parallel comparison technique! We'll compare against
         * the numbers 100000, 10000, 1000, 100, 10, and 1 all at the same
         * time by performing this subtraction, which is a parallel compare:
         *
         *    1aaaaaa1aaaaaa1aaaaaa1aaaaaa1aaaaaa1aaaaaa
         * -  010000000100000001000000010000000100000001
         *
         * That bottom string is the concatenation of all the powers of two listed
         * above, padded with zeros between the elements.
         */

        //         5         4         3         2         1
        //  7654321098765432109876543210987654321098765432109876543210
        //                  1aaaaaa1aaaaaa1aaaaaa1aaaaaa1aaaaaa1aaaaaa
        //                   100000  10000   1000    100     10      1
        long kComparator = 0b10000000100000001000000010000000100000001L;

        /* We need to spray out 6 copies of the value so that we can compare
         * multiple copies of it in parallel. To do this, we essentially want to make
         * a bunch of shifts and add them all together:
         *
         *                                                     aaaaaaa
         *                                              aaaaaaa
         *                                       aaaaaaa
         *                                aaaaaaa
         *                         aaaaaaa
         *                  aaaaaaa
         *  + ---------------------------------------------------------
         *
         * This corresponds to a multiplication by 2^0 + 2^7 + 2^14 + 2^21 + ...
         * It is (1<<0) | (1<<7) | (1<<14) | (1<<21) | (1<<28) | (1<<35)
         */
        //                  6         5         4         3         2         1
        //                 10987654321098765432109876543210987654321098765432109876543210
        long kSpreader = 0b00000000000000000000000000100000010000001000000100000010000001L;

        /* As before, to make a parallel comparison, we're going to force a 1 bit at
         * the start of each block. (This is why we special-cased away the top bit of
         * the byte - we need to recycle that bit for other purposes.)
         * (1<<6) | (1<<13) | (1<<20) | (1<<27) | (1<<34) | (1<<41)
         */
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        long kMask = 0b00000000000000000000100000010000001000000100000010000001000000L;

        /* Perform the parallel comparison:
         *
         *  1. Spray out multiple copies of the value.
         *  2. Put 1 bits at the front of each block.
         *  3. Do the subtraction with the comparator to see which powers of two
         *     we're bigger than.
         *  4. Mask off everything except the flag bits.
         */
        long comparison = (((kSpreader * value) | kMask) - kComparator) & kMask;

        /* We now have a flag integer holding bits indicating which powers of two
         * are smaller than us. Summing up those flags using a summation operation
         * gives us the number we want!
         */
        int nTiles = 6;
        int tileBitLength = 6;
        return parallelSum(comparison, nTiles, tileBitLength) - 1;
    }

    /**
     * Given a 6-bit value, returns the index of the highest 1 bit within that
     * value.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value
     * @return
     */
     static long highestBitSetIn6(long value) {

        if ((value & 0b100000L) != 0) {
            return 5;
        }

        /*
         * Since we have a 5-bit number at this point, there are only 5
         * powers of two that we need to test, and we can test them all in parallel
         * using our lovely parallel comparison technique! We'll compare against
         * the numbers  10000, 1000, 100, 10, and 1 all at the same
         * time by performing this subtraction, which is a parallel compare:
         *
         *    1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa
         * -  000000000000010000001000000100000010000001
         *
         * That bottom string is the concatenation of all the powers of two listed
         * above, padded with zeros between the elements.
         */

        //         5         4         3         2         1
        //  7654321098765432109876543210987654321098765432109876543210
        //                  1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa1aaaaa
        //                - 000000000000010000001000000100000010000001
        long kComparator = 0b00000000000010000001000000100000010000001L;

        /* We need to spray out 5 copies of the value so that we can compare
         * multiple copies of it in parallel. To do this, we essentially want to make
         * a bunch of shifts and add them all together:

                                         aaaaaa
                                   aaaaaa
                             aaaaaa
                       aaaaaa
                 aaaaaa
           aaaaaa
         ---------------------------------------------------------

         * This corresponds to a multiplication by 2^0 + 2^6 + 2^12 + 2^18 + ...
         * It is (1<<0) | (1<<6) | (1<<12) | (1<<18) | (1<<24)
         */
        //                  6         5         4         3         2         1
        //                 10987654321098765432109876543210987654321098765432109876543210
        long kSpreader = 0b00000000000000000000000000000000000001000001000001000001000001L;

        /* As before, to make a parallel comparison, we're going to force a 1 bit at
         * the start of each block. (This is why we special-cased away the top bit of
         * the byte - we need to recycle that bit for other purposes.)
         * (1<<5) | (1<<11) | (1<<17) | (1<<23) | (1<<29)
         */
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        long kMask = 0b00000000000000000000000000000000100000100000100000100000100000L;

        /* Perform the parallel comparison:
         *
         *  1. Spray out multiple copies of the value.
         *  2. Put 1 bits at the front of each block.
         *  3. Do the subtraction with the comparator to see which powers of two
         *     we're bigger than.
         *  4. Mask off everything except the flag bits.
         */
        long comparison = (((kSpreader * value) | kMask) - kComparator) & kMask;

        /* We now have a flag integer holding bits indicating which powers of two
         * are smaller than us. Summing up those flags using a summation operation
         * gives us the number we want!
         */
        int nTiles = 5;
        int tileBitLength = 5;
        return parallelSum(comparison, nTiles, tileBitLength) - 1;
    }

    /**
     * Given a 5-bit value, returns the index of the highest 1 bit within that
     * value.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value
     * @return
     */
     static long highestBitSetIn5(long value) {

        if ((value & 0b10000L) != 0) {
            return 4;
        }

        /*
         * Since we have a 4-bit number at this point, there are only 4
         * powers of two that we need to test, and we can test them all in parallel
         * using our lovely parallel comparison technique! We'll compare against
         * the numbers  1000, 100, 10, and 1 all at the same
         * time by performing this subtraction, which is a parallel compare:
         *
         *    1aaaa1aaaa1aaaa1aaaa1aaaa1aaaa1aaaa
         * -  00000000000000001000001000001000001
         *
         * That bottom string is the concatenation of all the powers of two listed
         * above, padded with zeros between the elements.
         */

        //         5         4         3         2         1
        //  7654321098765432109876543210987654321098765432109876543210
        //                         1aaaa1aaaa1aaaa1aaaa1aaaa1aaaa1aaaa
        //                       - 00000000000000001000001000001000001
        long kComparator = 0b00000000000000000000001000001000001000001L;

        /* We need to spray out 4 copies of the value so that we can compare
         * multiple copies of it in parallel. To do this, we essentially want to make
         * a bunch of shifts and add them all together:

                                    aaaaa
                               aaaaa
                          aaaaa
                     aaaaa
                aaaaa
           aaaaa
         ---------------------------------------------------------

         * This corresponds to a multiplication by 2^0 + 2^5 + 2^10 + 2^15
         * It is (1<<0) | (1<<5) | (1<<10) | (1<<15)
         */
        //                  6         5         4         3         2         1
        //                 10987654321098765432109876543210987654321098765432109876543210
        long kSpreader = 0b00000000000000000000000000000000000000000000001000010000100001L;

        /* As before, to make a parallel comparison, we're going to force a 1 bit at
         * the start of each block. (This is why we special-cased away the top bit of
         * the byte - we need to recycle that bit for other purposes.)
         * (1<<4) | (1<<9) | (1<<14) | (1<<19)
         */
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        long kMask = 0b00000000000000000000000000000000000000000010000100001000010000L;

        /* Perform the parallel comparison:
         *
         *  1. Spray out multiple copies of the value.
         *  2. Put 1 bits at the front of each block.
         *  3. Do the subtraction with the comparator to see which powers of two
         *     we're bigger than.
         *  4. Mask off everything except the flag bits.
         */
        long comparison = (((kSpreader * value) | kMask) - kComparator) & kMask;

        /* We now have a flag integer holding bits indicating which powers of two
         * are smaller than us. Summing up those flags using a summation operation
         * gives us the number we want!
         */
        int nTiles = 4;
        int tileBitLength = 4;
        return parallelSum(comparison, nTiles, tileBitLength) - 1;
    }

    /**
     * Given a 4-bit value, returns the index of the highest 1 bit within that
     * value.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     * @param value
     * @return
     */
     static long highestBitSetIn4(long value) {

        if ((value & 0b1000L) != 0) {
            return 3;
        }

        /*
         * Since we have a 3-bit number at this point, there are only 3
         * powers of two that we need to test, and we can test them all in parallel
         * using our lovely parallel comparison technique! We'll compare against
         * the numbers  100, 10, and 1 all at the same
         * time by performing this subtraction, which is a parallel compare:
         *
         *    1aaa1aaa1aaa1aaa1aaa1aaa1aaa
         * -     0000000000000010000100001
         *
         * That bottom string is the concatenation of all the powers of two listed
         * above, padded with zeros between the elements.
         */

        //         5         4         3         2         1
        //  7654321098765432109876543210987654321098765432109876543210
        //                                1aaa1aaa1aaa1aaa1aaa1aaa1aaa
        //                                 - 0000000000000010000100001
        long kComparator = 0b00000000000000000000000000000010000100001L;

        /* We need to spray out 4 copies of the value so that we can compare
         * multiple copies of it in parallel. To do this, we essentially want to make
         * a bunch of shifts and add them all together:

                                    aaaa
                                aaaa
                            aaaa
                        aaaa
                    aaaa
                aaaa
         ---------------------------------------------------------

         * This corresponds to a multiplication by 2^0 + 2^4 + 2^8
         * It is (1<<0) | (1<<4) | (1<<8)
         */
        //                   6         5         4         3         2         1
        //                 10987654321098765432109876543210987654321098765432109876543210
        long kSpreader = 0b00000000000000000000000000000000000000000000000000000100010001L;

        /* As before, to make a parallel comparison, we're going to force a 1 bit at
         * the start of each block. (This is why we special-cased away the top bit of
         * the byte - we need to recycle that bit for other purposes.)
         * (1<<3) | (1<<7) | (1<<11)
         */
        //              6         5         4         3         2         1
        //             10987654321098765432109876543210987654321098765432109876543210
        long kMask = 0b00000000000000000000000000000000000000000000000000100010001000L;

        /* Perform the parallel comparison:
         *
         *  1. Spray out multiple copies of the value.
         *  2. Put 1 bits at the front of each block.
         *  3. Do the subtraction with the comparator to see which powers of two
         *     we're bigger than.
         *  4. Mask off everything except the flag bits.
         */
        long comparison = (((kSpreader * value) | kMask) - kComparator) & kMask;

        /* We now have a flag integer holding bits indicating which powers of two
         * are smaller than us. Summing up those flags using a summation operation
         * gives us the number we want!
         */
        int nTiles = 3;
        int tileBitLength = 3;
        return parallelSum(comparison, nTiles, tileBitLength) - 1;
    }

    /**
     * Given a 3-bit value, returns the index of the highest 1 bit within that
     * value.
     * @param value
     * @return
     */
     static long highestBitSetIn3(long value) {
        if ((value & 0b100L) != 0) {
            return 2;
        } else if ((value & 0b10L) != 0) {
            return 1;
        } else if ((value & 0b1L) != 0) {
            return 0;
        }
        throw new IllegalArgumentException("no bits below 3 are set in value");
    }
     static long highestBitSetIn2(long value) {
        if ((value & 0b10L) != 0) {
            return 1;
        } else if ((value & 0b1L) != 0) {
            return 0;
        }
        throw new IllegalArgumentException("no bits below 3 are set in value");
    }
     static long highestBitSetIn1(long value) {
        if ((value & 0b1L) != 0) {
            return 0;
        }
        throw new IllegalArgumentException("no bits below 3 are set in value");
    }

    /**
     *  Returns a bitmask where each block's high bit is 1 if that block contains a
     * 1 bit and is 0 otherwise. All remaining bits are 0.
     *
     * Stated differently, given the input
     *
     *   aaaaaaaabbbbbbbbccccccccddddddddeeeeeeeeffffffffgggggggghhhhhhhh
     *
     * We'll return a 64-bit flag integer
     *
     *   A0000000B0000000C0000000D0000000E0000000F0000000G0000000H0000000
     *
     * where each letter is 1 if any of the bits in the block were set and is
     * 0 otherwise.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     */
    public static long usedBlocksIn(long value, int blockSize) {
        switch(blockSize) {
            case 8:
                return usedBlocksIn8(value);
            case 7:
                return usedBlocksIn7(value);
            case 6:
                return usedBlocksIn6(value);
            default:
                throw new UnsupportedOperationException("not yet implemented");
        }
    }

    /**
     * For a tiled value whose block size = 8 bits (and hence, the embedded tile size is 7 bits in between
     * flags of size 1), returns
     * a bitmask where each block's high bit is 1 if that block contains a
     * 1 bit and is 0 otherwise. All remaining bits are 0.
     *
     * Stated differently, given the input
     *
     *   aaaaaaaabbbbbbbbccccccccddddddddeeeeeeeeffffffffgggggggghhhhhhhh
     *
     * We'll return a 64-bit flag integer
     *
     *   A0000000B0000000C0000000D0000000E0000000F0000000G0000000H0000000
     *
     * where each letter is 1 if any of the bits in the block were set and is
     * 0 otherwise.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     */
    public static long usedBlocksIn8(long value) {

        /* Every block with a 1 bit set in it either
         *  1. has its highest bit set, or
         *  2. when that bit is cleared, has a numeric value of 1 or greater.
         *
         * We can check for ths first part by using a bitmask to extract the high bits
         * from each of the blocks. The remainder can be identified by using the
         * parallel comparison technique of comparing each block against 1.
         */
        // Positions of all the high bits within each block. (1<<7)|(1<<15)|(1<<23)|(1<<31)|(1<<39)|(1<<47)|(1<<55)|(1<<64)
        final long kHighBits = 0b000000010000000100000001000000010000000100000001000000010000000L;
    //    final long kHighBits =0b1000000010000000100000001000000010000000100000001000000010000000L;
        //                       6         5         4         3         2         1
        //                      10987654321098765432109876543210987654321098765432109876543210
        long highBitsSet = value & kHighBits;

        /* Now, do a parallel comparison on the 7-bit remainders of each block to
         * identify all the blocks with a nonzero bit set in them.
         *
         * The parallel comparison works as follows. We begin by reshaping the blocks
         * so that each block starts with a 1:
         *
         *   1aaaaaaa1bbbbbbb1ccccccc1ddddddd1eeeeeee1fffffff1ggggggg1hhhhhhh
         *
         * Now, subtract out the value with 1's at the bottom of each block:
         *
         *   1aaaaaaa 1bbbbbbb 1ccccccc 1ddddddd 1eeeeeee 1fffffff 1ggggggg 1hhhhhhh
         * - 00000001 00000001 00000001 00000001 00000001 00000001 00000001 00000001
         *
         * If a block is nonempty, then the subtraction will stop before hitting the
         * special 1 bit we placed at the front of the block. That 1 bit then means
         * "yes, there was some bit set here." Otherwise, if the block is empty, then
         * the subtraction within that block will be forced to borrow the 1 bit from
         * the flag, which means that the resulting 0 bit means "no, there was no bit
         * set here."
         *
         * We can therefore perform the subtraction, mask off all the bits except for
         * the flags, and we end up with what we're looking for.
         */
        // (1<<0) | (1<<8) | (1<<16) | (1<<24)| (1<<32) | (1<<40) | (1<<48) | (1<<56)
        //                   6         5         4         3         2         1
        //                  10987654321098765432109876543210987654321098765432109876543210
         long kLowBits =  0b00000100000001000000010000000100000001000000010000000100000001L;
         long lowBitsSet  = ((value | kHighBits) - kLowBits) & kHighBits;

        /* Combine them together to find nonempty blocks. */
        return highBitsSet | lowBitsSet;
    }

    /**
     * For a tiled value whose block size = 7 bits (and hence, the embedded tile size is 6 bits in between
     * flags of size 1), returns
     * a bitmask where each block's high bit is 1 if that block contains a
     * 1 bit and is 0 otherwise. All remaining bits are 0.
     *
     * Stated differently, given the input
     *
     *   aaaaaaabbbbbbbcccccccdddddddeeeeeeefffffffggggggghhhhhhh
     *
     * We'll return a 64-bit flag integer
     *
     *   A000000B000000C000000D000000E000000F000000G000000H000000
     *
     * where each letter is 1 if any of the bits in the block were set and is
     * 0 otherwise.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     */
    public static long usedBlocksIn7(long value) {

        /* Every block with a 1 bit set in it either
         *  1. has its highest bit set, or
         *  2. when that bit is cleared, has a numeric value of 1 or greater.
         *
         * We can check for ths first part by using a bitmask to extract the high bits
         * from each of the blocks. The remainder can be identified by using the
         * parallel comparison technique of comparing each block against 1.
         */
        // Positions of all the high bits within each block.
        // (1<<6)|(1<<13)|(1<<20)|(1<<27)|(1<<34)|(1<<41)|(1<<48)|(1<<55)|(1<<62)
        final long kHighBits = 0b100000010000001000000100000010000001000000100000010000001000000L;
        //                        6         5         4         3         2         1
        //                       10987654321098765432109876543210987654321098765432109876543210
        long highBitsSet = value & kHighBits;

        /* Now, do a parallel comparison on the 7-bit remainders of each block to
         * identify all the blocks with a nonzero bit set in them.
         *
         * The parallel comparison works as follows. We begin by reshaping the blocks
         * so that each block starts with a 1:
         *
         *   1aaaaaa1bbbbbb1cccccc1dddddd1eeeeee1ffffff1gggggg1hhhhhh
         *
         * Now, subtract out the value with 1's at the bottom of each block:
         *
         *   1aaaaaa 1bbbbbb 1cccccc 1dddddd 1eeeeee 1ffffff 1gggggg 1hhhhhh
         * - 0000001 0000001 0000001 0000001 0000001 0000001 0000001 0000001
         *
         * If a block is nonempty, then the subtraction will stop before hitting the
         * special 1 bit we placed at the front of the block. That 1 bit then means
         * "yes, there was some bit set here." Otherwise, if the block is empty, then
         * the subtraction within that block will be forced to borrow the 1 bit from
         * the flag, which means that the resulting 0 bit means "no, there was no bit
         * set here."
         *
         * We can therefore perform the subtraction, mask off all the bits except for
         * the flags, and we end up with what we're looking for.
         */
        // (1<<0) | (1<<7) | (1<<14) | (1<<21)| (1<<28) | (1<<35) | (1<<42) | (1<<49) | (1<<56)
        //                    6         5         4         3         2         1
        //                   10987654321098765432109876543210987654321098765432109876543210
        long kLowBits =    0b00000100000010000001000000100000010000001000000100000010000001L;
        long lowBitsSet  = ((value | kHighBits) - kLowBits) & kHighBits;

        /* Combine them together to find nonempty blocks. */
        return highBitsSet | lowBitsSet;
    }

    /**
     * For a tiled value whose block size = 6 bits (and hence, the embedded tile size is 5 bits in between
     * flags of size 1), returns
     * a bitmask where each block's high bit is 1 if that block contains a
     * 1 bit and is 0 otherwise. All remaining bits are 0.
     *
     * Stated differently, given the input
     *
     *   aaaaaabbbbbbccccccddddddeeeeeeffffffgggggghhhhhh
     *
     * We'll return a 64-bit flag integer
     *
     *   A00000B00000C00000D00000E00000F00000G00000H00000
     *
     * where each letter is 1 if any of the bits in the block were set and is
     * 0 otherwise.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     and code in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     Then edited here to allow block sizes other than 8.
     </pre>
     */
    public static long usedBlocksIn6(long value) {

        /* Every block with a 1 bit set in it either
         *  1. has its highest bit set, or
         *  2. when that bit is cleared, has a numeric value of 1 or greater.
         *
         * We can check for ths first part by using a bitmask to extract the high bits
         * from each of the blocks. The remainder can be identified by using the
         * parallel comparison technique of comparing each block against 1.
         */
        // Positions of all the high bits within each block.
        // (1<<5)|(1<<11)|(1<<17)|(1<<23)|(1<<29)|(1<<35)|(1<<41)|(1<<47)|(1<<53) |(1<<59)
        final long kHighBits = 0b00100000100000100000100000100000100000100000100000100000100000L;
        //                       10987654321098765432109876543210987654321098765432109876543210
        //                        6         5         4         3         2         1

        long highBitsSet = value & kHighBits;

        /* Now, do a parallel comparison on the 7-bit remainders of each block to
         * identify all the blocks with a nonzero bit set in them.
         *
         * The parallel comparison works as follows. We begin by reshaping the blocks
         * so that each block starts with a 1:
         *
         *   1aaaaa1bbbbb1ccccc1ddddd1eeeee1fffff1ggggg1hhhhh
         *
         * Now, subtract out the value with 1's at the bottom of each block:
         *
         *   1aaaaa 1bbbbb 1ccccc 1ddddd 1eeeee 1fffff 1ggggg 1hhhhh
         * - 000001 000001 000001 000001 000001 000001 000001 000001
         *
         * If a block is nonempty, then the subtraction will stop before hitting the
         * special 1 bit we placed at the front of the block. That 1 bit then means
         * "yes, there was some bit set here." Otherwise, if the block is empty, then
         * the subtraction within that block will be forced to borrow the 1 bit from
         * the flag, which means that the resulting 0 bit means "no, there was no bit
         * set here."
         *
         * We can therefore perform the subtraction, mask off all the bits except for
         * the flags, and we end up with what we're looking for.
         */
        // (1<<0) | (1<<6) | (1<<12) | (1<<18)| (1<<24) | (1<<30) | (1<<36) | (1<<42) | (1<<48) | (1<<54) | (1<<60)
        //                   6         5         4         3         2         1
        //                  10987654321098765432109876543210987654321098765432109876543210
        long kLowBits =   0b01000001000001000001000001000001000001000001000001000001000001L;
        long lowBitsSet  = ((value | kHighBits) - kLowBits) & kHighBits;

        /* Combine them together to find nonempty blocks. */
        return highBitsSet | lowBitsSet;
    }
}
