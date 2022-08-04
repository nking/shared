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

    //TODO: optimize code to reuse bit masks and multipliers.  consider a static array of
    //      masks and multipliers for block sizes 2 through 31.

    /**
     *
     * given an array of bitstringLength values, concatenate them and insert 0's on the high
     * end of each value.
     * e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     *     tiled1 is 0b0010010001100111, where 0's have been concatenated onto the high end of
     *     each tileBitLength bitstring, making a bitstring of length 16.

     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (values.length * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent values.length is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from values array in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.

     * @param values array of bitstrings, each of length bitstringLength
     * @param bitstringLength the bitlength of each value in values.  the tile for each will be bitstringLength + 1
     *                        bits long.  the total tiled result will be values.length * (bitstringLength + 1) bits.
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
     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (nTiles * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.
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
     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (nTiles * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.
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
        //System.out.printf("value=%s, tiled=%s\n", Integer.toBinaryString(value), Long.toBinaryString(kTiled));
        return kTiled;
    }

    /**
    parallel compare of tiled1 to tiled2 and return a masked bit array whose set bits indicate which
    tiles of tiled1 are .gte. the tiles of tiled2 in the same position.
     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     </pre>
     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (nTiles * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.
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
     Then edited here to allow block sizes other than 8.
     </pre>
     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (nTiles * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.
     @param tiled1 a bit array holding numbers of length tileBitLength (called tiles) separated by 1's.
     e.g. For bitstrings 0b0100100 and 0b1100111 which are 7 bits long,
     tiled1 is 0b1010010011100111, where 1's have been concatenated onto the high end of
     each tileBitLength bitstring, making a bitstring of length 16.
     @param tiled2 a bit array holding numbers of length tileBitLength separated by 0's.
     @param tileBitLength the length of each tile in the bit arrays.  the block size is tileBitLength + 1
     because it includes the gap bit between tiles.
     @param mask1 the 1's mask (same used in setting the gap bits in tiled1)
     @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     in tiled2.
     */
    public static long parallelCompare10(long tiled1, long tiled2, int nTiles, int tileBitLength, long mask1) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        if (nTiles < 1) {
            throw new IllegalArgumentException("nTiles must be > 0");
        }

        final int bSz = tileBitLength + 1;

        final int nMaskBits = (int)Math.ceil(Math.log(nTiles)/Math.log(2));

        // by default the number of blacks used to hold the number nTiles is 1.
        // if number of bits in nTiles > bSz, nBExtra is the number of blacks to add to the default '1 reserved block'
        // to hold the number of bits in nTiles.
        // e.g. nTiles=7 can be held in a block of size 3 bits. if bSz=2, need 1 extra block to hold nTiles.
        //      nTiles=16 can be held in a block of size 5 bits. if bSz=2, need 2 extra blocks to hold nTiles.
        int nBExtra = 0;

        // assert that there is enough space to hold the bits to represent nTiles
        if (nMaskBits > bSz) {
            // how many blocks needed to store nMaskBits?  then subtract 1 which is already reserved for it.
            nBExtra = ((int)Math.ceil((double)nMaskBits/bSz)) - 1;
            int tiledLength = nTiles * bSz;
            if ((tiledLength + nBExtra*bSz) > 63) {
                throw new IllegalArgumentException(String.format("nTiles needs %d blocks of size %d bits above the " +
                        "total tiled bit length =%d.\n  That total %d must fit within 63 bits.",
                        nBExtra, bSz, tiledLength, (tiledLength + nBExtra*bSz)));
            }
        }

        //3. Compute X – Y. The bit preceding xi – yi is 1 if xi ≥ yi and 0 otherwise.
        long diff = tiled1 - tiled2;

        long comparison = diff & mask1;

        //System.out.printf("tiled1=%30s\ntiled2=%30s\ndiff=%32s\ncomp=%32s\n", Long.toBinaryString(tiled1),
        //        Long.toBinaryString(tiled2), Long.toBinaryString(diff), Long.toBinaryString(comparison));

        return parallelSum(comparison, nTiles, tileBitLength);
    }

    /**
     sum the set bits of bitstring comparison.  the flags that may have set bits are the MSB if each block.

     <pre>
     following lecture notes http://web.stanford.edu/class/cs166/lectures/16/Small16.pdf
     Then edited here to allow block sizes other than 8.
     </pre>
     NOTE that there are some size restrictions to the packing especially in context of further use such as the compare
     operations.
     Let block size = (bistringLength + 1).
     The unsigned long restricts the total bit length of the tiled result of this method to 63 bits,
     and so (nTiles * block) must be less than or equal to 63.
     Also, regarding the number of values to be tiled: the compare operation has to be able to store the bit
     representation of the number of tiles into the highest blocks of a mask that is the same size as the
     total tiled bit length.  If the number of bits needed to represent nTiles is not less than or equal to
     block size, then more blocks are needed to hold that number and that number of extra blocks may need to be subtracted
     from nTiles in order for the compare bitMask to fit within the limits of the tiled bit length
     and the 63 bit limit.
     @param comparison a bit array with flags at the MSB of each block.  The flags that are set bits
     are summed in this method.
     @param tileBitLength the length of each tile in the bit arrays.  the block size is tileBitLength + 1
     because it includes the gap bit between tiles.
     @return a bit array of same size as tiled1 and tiled2 in which the bit of each
     tile is 1 if the tile in tiled1 1 is greater than or equal to the tile at the same position
     in tiled2.
     */
    static long parallelSum(long comparison, int nTiles, int tileBitLength) {

        //following sumOf in http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
        // then edited to make the block size variable

        final int bSz = tileBitLength + 1;

        final int nMaskBits = (int)Math.ceil(Math.log(nTiles)/Math.log(2));

        // by default the number of blacks used to hold the number nTiles is 1.
        // if number of bits in nTiles > bSz, nBExtra is the number of blacks to add to the default '1 reserved block'
        // to hold the number of bits in nTiles.
        // e.g. nTiles=7 can be held in a block of size 3 bits. if bSz=2, need 1 extra block to hold nTiles.
        //      nTiles=16 can be held in a block of size 5 bits. if bSz=2, need 2 extra blocks to hold nTiles.
        int nBExtra = 0;

        // assert that there is enough space to hold the bits to represent nTiles
        if (nMaskBits > bSz) {
            // how many blocks needed to store nMaskBits?  then subtract 1 which is already reserved for it.
            nBExtra = ((int)Math.ceil((double)nMaskBits/bSz)) - 1;
            int tiledLength = nTiles * bSz;
            if ((tiledLength + nBExtra*bSz) > 63) {
                throw new IllegalArgumentException(String.format("nTiles needs %d blocks of size %d bits above the " +
                                "total tiled bit length =%d.\n  That total %d must fit within 63 bits.",
                        nBExtra, bSz, tiledLength, (tiledLength + nBExtra*bSz)));
            }
        }

        int i0 = 0;
        // e.g. for bSz=8, kMult=(1<<8)|(1<<0) etc
        long kMult = 0;
        int i;
        for (i = 0; i < (nTiles-1); ++i) {
            kMult |= (1L << i0);
            i0 += bSz;
        }

        int nBitsInMask = (int)Math.ceil(Math.log(nTiles)/Math.log(2));
        long kMask = 0;
        int b2 = ((nTiles - (1 + nBExtra)) * bSz) - 1;
        for (i = 0; i < nBitsInMask; ++i) {
            kMask |= (1L << (b2 + i));
        }
        int kShift = (bSz * nTiles) - (bSz * (nBExtra + 1)) - 1;
        int kShift2 = (bSz * nTiles) - 1 - (nBExtra * bSz);

        long s1 = (((comparison * kMult) & kMask) >> kShift);
        long s2 = (comparison >> kShift2);

        /*System.out.printf("\nkMask= %30s\nkMult= %30s\nkShift=%d\nkShift2=%d\n" +
                "nBExtra=%d\n" +
                "(((comparison * kMult) & kMask) >> kShift)=\n%37s\n=%d" +
                "\n(comparison >> kShift2)=\n%37s\n=%d\n",
                Long.toBinaryString(kMask), Long.toBinaryString(kMult), kShift, kShift2, nBExtra,
                Long.toBinaryString(s1), s1, Long.toBinaryString(s2), s2);
         */

        long sum = s1 + s2;

        return sum;

        /*
        case 1:
            caveat: the first example is for a 64-bit string, and in java we're limited to unsigned 63 bit length.
            sumOf for 8 bit tiling of 7 bit bitstrings:
                       6         5         4         3         2         1
                    3210987654321098765432109876543210987654321098765432109876543210
            value=0b1000000010000000100000000000000000000000000000000000000000000000#<== 8 tiles comparison
            kMult=0b0000000000000001000000010000000100000001000000010000000100000001#<== 7 set bits
            kMask=0b0000001110000000000000000000000000000000000000000000000000000000#<== 3 bit mask
            kShift = 64 - 8 - 1
            sum=(((value * kMult) & kMask) >> kShift) + (value >> 63)

        case 2:
            sumOf for 6 bit tiling of 5 bit bitstrings, 10 tiles (total was 64 bits, is now 60 bits):
                       6         5         4         3         2         1
                    3210987654321098765432109876543210987654321098765432109876543210
            value=0b0000100000100000100000100000100000100000100000100000100000100000 #<== 10 tiles comparison
            kMult=0b0000000000000001000001000001000001000001000001000001000001000001 #<== 9 set bits
            kMask=0b0000000111100000000000000000000000000000000000000000000000000000 #<== 4 bit mask
                     # the number of tiles fits in 4 bits, so set the 4 bits at start of 2nd block
            kShift = 60 - 6 - 1
            sum=(((value * kMult) & kMask) >> kShift) + (value >> 59)

        case 3:
            sumOf for 3 bit tiling of 2 bit bitstrings.  21 tiles (total was 64 bits, is now 63 bits.
            because 21 is 5 bits, need to reserve an extra bit at top of array, so can only
            pack 20 tiles into the tiled bit array.
                       6         5         4         3         2         1
                    3210987654321098765432109876543210987654321098765432109876543210
            value=0b0000100100100100100100100100100100100100100100100100000000000000 # 20 tiles comparison, 16 set
            kMult=0b0000000001001001001001001001001001001001001001001001001001001001 # 19 set bits
            kMask=0b0000001111100000000000000000000000000000000000000000000000000000 # 5 bit mask, set from 3rd block?
            # need an additional block shift for 5 bit mask:
            kShift = 60 - 6 - 1 # previously, for 3 bit mask, kShift was tiledBitLength - (tileBitLength + 1) - 1
            kShift2 = 60 - 1 - 3 # previously, for 3 bit mask, kShift2 was tiledBitLength - 1
            sum=(((value * kMult) & kMask) >> kShift) + (value >> kShift2)
        */
    }

    /**
     * given a bitarray packed full of tiles separated by flags, extract and return the flags.
     * e.g. if tiled were A0000000B0000000C0000000D0000000, this method would return ABCD.
     * <pre>
     *     reference http://web.stanford.edu/class/cs166/lectures/16/code/msb64/MSB64.cpp
     *     then edited here for variable block size and number of tiles packed into tiled.
     * </pre>
     * @param tiled a bitarray of concatenated bitstrings of length tileBiltLength sepearated by flag bits.
     *              the portion of tiled read is the first nTiles * (tileBitLength + 1) bits.
     * @param nTiles the number of tiles packed into the bitarray tiled.
     * @param tileBitLength the size of a tile before a gap is appended to it.  the block size is tileBitlength + 1.
     * @return
     */
    public static long sketch(long tiled, int nTiles, int tileBitLength) {

        /*
        Case 1:
           example in the MSB64.cpp sketchOf comments:
           8 bit tiling, 7 bit bitstrings, nTiles=8
                           6         5         4         3         2         1
                        3210987654321098765432109876543210987654321098765432109876543210
        value1      = 0b1000000010000000100000001000000010000000000000000000000000000000;#8 positions, 5 are set
        #               1       2       3       4       5       6       7       8
        kMult1      = 0b0000000000000010000001000000100000010000001000000100000010000001;#8 flags set: 0,7,14,21,...<== smallest interval of 7 gives maske size 8 which is enough to hold the 1 bit each to represent a set tile ( interval <= nTiles)
        #                             8      7      6      5      4      3      2      1
        kMask1      = 0b1111111100000000000000000000000000000000000000000000000000000000;#8 bits masked
        kShift1     = 64 - 8;
        sketch1=((value1 * kMult1) & kMask1) >> kShift1;
        bin((value1 * kMult1)); bin(kMask1); bin((value1 * kMult1) & kMask1)

        the multiplication:
                                                        '0b1000000010000000100000001000000010000000000000000000000000000000'
                                                  0b1000000010000000100000001000000010000000000000000000000000000000'
                                           0b1000000010000000100000001000000010000000000000000000000000000000'
                                    0b1000000010000000100000001000000010000000000000000000000000000000'
                             0b1000000010000000100000001000000010000000000000000000000000000000'
                      0b1000000010000000100000001000000010000000000000000000000000000000'
               0b1000000010000000100000001000000010000000000000000000000000000000'
        0b1000000010000000100000001000000010000000000000000000000000000000'
        =
        0b10000001100000111000011110001111100111110011111001111100011110000111000001100000010000000000000000000000000000000'
                                           kMask1      = 0b1111111100000000000000000000000000000000000000000000000000000000;

        Case 2:
            let block size = 7, and the bitstrings in between the flags are 6 bits in length.
            the packing is such that there are 9 tiles in the tiled bitstring of size 63.

            The 9 tiles are 9 set bits and so cannot fit within a block size of 7 using the
            multiplier above which uses (block size - 1) intervals.
            So one has to use a large enough interval in the multiplier so that the stacked 9 tiles
            are sequential and can be masked by a 9 bit mask.
            Therefore, the multiplier mask should have intervals 0, 8, 16, 24, 32, 40, 48, 56, 64.
            On Java platform, we have a limit for unsigned long of 63 bits.
            So we could use BigInteger instead of long, but then we incur an increasingly large cost
            with each operation because that class creates a new object every time instead of
            modifying the properties of one instance.  (I have a class to do that to replace BigInteger,
            but will continue using java primitives for this sketch method instead).

            Sticking with java's long,
            we then have the problem that the original tiled number can only have contained nTiles=8,
            which changes our multiplier one more time in order to not overflow and to be able
            to sketch 8 bits sequentially.

            block size=7, bitstring length of each tile in between flags is 6, and nTiles <= 8.
            then the multiplier is 0, 7, 14, 21, 28, 35, 42, 49
            which is the same multiplier as above in Case 1.
            the mask location and shift may need edits (include details here).

        Case 3:
            3 bit tiling of 2 bit bitstrings.  21 tiles (total was 64 bits, is now 63 bits.

                       6         5         4         3         2         1
                    3210987654321098765432109876543210987654321098765432109876543210
            value=0b0000100100100100100100100100100100100100100100100100000000000000 # 20 tiles comparison, 16 set
            kMult=0b0000000001001001001001001001001001001001001001001001001001001001 # 19 set bits
         */
        throw new UnsupportedOperationException("not yet implemented");
    }
}
