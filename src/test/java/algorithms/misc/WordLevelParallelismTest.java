package algorithms.misc;

import algorithms.SubsetChooser;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class WordLevelParallelismTest extends TestCase {

    public WordLevelParallelismTest(String testName) {
        super(testName);
    }


    @SuppressWarnings("fallthrough")
    public void estHighestBitIn() {
        long h;
        long v;
        int j;
        int bl;
        for (int i = 0; i < 8; ++i) {
            v = 0 | (1L << i);
            bl = i + 1;
            h = WordLevelParallelism.highestBitSetIn(v, bl);
            assertEquals(i, h);
        }

    }

    public void estHighestOneBitIn() {

        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        seed = 232949844799850L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);

        long tiled;
        int nTiles;
        int tileBitLength;
        long h;
        long v;
        // b is block size in bits
        // bn is block number
        // hbn is highest block number having any set bits in it.  it's a zero-based number.
        // v is the tiled bitarray with blocks of size bn
        // nb is the number of blocks in the tiled bitarray.  nb*b <= 62
        // h is the calculation of highestOneBitIn which should be == to hbn
        int nb;
        int r;
        long t;
        System.out.println("   6         5         4         3         2         1");
        System.out.println("  10987654321098765432109876543210987654321098765432109876543210");
        //                        0b10000000000000000000000000000000000000000000000000000000
        //                          *       *       *       *       *       *       *
        for (int b = 8; b < 9; ++b) {
            nb = (int)Math.floor(62./b);
            for (long hbn = nb-1; hbn >= 0; hbn--) {
                // set the bits in v:
                //     randomly choose a number between 0 and b, exclusive and set that within block hbn.
                //     do the same for the lower blocks in v.
                r = rand.nextInt(b);
                t = 0L | (1L << r);
                v = (t << (hbn*b));
                System.out.printf("r=%d, b=%d, hbn=%d\n", r, b, hbn);
                System.out.flush();
                //System.out.printf("v=%62s\n", Long.toBinaryString(v));
                h = WordLevelParallelism.highestOneBitIn(v, nb, b-1);
                assertEquals(hbn, h);
            }
        }

        //

    }

    public void estsHighestBlockSetIn() {
        //               6         5         4         3         2         1
        //             210987654321098765432109876543210987654321098765432109876543210
        long tiled = 0b000000010000000100000001000000010000000100000001000000010000000L;
        int nTiles = 8;
        int tileBitLength = 7;

        long h = WordLevelParallelism.highestBlockSetIn(tiled, nTiles, tileBitLength);
        assertEquals(6L, h);
    }

    public void estUsedBlocksIn() {
        //                   6         5         4         3         2         1
        //                 210987654321098765432109876543210987654321098765432109876543210
        long expected = 0b0000000010000000100000001000000010000000100000001000000010000000L;
        long value    = 0b0000000010000000100000001000000010000000100000001000000010000000L;
        long u = WordLevelParallelism.usedBlocksIn(value, 8);
        assertEquals(expected, u);

        //          6         5         4         3         2         1
        //        210987654321098765432109876543210987654321098765432109876543210
        value = 0b100000001000000010000001000000010000000100000000000100000000011L;
        u = WordLevelParallelism.usedBlocksIn(value, 8);
        assertEquals(expected, u);

        expected = 0b0000000010000001000000100000010000001000000100000010000001000000L;
        value = 0b0000000010000001000000100000010000001000000100000010000001000000L;
        u = WordLevelParallelism.usedBlocksIn(value, 7);
        assertEquals(expected, u);

        value = 0b0000000001000000100000110001000001100000100000010000000010000001L;
        u = WordLevelParallelism.usedBlocksIn(value, 7);
        assertEquals(expected, u);

        expected = 0b0100000010000001000000100000010000001000000100000010000001000000L;
        value    = 0b0100000010000001000000100000010000001000000100000010000001000000L;
        u = WordLevelParallelism.usedBlocksIn(value, 7);
        assertEquals(expected, u);

        value    = 0b0000000100000010000001000000100000010000001000000100000010000001L;
        u = WordLevelParallelism.usedBlocksIn(value, 7);
        assertEquals(expected, u);

        expected = 0b0000100000100000100000100000100000100000100000100000100000100000L;
        value    = 0b0000100000100000100000100000100000100000100000100000100000100000L;
        u = WordLevelParallelism.usedBlocksIn(value, 6);
        assertEquals(expected, u);

        value    = 0b0001000001000001000001000001000001000001000001000001000001000001L;
        u = WordLevelParallelism.usedBlocksIn(value, 6);
        assertEquals(expected, u);

        value    = 0b0010000010000010000010000010000010000010000010000010000010000010L;
        u = WordLevelParallelism.usedBlocksIn(value, 6);
        assertEquals(expected, u);
    }

    public void testSketch() {

        // the sketch input is a tiled bit array of all 0's except the flag bits.
        //  it extracts the flag bits and concatenates them and returns that.

        SubsetChooser chooser;
        long comp;
        // b is block size in bits
        // nb is the number of blocks in the tiled bitarray.  nb*b <= 62.  this is nTiles.
        // tileBitLength is b-1
        int nb;
        long sketch;
        long expectedSketch;
        int[] selectedIndexes;
        for (int b = 8; b > 6; --b) {//8

            nb = (int) Math.floor(62. / b);

            // test all combinations of the blocks' flag bits for  [0, nb-1] inclusive
            for (int k = 0; k < nb; ++k) {
                //System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                comp = 0L;
                if (k == 0) {
                    sketch = WordLevelParallelism.sketch(comp, nb, b-1);
                    //System.out.printf("comp=%62s\n", Long.toBinaryString(comp));
                    //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
                    continue;
                }
                selectedIndexes = new int[k];
                chooser = new SubsetChooser(nb-1, k);
                while (chooser.getNextSubset(selectedIndexes) != -1) {
                    comp = 0L;
                    expectedSketch = 0L;
                    // set the high bit of each block in selectedIndexes
                    for (int i = 0; i < selectedIndexes.length; ++i) {
                        //         6         5         4         3         2         1
                        //       210987654321098765432109876543210987654321098765432109876543210
                        //              1_______1_______1_______1_______1_______1_______1_______
                        // block 0, b=8 =>   1_______  bit 7
                        // block 1, b=8 =>   1_______0_______ bit 15

                        comp |= (1L << ((b-1) + selectedIndexes[i] * b));
                        expectedSketch |= (1L << selectedIndexes[i]);
                    }

                    // test that the sketch finds the same set bits
                    sketch = WordLevelParallelism.sketch(comp, nb, b-1);

                    if (expectedSketch != sketch) {
                        System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                        System.out.printf("comp=%62s\n", Long.toBinaryString(comp));
                        System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        System.out.printf("  sketch=%8s\n", Long.toBinaryString(sketch));
                        System.out.printf("e sketch=%8s\n", Long.toBinaryString(expectedSketch));
                    }
                    assertEquals(expectedSketch, sketch);
                }
            }
        }

        //#7 tiles, 4 have set bit flags
        long value = 0b10000000100000001000000010000000000000000000000000000000L;

        sketch = WordLevelParallelism.sketch(value, 7, 7);
        System.out.printf("sketch=%62s\n", Long.toBinaryString(sketch));
        long expected = 0b1111000L;
        //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
        assertEquals(expected, sketch);
    }

    public void est10Compare() {
        int k = 0b1100111;
        //System.out.println("k=" + Integer.toBinaryString(k));
        int kBits = 7;
        int nTiles = 2;
        long tiled1 = WordLevelParallelism.createTiledBitstring1(k, nTiles, kBits);

        long mask1 = WordLevelParallelism.createTiledBitMask1(nTiles, kBits);

        int[] k2 = new int[]{0b0100100,0b1100111};
        long tiled2 = WordLevelParallelism.createTiledBitstring0(k2, kBits);

        long comp = WordLevelParallelism.parallelCompare10(tiled1, tiled2, nTiles, kBits, mask1);

        assertEquals(2, comp);
    }

    public void est6BitTilesSum() {
        /*
        6 bit tiling of 5 bit bitstrings
        value=0b0000100000100000100000100000100000100000100000100000100000100000 #<== 10 tiles comparison, all 10 are set
        kMult=0b0000000000000001000001000001000001000001000001000001000001000001 #<== 9 set bits
        kMask=0b0000000111100000000000000000000000000000000000000000000000000000 #<== 4 bit mask
                     # the number of tiles fits in 4 bits, so set the 4 bits at start of 2nd block
        kShift = 60 - 6 - 1
        sum=(((value * kMult) & kMask) >> kShift) + (value >> 59)
         */
        long comparison = 0b0000100000100000100000100000100000100000100000100000100000100000L;
        //System.out.println("comparison=" + Long.toBinaryString(comparison));

        int nTiles = 10;
        int tileBitLength = 5;

        long sum = WordLevelParallelism.parallelSum(comparison, nTiles, tileBitLength);
        assertEquals(10, sum);

        // change compare to have 7 set bits
        comparison = 0b0000100000100000100000100000100000100000100000000000000000000000L;
        sum = WordLevelParallelism.parallelSum(comparison, nTiles, tileBitLength);
        assertEquals(7, sum);

        // change compare to have 2 set bits
        comparison = 0b0000100000100000000000000000000000000000000000000000000000000000L;
        sum = WordLevelParallelism.parallelSum(comparison, nTiles, tileBitLength);
        assertEquals(2, sum);
    }

    public void est3BitTilesSum() {
        /*
        sumOf for 3 bit tiling of 2 bit bitstrings.  21 tiles (63 bits total).
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
        */
        long comparison = 0b0000100100100100100100100100100100100100100100100100000000000000L;
        //System.out.println("comparison=" + Long.toBinaryString(comparison));

        int nTiles = 20;
        int tileBitLength = 2;

        long sum = WordLevelParallelism.parallelSum(comparison, nTiles, tileBitLength);
        assertEquals(16, sum);
    }
}
