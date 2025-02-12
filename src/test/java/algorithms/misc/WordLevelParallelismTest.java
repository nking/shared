package algorithms.misc;

import algorithms.combPerm.SubsetChooser;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.Random;

/**
 *
 * @author nichole
 */
public class WordLevelParallelismTest extends TestCase {

    final static Random rand;
    static {
        rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //seed = 377949482646163L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
    }

    public WordLevelParallelismTest(String testName) {
        super(testName);
    }

    public void testParallelCompare() {

        /*
          11100111 11100111 11100111 11100111 11100111 11100111 11100111 11100111
        – 00101001 01011101 01100111 01101010 01101011 01101101 01101110 01111111

          10111110 10001010 10000000 01111101 01111100 01111010 01111001 01101000
          1        1        1        0        0        0        0        0
         */
        int nTiles = 7;
        int bitstringLength = 7;

        long tiledQ = WordLevelParallelism.createTiledBitstring1(0b1100111, nTiles, bitstringLength);
        //                 6         5         4         3         2         1
        //               210987654321098765432109876543210987654321098765432109876543210
        assertEquals(0b11100111111001111110011111100111111001111110011111100111L, tiledQ);

        long tiledKeys = WordLevelParallelism.createTiledBitstring0(
                new int[]{0b1011101, 0b1100111, 0b1101010, 0b1101011, 0b1101101, 0b1101110, 0b1111111},
                bitstringLength);

        //          6         5         4         3         2         1
        //        210987654321098765432109876543210987654321098765432109876543210
        assertEquals(0b01011101_01100111_01101010_01101011_01101101_01101110_01111111L, tiledKeys);

        long rank = WordLevelParallelism.rank(tiledQ, tiledKeys, nTiles, bitstringLength);

        assertEquals(2, rank);

    }

    public void testHighestOneBitIn_MSB() {

        //TODO: finish tests

        //           6         5         4         3         2         1
        //         210987654321098765432109876543210987654321098765432109876543210
        //  1000000100000010000001000000100000010000001000000100000010000001000000
        long t = 0b000001000000000010000000000010000000000010000000000010000000000L;
        //            0b1000000000010000000000010000000000010000000000010000000000

        //           6         5         4         3         2         1
        //         210987654321098765432109876543210987654321098765432109876543210
        //  1000000100000010000001000000100000010000001000000100000010000001000000
        // used blocks in
        //       0b100000000000001000000000000010000001000000000000010000000000000
        // sketch of used blocks in:
        //                                                             0b101011010
        int nTiles = 9;
        int tileBitLength = 7-1;
        long msb = WordLevelParallelism.highestOneBitIn(t, nTiles, tileBitLength);
        assertEquals(57, msb);
    }

    public void testHighestBitSetIn() {
        long value;
        int highBit;
        int r;
        int nb;
        long hb;
        int nRTests = 100;
        long mb;
        for (int b = 63;  b > 2; --b) {

            // test no set bits
            value = 0L;
            hb = WordLevelParallelism.highestBitSetIn(value, b);
            assertTrue(hb < 0);

            // test random combinations of the blocks' flag bits for  [0, b-1] inclusive
            for (int kBits = 1; kBits <= b; ++kBits) {
                // randomly choose kBits indexes in the range [0, b] inclusive
                // repeat nTest times
                for (int m = 0; m < nRTests; ++m) {
                    value = 0L;
                    mb = 0L;
                    // set k bits randomly
                    for (int n = 0; n < kBits; ++n) {
                        r = rand.nextInt(b);
                        value |= (1L << r);
                        mb = Math.max(mb, r);
                    }
                    hb = WordLevelParallelism.highestBitSetIn(value, b);

                    if (mb != hb) {
                        System.out.printf("b=%d, kBits=%d\n", b, kBits);
                        System.out.printf("value=%63s\n", Long.toBinaryString(value));
                        //System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        hb = WordLevelParallelism.highestBitSetIn(value, b);
                        System.out.printf("  hb=%d\n", hb);
                        System.out.printf("e hb=%d\n", mb);
                    }
                    assertEquals(mb, hb);
                }
            }
        }
    }

    public void testUsedBlocksIn_highestBlockSetIn() {

        /*
        for each block size b=[1,8], inclusive
           calc nb = the number of blocks in 63 bits
           # testing:
           no bits set within all nb blocks: result should 0
              set all bits within all nb blocks: result should be nb blocks with highest bit set
              randomly set [0,b-1] bits for each of the nb blocks
         */

        long tiled;
        // b is block size in bits
        // nb is the number of blocks in the tiled bitarray.  nb*b <= 63.  this is nTiles.
        // u is the used blocks which is all zeroes except the high bits of blocks which had any bits set in them.
        int nb;
        long u;
        int bn;
        long expectedU;
        int r;
        int kBits;
        int nRTests = 100;
        long hb;
        long expectedHB;

        for (int b = 8; b > 7; --b) {

            nb = (int) Math.floor(63. / b);

            // test all blocks empty:
            tiled = 0L;
            u = WordLevelParallelism.usedBlocksIn(tiled, b);
            assertEquals(0, u);
            hb = WordLevelParallelism.highestBlockSetIn(tiled, nb, b-1);
            assertTrue(hb < 0);

            // set all bits within all nb blocks: result should be nb blocks with highest bit set
            tiled = (1L << (nb*b)) - 1L;
            expectedU = 0;
            for (int i = 0; i < nb; ++i) {
                expectedU |= (1L << ((b-1) + i * b));
            }
            u = WordLevelParallelism.usedBlocksIn(tiled, b);
            hb = WordLevelParallelism.highestBlockSetIn(tiled, nb, b-1);
            if (expectedU != u) {
                System.out.printf("(2) b=%d, nb=%d\n", b, nb);
                System.out.printf("tiled=%63s\n", Long.toBinaryString(tiled));
                u = WordLevelParallelism.usedBlocksIn(tiled, b);
                System.out.printf("    u=%63s\n", Long.toBinaryString(u));
                System.out.printf("e   u=%8s\n", Long.toBinaryString(expectedU));
            }
            assertEquals(expectedU, u);
            assertEquals(nb-1, hb);

            // randomly set [0,b-1] bits for each of the nb blocks
            // repeat nRTests times
            for (int m = 0; m < nRTests; ++m) {
                tiled = 0L;
                expectedU = 0L;
                expectedHB = -1;
                for (bn = 0; bn < nb; ++bn) {
                    // set kBits in each of the nb blocks
                    kBits = rand.nextInt(b);
                    for (int k = 0; k < kBits; ++k) {
                        r = rand.nextInt(b);
                        // set bit r in block nb: nb*b + r
                        // 1_______1_______4_______3_______2_______1_______0_______
                        //                         10987654321098765432109876543210
                        tiled |= (1L << (bn*b + r));
                    }
                    if (kBits > 0) {
                        expectedU |= (1L << ((bn*b) + (b-1)));
                        expectedHB = bn;
                    }
                }
                // test that the sketch finds the same set bits
                u = WordLevelParallelism.usedBlocksIn(tiled, b);
                hb = WordLevelParallelism.highestBlockSetIn(tiled, nb, b-1);
                if (expectedU != u) {
                    System.out.printf("(3)b=%d, nb=%d\n", b, nb);
                    System.out.printf("tiled=%63s\n", Long.toBinaryString(tiled));
                    //System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                    u = WordLevelParallelism.usedBlocksIn(tiled, b);
                    System.out.printf("    u=%8s\n", Long.toBinaryString(u));
                    System.out.printf("e   u=%8s\n", Long.toBinaryString(expectedU));
                }
                assertEquals(expectedU, u);
                assertEquals(expectedHB, hb);
            }
        }
    }

    public void testSketch() {

        // the sketch input is a tiled bit array of all 0's except the flag bits.
        //  it extracts the flag bits and concatenates them and returns that.

        int nRTests = 100;

        SubsetChooser chooser;
        long comp;
        // b is block size in bits
        // nb is the number of blocks in the tiled bitarray.  nb*b <= 63.  this is nTiles.
        // tileBitLength is b-1
        int nb;
        long sketch;
        long expectedSketch;
        int[] selectedIndexes;
        // the block size 2 and 1 sketches take a long time because of the test use of subset chooser
        // so will use random tests for those below this block
        for (int b = 8; b > 2; --b) {//8

            nb = (int) Math.floor(63. / b);

            // test all combinations of the blocks' flag bits for  [0, nb-1] inclusive
            for (int k = 0; k < nb; ++k) {
                //System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                comp = 0L;
                if (k == 0) {
                    sketch = WordLevelParallelism.sketch(comp, nb, b-1);
                    //System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                    //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
                    assertEquals(0b0, sketch);
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
                        System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                        System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        sketch = WordLevelParallelism.sketch(comp, nb, b-1);
                        System.out.printf("  sketch=%8s\n", Long.toBinaryString(sketch));
                        System.out.printf("e sketch=%8s\n", Long.toBinaryString(expectedSketch));
                    }
                    assertEquals(expectedSketch, sketch);
                }
            }
        }

        // for block sizes of 2 and 1
        int r;
        for (int b = 2; b > 0; --b) {
            nb = (int) Math.floor(63. / b);
            // test random combinations of the blocks' flag bits for  [0, nb-1] inclusive
            for (int k = 0; k < nb; ++k) {
                //System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                comp = 0L;
                if (k == 0) {
                    sketch = WordLevelParallelism.sketch(comp, nb, b - 1);
                    //System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                    //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
                    assertEquals(0b0, sketch);
                    continue;
                }
                // randomly choose k indexes in the range [0, nb-1] inclusive
                // repeat nTest times
                for (int m = 0; m < nRTests; ++m) {
                    comp = 0L;
                    expectedSketch = 0L;
                    // set k bits randomly
                    for (int n = 0; n < k; ++n) {
                        r = rand.nextInt(nb);
                        comp |= (1L << ((b - 1) + r * b));
                        expectedSketch |= (1L << r);
                    }
                    // test that the sketch finds the same set bits
                    sketch = WordLevelParallelism.sketch(comp, nb, b - 1);

                    if (expectedSketch != sketch) {
                        System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                        System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                        //System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        sketch = WordLevelParallelism.sketch(comp, nb, b - 1);
                        System.out.printf("  sketch=%8s\n", Long.toBinaryString(sketch));
                        System.out.printf("e sketch=%8s\n", Long.toBinaryString(expectedSketch));
                    }
                    assertEquals(expectedSketch, sketch);
                }
            }
        }
    }

    public void testParallelSum() {

        int nRTests = 100;

        SubsetChooser chooser;
        long comp;
        // b is block size in bits
        // nb is the number of blocks in the tiled bitarray.  nb*b <= 63.  this is nTiles.
        // tileBitLength is b-1
        int nb;
        long sum;
        long expectedSum;
        int[] selectedIndexes;
        // the block size 2 and 1 sketches take a long time because of the test use of subset chooser
        // so will use random tests for those below this block
        int[] bs = new int[]{32, 24, 16, 8,7,6,5,4,3};
        for (int b : bs) {

            nb = (int) Math.floor(63. / b);

            // test all combinations of the blocks' flag bits for  [0, nb-1] inclusive
            for (int k = 0; k < nb; ++k) {
                //System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                comp = 0L;
                if (k == 0) {
                    sum = WordLevelParallelism.parallelSum(comp, nb, b-1);
                    //System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                    //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
                    assertEquals(0L, sum);
                    continue;
                }
                expectedSum = k;
                selectedIndexes = new int[k];
                chooser = new SubsetChooser(nb-1, k);
                while (chooser.getNextSubset(selectedIndexes) != -1) {
                    comp = 0L;
                    // set the high bit of each block in selectedIndexes
                    for (int i = 0; i < selectedIndexes.length; ++i) {
                        //         6         5         4         3         2         1
                        //       210987654321098765432109876543210987654321098765432109876543210
                        //              1_______1_______1_______1_______1_______1_______1_______
                        // block 0, b=8 =>   1_______  bit 7
                        // block 1, b=8 =>   1_______0_______ bit 15

                        comp |= (1L << ((b-1) + selectedIndexes[i] * b));
                    }
                    sum = WordLevelParallelism.parallelSum(comp, nb, b-1);

                    if (expectedSum != sum) {
                        System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                        System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                        System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        System.out.printf("  sum=%d\n", sum);
                        System.out.printf("e sum=%d\n", expectedSum);
                    }
                    assertEquals(expectedSum, sum);
                }
            }
        }

        // for block sizes of 2 and 1
        int r;
        TIntSet set = new TIntHashSet();
        for (int b = 2; b > 1; --b) {
            nb = (int) Math.floor(63. / b);
            // test random combinations of the blocks' flag bits for  [0, nb-1] inclusive
            for (int k = 0; k < nb; ++k) {
                //System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                comp = 0L;
                if (k == 0) {
                    sum = WordLevelParallelism.parallelSum(comp, nb, b - 1);
                    //System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                    //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));
                    assertEquals(0L, sum);
                    continue;
                }
                // randomly choose k indexes in the range [0, nb-1] inclusive
                // repeat nTest times
                for (int m = 0; m < nRTests; ++m) {
                    comp = 0L;
                    set.clear();
                    // set k bits randomly
                    for (int n = 0; n < k; ++n) {
                        r = rand.nextInt(nb);
                        set.add(r);
                        comp |= (1L << ((b - 1) + r * b));
                    }
                    expectedSum = set.size();

                    sum = WordLevelParallelism.parallelSum(comp, nb, b - 1);

                    if (expectedSum != sum) {
                        System.out.printf("b=%d, nb=%d, k=%d\n", b, nb, k);
                        System.out.printf("comp=%63s\n", Long.toBinaryString(comp));
                        //System.out.printf("selected=%s\n", Arrays.toString(selectedIndexes));
                        System.out.printf("  sum=%d\n", sum);
                        System.out.printf("e sum=%d\n", expectedSum);
                    }
                    assertEquals(expectedSum, sum);
                }
            }
        }
    }

}
