package algorithms.misc;

import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.TestCase;

import java.io.IOException;
import java.math.BigInteger;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.*;

/**
 *
 * @author nichole
 */
public class WordLevelParallelismTest extends TestCase {

    public WordLevelParallelismTest(String testName) {
        super(testName);
    }

    public void testHighestBitIn() {

        long value = 0b01000000;
        long h = WordLevelParallelism.highestBitSet(value, 8);
        assertEquals(6L, h);

        h = WordLevelParallelism.highestBitSet(value, 7);
        assertEquals(6L, h);

        //        76543210
        value = 0b00010000;
        h = WordLevelParallelism.highestBitSet(value, 8);
        assertEquals(4L, h);

        h = WordLevelParallelism.highestBitSet(value, 7);
        assertEquals(4L, h);

    }
    public void testSketch() {

        //#7 tiles, 4 have set bit flags
        long value = 0b10000000100000001000000010000000000000000000000000000000L;

        long sketch = WordLevelParallelism.sketch(value, 7, 7);

        long expected = Long.parseLong("1111000", 2);

        //System.out.printf("sketch=%8s\n", Long.toBinaryString(sketch));

        assertEquals(expected, sketch);
    }

    public void test10Compare() {
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

    public void test6BitTilesSum() {
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

    public void test3BitTilesSum() {
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
