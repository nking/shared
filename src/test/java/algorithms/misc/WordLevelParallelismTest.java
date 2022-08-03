package algorithms.misc;

import algorithms.util.FormatArray;
import algorithms.util.PolygonAndPointPlotter;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.TestCase;

import java.io.IOException;
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
    
    public void test0() throws NoSuchAlgorithmException, IOException {
        int k = Integer.parseInt("1100111", 2);
        System.out.println("k=" + Integer.toBinaryString(k));
        int kBits = 7;
        int nTiles = 2;
        long tiled1 = WordLevelParallelism.createTiledBitstring1(k, nTiles, kBits);

        long mask1 = WordLevelParallelism.createTiledBitMask1(nTiles, kBits);

        int[] k2 = new int[]{Integer.parseInt("0100100", 2),
            Integer.parseInt("1100111", 2)};
        long tiled2 = WordLevelParallelism.createTiledBitstring0(k2, kBits);

        String tiled2Expected = "0110011100100100";

        long comp = WordLevelParallelism.parallelCompare10(tiled1, tiled2, nTiles, kBits, mask1);
    }
}
