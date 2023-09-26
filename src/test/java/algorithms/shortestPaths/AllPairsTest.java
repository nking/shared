package algorithms.shortestPaths;

import algorithms.util.SimpleLinkedListNode;
import gnu.trove.map.TIntIntMap;
import gnu.trove.map.hash.TIntIntHashMap;
import junit.framework.TestCase;

import java.io.IOException;
import java.util.Arrays;

/**
 *
 * @author nichole
 */
public class AllPairsTest extends TestCase {

    public AllPairsTest() {
    }
    
    public void test0() throws Exception {

        StanfordMilesReader reader = new StanfordMilesReader();
        reader.loadFile();
        long[] x = new long[reader.nCities];
        long[] y = new long[reader.nCities];
        reader.fillWithCoordinates(x, y);

    }
}
