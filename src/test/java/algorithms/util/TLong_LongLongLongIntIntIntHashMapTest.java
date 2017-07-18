package algorithms.util;

import algorithms.misc.Misc0;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TLong_LongLongLongIntIntIntHashMapTest extends TestCase {
    
    public TLong_LongLongLongIntIntIntHashMapTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        Random rng = Misc0.getSecureRandom();
        long seed = System.currentTimeMillis();
        System.out.println("SEED=" + seed);
        rng.setSeed(seed);
        
        int n = 100;
        
        TLong_LongLongLongIntIntIntHashMap map 
            = new TLong_LongLongLongIntIntIntHashMap();
    
        TLongList keys = new TLongArrayList();
        TLongList v0s = new TLongArrayList();
        TLongList v1s = new TLongArrayList();
        TLongList v2s = new TLongArrayList();
        TIntList v3s = new TIntArrayList();
        TIntList v4s = new TIntArrayList();
        TIntList v5s = new TIntArrayList();
        TIntList v6s = new TIntArrayList();
    
        while (keys.size() < n) {
            
            long key = rng.nextLong();
            while (keys.contains(key)) {
                key = rng.nextLong();
            }
            int val3 = rng.nextInt();
            while (v3s.contains(val3)) {
                val3 = rng.nextInt();
            }
            
            //TODO: could randomly choose from keys here to maintain left, right
            // and parent links, but should put that test in a wrapper method
            // for using the map as a tree
            long val0 = rng.nextLong();
            long val1 = rng.nextLong();
            long val2 = rng.nextLong();
            
            int val4 = rng.nextBoolean() ? 1 : 0;
            int val5 = rng.nextInt(n);
            
            assertFalse(map.contains(key));
            
            if ((map.size() & 1) == 1) {
                assertTrue(map.putIfAbsent(key, val0, val1, val2, val3, val4, val5));
                assertFalse(map.putIfAbsent(key, val3, val4, val5));
            } else {
                
                map.put(key, val0, val1, val2, val3, val4, val5);
                
                if ((map.size() % 5) == 0) {
                    val0 = Long.MIN_VALUE;
                    val1 = Long.MAX_VALUE;
                    val2 = Long.MIN_VALUE;
                    map.updateValue0(key, val0);
                    map.updateValue1(key, val1);
                    map.updateValue2(key, val2);
                }
                
            }
            
            assertTrue(map.contains(key));
            assertTrue(map.containsKey(key));
            assertEquals(val0, map.getValue0(key));
            assertEquals(val1, map.getValue1(key));
            assertEquals(val2, map.getValue2(key));
            assertEquals(val3, map.getValue3(key));
            assertEquals(val4, map.getValue4(key));
            assertEquals(val5, map.getValue5(key));
            
            keys.add(key);
            v0s.add(val0);
            v1s.add(val1);
            v2s.add(val2);
            v3s.add(val3);
            v4s.add(val4);
            v5s.add(val5);
            
            assertEquals(keys.size(), map.size());
            
        }
        
        int ndel = rng.nextInt(n/2);
        
        int endSize = map.size() - ndel;
        
        while (map.size() > endSize) {
            int idx = rng.nextInt(map.size());
            long key = keys.get(idx);
            
            assertTrue(map.remove(key));
            assertFalse(map.containsKey(key));
            keys.removeAt(idx);
            v0s.removeAt(idx);
            v1s.removeAt(idx);
            v2s.removeAt(idx);
            v3s.removeAt(idx);
            v4s.removeAt(idx);
            v5s.removeAt(idx);
            
            assertEquals(keys.size(), map.size());
        }
    }
}
