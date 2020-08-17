package algorithms.util;

import algorithms.misc.Misc0;
import gnu.trove.list.TIntList;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.list.array.TLongArrayList;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInput;
import java.io.ObjectInputStream;
import java.io.ObjectOutput;
import java.io.ObjectOutputStream;
import java.util.Map;
import java.util.Random;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class NodeMapTest extends TestCase {
    
    public NodeMapTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        Random rng = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rng.setSeed(seed);
        
        int n = 100;
        
        NodeMap map = new NodeMap();
    
        TLongList keys = new TLongArrayList();
        TLongList v0s = new TLongArrayList();
        TLongList v1s = new TLongArrayList();
        TLongList v2s = new TLongArrayList();
        TIntList v3s = new TIntArrayList();
        TIntList v4s = new TIntArrayList();
        TIntList v5s = new TIntArrayList();
        TIntList v6s = new TIntArrayList();

        populateData(rng, n, keys, v0s, v1s, v2s, v3s, v4s, v5s, v6s);
        
        for (int i = 0; i < n; ++i) {
            
            long key = keys.get(i);
            long val0 = v0s.get(i);
            long val1 = v1s.get(i);
            long val2 = v2s.get(i);
            int val3 = v3s.get(i);
            int val4 = v4s.get(i);
            int val5 = v5s.get(i);
          
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
                    map.updateParent(key, val0);
                    map.updateLeft(key, val1);
                    map.updateRight(key, val2);
                    v0s.set(i, val0);
                    v1s.set(i, val1);
                    v2s.set(i, val2);
                }
                
            }
            
            assertTrue(map.contains(key));
            assertTrue(map.containsKey(key));
            assertEquals(val0, map.getParent(key));
            assertEquals(val1, map.getLeft(key));
            assertEquals(val2, map.getRight(key));
            assertEquals(val3, map.getNodeValue(key));
            assertEquals(val4, map.getNodeColor(key));
            assertEquals(val5, map.getNodeSize(key));
            
            assertEquals(i + 1, map.size());
            
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

    public void test1() throws IOException {
        
        Random rng = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        System.out.println("SEED=" + seed);
        rng.setSeed(seed);
        
        int n = 100;
        
        TLongList keys = new TLongArrayList();
        TLongList v0s = new TLongArrayList();
        TLongList v1s = new TLongArrayList();
        TLongList v2s = new TLongArrayList();
        TIntList v3s = new TIntArrayList();
        TIntList v4s = new TIntArrayList();
        TIntList v5s = new TIntArrayList();
        TIntList v6s = new TIntArrayList();

        populateData(rng, n, keys, v0s, v1s, v2s, v3s, v4s, v5s, v6s);
        
        NodeMap map = new NodeMap(keys.toArray(new long[n]), 
            v0s.toArray(new long[n]), v1s.toArray(new long[n]), 
            v2s.toArray(new long[n]), 
            v3s.toArray(new int[n]), v4s.toArray(new int[n]), 
            v5s.toArray(new int[n]));
    
        for (int i = 0; i < n; ++i) {
            
            long key = keys.get(i);
            long val0 = v0s.get(i);
            long val1 = v1s.get(i);
            long val2 = v2s.get(i);
            int val3 = v3s.get(i);
            int val4 = v4s.get(i);
            int val5 = v5s.get(i);
          
            assertTrue(map.contains(key));
            
            if ((map.size() % 5) == 0) {
                val0 = Long.MIN_VALUE;
                val1 = Long.MAX_VALUE;
                val2 = Long.MIN_VALUE;
                map.updateParent(key, val0);
                map.updateLeft(key, val1);
                map.updateRight(key, val2);
                v0s.set(i, val0);
                v1s.set(i, val1);
                v2s.set(i, val2);
                
                val3++;
                val4++;
                val5++;
                map.updateNodeValue(key, val3);
                map.updateNodeColor(key, val4);
                map.updateNodeSize(key, val5);
            
                v3s.set(i, val3);
                v4s.set(i, val4);
                v5s.set(i, val5);
            } else {
                assertTrue(map.parentIsSet(key));
                assertTrue(map.leftIsSet(key));
                assertTrue(map.rightIsSet(key));
            }
            
            assertTrue(map.contains(key));
            assertTrue(map.containsKey(key));
            assertEquals(val0, map.getParent(key));
            assertEquals(val1, map.getLeft(key));
            assertEquals(val2, map.getRight(key));
            assertEquals(val3, map.getNodeValue(key));
            assertEquals(val4, map.getNodeColor(key));
            assertEquals(val5, map.getNodeSize(key));
            
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
        
        map.clear();
        assertTrue(map.isEmpty());
        
        map = new NodeMap(keys.toArray(new long[n]), 
            v0s.toArray(new long[n]), v1s.toArray(new long[n]), 
            v2s.toArray(new long[n]), 
            v3s.toArray(new int[n]), v4s.toArray(new int[n]), 
            v5s.toArray(new int[n]));
        
        NodeMap map2 = new NodeMap(map);
        
        assertTrue(map.equals(map2));
        
        map2.clear();
        
        map2.putAll(map);
        
        assertTrue(map.equals(map2));
        
        map2.clear();
        
        ByteArrayOutputStream bout = null;
        ObjectOutput oout = null;
        ObjectInput oin = null;
        ByteArrayInputStream bin = null;
        try {
            bout = new ByteArrayOutputStream();
            oout = new ObjectOutputStream(bout);
            
            map.writeExternal(oout);
            
            bin = new ByteArrayInputStream(bout.toByteArray());
            oin = new ObjectInputStream(bin);
            
            map2.readExternal(oin);
            
            assertTrue(map.equals(map2));
            
        } catch (Exception e) {
            fail(e.getMessage());
        } finally {
            if (oout != null) {
                oout.close();
            }
            if (bout != null) {
                bout.close();
            }
            if (oin != null) {
                oin.close();
            }
            if (bin != null) {
                bin.close();
            }
        }
        
        long[] keys1 = map.keys();
        assertNotNull(keys1);
        assertEquals(map.size(), keys1.length);
        
        TLongSet keySet2 = map2.keySet();
        for (int i = 0; i < keys1.length; ++i) {
            assertTrue(keySet2.contains(keys1[i]));
        }
        
        keys1 = map.keys(new long[map.size()]);
        assertNotNull(keys1);
        assertEquals(map.size(), keys1.length);
        for (int i = 0; i < keys1.length; ++i) {
            assertTrue(keySet2.contains(keys1[i]));
        }
        
        assertNotNull(map.toString());
        
        long missingKey = rng.nextLong();
        while (keySet2.contains(missingKey)) {
            missingKey = rng.nextLong();
        }
        
        assertTrue(map.putIfAbsent(missingKey, 1, 1, 1));
        
        map2 = new NodeMap(map.size());
                        
        map2.putAll(map);
        
        assertTrue(map.equals(map2));
    }
    
    private void populateData(Random rng, int n,
        TLongList keys, TLongList v0s, TLongList v1s, 
        TLongList v2s, TIntList v3s, TIntList v4s, TIntList v5s, TIntList v6s) {
        
        while (keys.size() < n) {
            
            long key = rng.nextLong();
            while (keys.contains(key) || key == Long.MIN_VALUE) {
                key = rng.nextLong();
            }
            int val3 = rng.nextInt();
            while (v3s.contains(val3)) {
                val3 = rng.nextInt();
            }
            
            long val0 = rng.nextLong();
            while (val0 == Long.MIN_VALUE) {
                val0 = rng.nextLong();
            }
            long val1 = rng.nextLong();
            while (val1 == Long.MIN_VALUE) {
                val1 = rng.nextLong();
            }
            long val2 = rng.nextLong();
            while (val2 == Long.MIN_VALUE) {
                val2 = rng.nextLong();
            }
            
            int val4 = rng.nextBoolean() ? 1 : 0;
            int val5 = rng.nextInt(n);
            
            
            keys.add(key);
            v0s.add(val0);
            v1s.add(val1);
            v2s.add(val2);
            v3s.add(val3);
            v4s.add(val4);
            v5s.add(val5);
            
        }
        
    }
    
}
