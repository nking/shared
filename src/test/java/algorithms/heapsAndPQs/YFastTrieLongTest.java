package algorithms.heapsAndPQs;

import algorithms.misc.Misc0;
import gnu.trove.list.TLongList;
import gnu.trove.list.array.TLongArrayList;
import java.util.Random;
import static junit.framework.Assert.assertEquals;
import static junit.framework.Assert.assertTrue;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class YFastTrieLongTest extends TestCase {
    
    public YFastTrieLongTest() {
    }
    
    public void testMemoryEstimates() {
        long[] est;
        
        long MB = 1024 * 1024;
        
        System.out.println("mem is in MB:");
        
        int[] x = new int[]{
              5000, 1024, 512, 256, 128, 64
        };
        int[] y = new int[]{
              7000, 1024, 512, 256, 128, 64
        };
        
        double rt;
        for (int lenIdx = 0; lenIdx < x.length; ++lenIdx) {
            
            float f = 1.f;
            
            for (int i = 0; i < 2; ++i) {
            
                if (i == 1) {
                    f = 0.1f;
                }
                
                int w = 62;
                
                int n = Math.round(f * x[lenIdx] * y[lenIdx]);
                
                est = YFastTrieLong.estimateSizeOnHeap(n, w);
                
                rt = Math.log(w)/Math.log(2);
                
                System.out.format(
                    "width=%4d, height=%4d n=%10d mem=%10d:%10d  w=%2d  rt=%3d\n",
                    x[lenIdx], y[lenIdx], n, est[1]/MB, est[0]/MB,
                    w, (int)Math.round(rt));
            
                w = (int)Math.round(Math.log(x[lenIdx] * y[lenIdx])/Math.log(2));
                
                est = YFastTrieLong.estimateSizeOnHeap(n, w);
                
                rt = Math.log(w)/Math.log(2);
                
                System.out.format(
                    "width=%4d, height=%4d n=%10d mem=%10d:%10d  w=%2d  rt=%3d\n",
                    x[lenIdx], y[lenIdx], n, est[1]/MB, est[0]/MB,
                    w, (int)Math.round(rt));
            }
        }

    }

    public void testKeyOperations0() throws Exception {
    
        //System.out.println("testKeyOperations");
        
        int n = 100;
        
        Random rand = Misc0.getSecureRandom();
        long seed = System.nanoTime();
        //seed = 1500230008845L;
        System.out.println("SEED=" + seed);
        rand.setSeed(seed);
        
        YFastTrieLong bt = new YFastTrieLong(8);
        
        int count = 0;
        
        TLongList nodes = new TLongArrayList(2*n);
        
        for (int i = 0; i < n/2; ++i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i += skip;
                continue;
            }
            nodes.add(i);
            bt.add(i);
            //System.out.println("add " + i);
            assertEquals(i, bt.find(i));
            assertEquals(nodes.get(0), bt.minimum());
            assertEquals(i, bt.maximum());
            count++;
        }
        for (int i = (n - 1); i >= (n/2); --i) {
            if (rand.nextBoolean()) {
                int skip = rand.nextInt(5);
                i -= skip;
                continue;
            }
            nodes.add(i);
            bt.add(i);
            //System.out.println("add " + i);
            assertEquals(i, bt.find(i));
            assertEquals(nodes.get(0), bt.minimum());
            count++;
        }
        nodes.sort();
        int n2 = bt.size();
        assertEquals(count, n2);
        assertEquals(count, nodes.size());
        
        //bt.debugPrint();
        
        for (int nIter = 0; nIter < 3; ++nIter){
            n2 = bt.size();
            assertEquals(n2, nodes.size());

            for (int i = 0; i < n2 - 1; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.find(idx);
                assertTrue(foundIndex > -1);

                long expected = nodes.get(i + 1);
                //System.out.println("\n* " + idx + " expected next=" + expected);

                long next = bt.successor(idx);
                assertTrue(next != -1);
                //System.out.println(idx + "   next=" + next);
                assertEquals(expected, next);

            }

            for (int i = 1; i < n2; ++i) {
                long idx = nodes.get(i);
                long foundIndex = bt.find(idx);
                assertTrue(foundIndex > -1);

                assertEquals(idx, foundIndex);
                
                long expected = nodes.get(i - 1);
                //System.out.println("\n* " + idx + " expected prev=" + expected);

                long prev = bt.predecessor(idx);
                assertTrue(prev != -1);
                //System.out.println(idx + "   prev=" + prev);
                
                assertEquals(expected, prev);

            }
        
            if ((nIter & 1) == 1) {
                //randomly remove some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = rand.nextInt(nodes.size());
                    long v = nodes.get(idx);
                    assertEquals(v, bt.find(v));
                    bt.remove(v);
                    assertEquals(-1, bt.find(v));
                    nodes.removeAt(idx);
                    assertEquals(nodes.size(), bt.size());
                }
            } else {
                //randomly add some nodes
                for (int i = 0; i < n2/4; ++i) {
                    int idx = nodes.size() + rand.nextInt(2*nodes.size());
                    // the nodes contains is linear search, so could use a temp
                    // set if this test gets large one day
                    if (!nodes.contains(idx)) {
                        bt.add(idx);
                        nodes.add(idx);
                        assertEquals(idx, bt.find(idx));
                    }
                }
            }
            
            n2 = nodes.size();
            nodes.sort();
                               
            long max = nodes.get(nodes.size() - 1);
            assertEquals(max, bt.find(max));
            //System.out.println("will delete max=" + max + " from this tree:");
            //bt.printPreOrderTraversal();
            long max2 = bt.extractMaximum();
            assertEquals(max, max2);
            //System.out.println(" after delete:");
            //bt.printPreOrderTraversal();
            assertEquals(-1, bt.find(max));
            nodes.removeAt(nodes.size() - 1);

            long min = nodes.get(0);
            assertEquals(min, bt.find(min));
            long min2 = bt.extractMinimum();
            assertEquals(min, min2);
            assertEquals(-1, bt.find(min));
            nodes.removeAt(0);
        }        
    }

}
