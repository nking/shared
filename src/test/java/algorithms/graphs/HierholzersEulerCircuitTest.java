package algorithms.graphs;

import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class HierholzersEulerCircuitTest extends TestCase {
    
    public HierholzersEulerCircuitTest(String testName) {
        super(testName);
    }
    
     public void test0() {
        // test is from:
        // https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
        // written by
        // https://in.linkedin.com/in/ashutosh-kumar-9527a7105
 
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        g.put(0, new TIntHashSet());
        g.put(1, new TIntHashSet());
        g.put(2, new TIntHashSet());
        g.get(0).add(1);
        g.get(1).add(2);
        g.get(2).add(0);
        
        HierholzersEulerCircuit ec = new HierholzersEulerCircuit();
        int[] circuit = ec.createCircuit(g);
        
        int[] expected = new int[]{0, 1, 2, 0};
        assertTrue(Arrays.equals(expected, circuit));
    }
     
    public void test1() {
        // test is from:
        // https://www.geeksforgeeks.org/hierholzers-algorithm-directed-graph/
        // written by
        // https://in.linkedin.com/in/ashutosh-kumar-9527a7105
 
        TIntObjectMap<TIntSet> g = new TIntObjectHashMap<TIntSet>();
        for (int i = 0; i < 7; ++i) {
            g.put(i, new TIntHashSet());
        }
        g.get(0).add(1);
        g.get(0).add(6);
        g.get(1).add(2);
        g.get(2).add(0);
        g.get(2).add(3);
        g.get(3).add(4);
        g.get(4).add(2);
        g.get(4).add(5);
        g.get(5).add(0);
        g.get(6).add(4);
        
        HierholzersEulerCircuit ec = new HierholzersEulerCircuit();
        int[] circuit = ec.createCircuit(g);
                
        int[] expected = new int[]{0, 6, 4, 5, 0, 1, 2, 3, 4, 2, 0};
        assertTrue(Arrays.equals(expected, circuit));
        
    }
    
}
