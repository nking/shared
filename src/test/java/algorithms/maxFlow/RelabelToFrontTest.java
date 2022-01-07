package algorithms.maxFlow;

import algorithms.VertexNode;
import algorithms.util.PairInt;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class RelabelToFrontTest extends TestCase {
    
    public RelabelToFrontTest(String testName) {
        super(testName);
    }
    
    public void testDischarge() {
        
        // Figure 26.10 from Cormen et al. Introduction to Algorithms"
        
        TIntObjectMap<TIntSet> adj = new TIntObjectHashMap<TIntSet>();
        int nV = 5;
        int i;
        for (i = 0; i < nV; ++i) {
            adj.put(i, new TIntHashSet());
        }
        adj.get(0).add(1);
        adj.get(0).add(2);
        adj.get(1).add(2);
        adj.get(1).add(4);
        adj.get(2).add(3);
        adj.get(3).add(1);
        adj.get(3).add(4);
        int srcIdx = 0;
        int sinkIdx = 4;
        
        TObjectDoubleMap<PairInt> cap = new TObjectDoubleHashMap<PairInt>();
        cap.put(new PairInt(0, 1), 12);
        cap.put(new PairInt(0, 2), 14);
        cap.put(new PairInt(1, 2), 5);
        cap.put(new PairInt(1, 4), 16);
        cap.put(new PairInt(2, 3), 8);
        cap.put(new PairInt(3, 1), 7);
        cap.put(new PairInt(3, 4), 10);
        
        RelabelToFront rTF = new RelabelToFront(adj, cap, srcIdx, sinkIdx);
        System.out.println("rTF initialized");
        rTF.print();
        
        assertEquals(-26., rTF.eF[0]);
        assertEquals(12., rTF.eF[1]);
        assertEquals(14., rTF.eF[2]);
        assertEquals(0., rTF.eF[3]);
        assertEquals(0., rTF.eF[4]);
        assertEquals(5, rTF.h[0]);
        assertEquals(0, rTF.h[1]);
        assertEquals(0, rTF.h[2]);
        assertEquals(0, rTF.h[3]);
        assertEquals(0, rTF.h[4]);
        
        // test expects initial L=1,2,3
        assertEquals(1, rTF.ell.peekFirst().vertex);
        assertEquals(2, ((VertexNode)rTF.ell.peekFirst().next).vertex);
        assertEquals(3, ((VertexNode)rTF.ell.peekFirst().next.next).vertex);
        // assert
        
        int count = 0;
        VertexNode uNode = rTF.ell.peekFirst();
        while (uNode != null) {
            
            System.out.printf("count=%d\n", count);
            
            // assert that all nodes befor u in L have excess==0
            VertexNode pNode = (VertexNode) uNode.prev;
            while (pNode != null) {
                assertTrue(Math.abs(rTF.eF[pNode.vertex]) < 1e-7);
                pNode = (VertexNode) pNode.prev;
            }
            
            rTF.dischargeLoop(uNode);
                        
            switch(count) {
                case 0: {
                    //asserting Fig 26.10(b)
                    assertEquals(-26., rTF.eF[0]);
                    assertEquals(0., rTF.eF[1]);
                    assertEquals(19., rTF.eF[2]);
                    assertEquals(0., rTF.eF[3]);
                    assertEquals(7., rTF.eF[4]);
                    PairInt p = new PairInt(0, 1);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(12., rTF.c.get(p));
                    p = new PairInt(0, 2);
                    assertEquals(14., rTF.f.get(p));
                    assertEquals(14., rTF.c.get(p));
                    p = new PairInt(1, 2);
                    assertEquals(5., rTF.f.get(p));
                    assertEquals(5., rTF.c.get(p));
                    p = new PairInt(1, 4);
                    assertEquals(7., rTF.f.get(p));
                    assertEquals(16., rTF.c.get(p));
                    p = new PairInt(2, 3);
                    assertEquals(8., rTF.c.get(p));
                    p = new PairInt(3, 1);
                    assertEquals(7., rTF.c.get(p));
                    p = new PairInt(3, 4);
                    assertEquals(10., rTF.c.get(p));
                    break;
                }
                case 1: {
                    //asserting Fig 26.10(c)
                    assertEquals(-20., rTF.eF[0]);
                    assertEquals(5., rTF.eF[1]);
                    assertEquals(0., rTF.eF[2]);
                    assertEquals(8., rTF.eF[3]);
                    assertEquals(7., rTF.eF[4]);
                    PairInt p = new PairInt(0, 1);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(12., rTF.c.get(p));
                    p = new PairInt(0, 2);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(14., rTF.c.get(p));
                    p = new PairInt(1, 2);
                    assertEquals(5., rTF.c.get(p));
                    p = new PairInt(1, 4);
                    assertEquals(7., rTF.f.get(p));
                    assertEquals(16., rTF.c.get(p));
                    p = new PairInt(2, 3);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(8., rTF.c.get(p));
                    p = new PairInt(3, 1);
                    assertEquals(7., rTF.c.get(p));
                    p = new PairInt(3, 4);
                    assertEquals(10., rTF.c.get(p));
                    break;
                }
                case 2: {
                    //asserting Fig 26.10(d)
                    assertEquals(-20., rTF.eF[0]);
                    assertEquals(0., rTF.eF[1]);
                    assertEquals(0., rTF.eF[2]);
                    assertEquals(8., rTF.eF[3]);
                    assertEquals(12., rTF.eF[4]);
                    PairInt p = new PairInt(0, 1);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(12., rTF.c.get(p));
                    p = new PairInt(0, 2);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(14., rTF.c.get(p));
                    p = new PairInt(1, 2);
                    assertEquals(5., rTF.c.get(p));
                    p = new PairInt(1, 4);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(16., rTF.c.get(p));
                    p = new PairInt(2, 3);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(8., rTF.c.get(p));
                    p = new PairInt(3, 1);
                    assertEquals(7., rTF.c.get(p));
                    p = new PairInt(3, 4);
                    assertEquals(10., rTF.c.get(p));
                    break;
                }
                case 3: {
                    //asserting Fig 26.10(e)
                    assertEquals(-20., rTF.eF[0]);
                    assertEquals(0., rTF.eF[1]);
                    assertEquals(0., rTF.eF[2]);
                    assertEquals(0., rTF.eF[3]);
                    assertEquals(20., rTF.eF[4]);
                    PairInt p = new PairInt(0, 1);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(12., rTF.c.get(p));
                    p = new PairInt(0, 2);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(14., rTF.c.get(p));
                    p = new PairInt(1, 2);
                    assertEquals(5., rTF.c.get(p));
                    p = new PairInt(1, 4);
                    assertEquals(12., rTF.f.get(p));
                    assertEquals(16., rTF.c.get(p));
                    p = new PairInt(2, 3);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(8., rTF.c.get(p));
                    p = new PairInt(3, 1);
                    assertEquals(7., rTF.c.get(p));
                    p = new PairInt(3, 4);
                    assertEquals(8., rTF.f.get(p));
                    assertEquals(10., rTF.c.get(p));
                    break;
                }
                default:
                    break;
            }
            
            uNode = (VertexNode) uNode.next;
            count++;
        }
    }
}
