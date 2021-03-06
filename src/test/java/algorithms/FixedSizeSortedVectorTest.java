package algorithms;

import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import junit.framework.TestCase;
import static junit.framework.TestCase.assertEquals;

/**
 *
 * @author nichole
 */
public class FixedSizeSortedVectorTest extends TestCase {
    
    public void testAdd() throws Exception {
        
        FixedSizeSortedVector<Integer> sortedVector = 
            new FixedSizeSortedVector<>(4, Integer.class);
        
        assertTrue(sortedVector.add(Integer.valueOf(7)));
        assertTrue(sortedVector.add(Integer.valueOf(6)));
        assertTrue(sortedVector.add(Integer.valueOf(5)));
        assertTrue(sortedVector.add(Integer.valueOf(4)));
        
        Integer[] values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        assertFalse(sortedVector.add(Integer.valueOf(10)));

        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        assertTrue(sortedVector.add(Integer.valueOf(3)));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);

        assertTrue(values[0].intValue() == 3);
        assertTrue(values[1].intValue() == 4);
        assertTrue(values[2].intValue() == 5);
        assertTrue(values[3].intValue() == 6);
        
        //----
        sortedVector = new FixedSizeSortedVector<>(4, Integer.class);
        
        sortedVector.add(Integer.valueOf(6));
        sortedVector.add(Integer.valueOf(4));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 6);
        
        //test whether can have an array with nulls in it and still use
        //Arrays.binarySearch for non-null indexes
   
        Integer[] a = new Integer[4];
        a[0] = Integer.valueOf(6);
        a[1] = Integer.valueOf(7);
        int idx = Arrays.binarySearch(a, 0, 2, Integer.valueOf(6));
        assertTrue(idx == 0);
        
        sortedVector = 
            new FixedSizeSortedVector<>(4, Integer.class);
        
        sortedVector.add(Integer.valueOf(7));
        sortedVector.add(Integer.valueOf(6));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertTrue(sortedVector.getNumberOfItems() == 2);
        
        assertTrue(values[0].intValue() == 6);
        assertTrue(values[1].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(4));
        sortedVector.add(Integer.valueOf(5));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertTrue(sortedVector.getNumberOfItems() == 4);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
        
        sortedVector.add(Integer.valueOf(10));
        
        values = sortedVector.getArray();
        
        assertNotNull(values);
        assertTrue(sortedVector.getNumberOfItems() == 4);
        
        assertTrue(values[0].intValue() == 4);
        assertTrue(values[1].intValue() == 5);
        assertTrue(values[2].intValue() == 6);
        assertTrue(values[3].intValue() == 7);
    }
    
    public void testTypes() throws Exception {
        
        FixedSizeSortedVector<A> vector = new FixedSizeSortedVector<>(2, A.class);
        
        
        B b = new B(1);
        A a = new A(1);
        C c = new C(0);
        
        vector.add(a);
        vector.add(b);
        vector.add(c);
        
        A[] values = vector.getArray();
        
        assertNotNull(values);        
    }
    
    private static class A implements Comparable<A> {
        int v = 0;
        public A(int v) {
            this.v = v;
        }
        @Override
        public int compareTo(A o) {
            if (o == null) {
                return -1;
            }
            return Integer.compare(v, o.v);
        }
    }
    
    private static class B extends A {
        public B(int v) {
            super(v);
        }
    }
    
    private static class C extends A {
        public C(int v) {
            super(v);
        }
    }
    
    public void testRandom() throws Exception {
        
        SecureRandom sr = SecureRandom.getInstance("SHA1PRNG");
        long seed = System.nanoTime();
        //seed=1393807938003554000l;
        sr.setSeed( seed );
        
        int n = 1000;
        int sz = 100;
                
        List<Integer> sorted = new ArrayList<Integer>();
        
        FixedSizeSortedVector<Integer> sVec = new FixedSizeSortedVector<Integer>(sz, Integer.class);
                
        for (int i = 0; i < n; ++i) {
            
            int number = sr.nextInt(sz);
             
            sorted.add(Integer.valueOf(number));
            
            sVec.add(Integer.valueOf(number));
            
            if (i >= sz) {
                
                Collections.sort(sorted);
                sorted.remove(sorted.size() - 1);
                
                // compare contents
                for (int j = 0; j < sorted.size(); ++j) {
                    Integer number2 = sorted.get(j);
                    assertEquals(number2.intValue(), sVec.getArray()[j].intValue());
                }
            }
        }
    }
}
