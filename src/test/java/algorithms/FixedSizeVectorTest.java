package algorithms;

import java.util.Arrays;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class FixedSizeVectorTest extends TestCase {
    
    public FixedSizeVectorTest(String testName) {
        super(testName);
    }

    public void testAdd() {
        
        FixedSizeVector<Integer> instance = new 
            FixedSizeVector<Integer>(3, Integer.class);
        
        for (int i = 0; i < 7; ++i) {
            instance.add(i);
            if (i < 3) {
                assertEquals(i + 1, instance.size());
            }
        }
        assertEquals(6, instance.get(2).intValue());
        assertEquals(5, instance.get(1).intValue());
        assertEquals(4, instance.get(0).intValue());

        Integer[] a = instance.getArray();
        Arrays.sort(a);
        assertEquals(6, instance.get(2).intValue());
        assertEquals(5, instance.get(1).intValue());
        assertEquals(4, instance.get(0).intValue());
        
        assertEquals(3, instance.size());
        assertEquals(3, instance.capacity);
    }

}
