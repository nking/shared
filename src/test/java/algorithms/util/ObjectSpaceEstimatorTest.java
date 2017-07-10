package algorithms.util;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ObjectSpaceEstimatorTest extends TestCase {
    
    public ObjectSpaceEstimatorTest(String testName) {
        super(testName);
    }
    
    public void test0() {
        
        /*
        test creation of a large object.
        
        object with array of 1000 primitive longs
               and  array of 1000 primitive ints
        */
        final int nL = 100000;
        final int nI = 200000;
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNArrayRefsFields(2);
        est.setNLongFields(nL);
        est.setNIntFields(nI);
        long heapEstimate = est.estimateSizeOnHeap();
        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        Object obj = new Object() {
            long[] a = new long[nL];
            int[] b = new int[nI];
        };
       
        long heapUsage2 = mbean.getHeapMemoryUsage().getUsed();
        long avail2 = totalMemory - heapUsage2;
        long used = avail - avail2;
    
        System.out.println("avail - avail2=" + used + " est=" + 
            heapEstimate);
        
        // seems to be off by 16 Bytes.
        // is the object overhead larger than 16 Bytes on 64 bit platforms?
        long eps = 16;
        
        assertTrue(Math.abs(used - heapEstimate) <= eps);
        
    }
}
