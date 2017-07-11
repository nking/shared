package algorithms.util;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.Map.Entry;
import java.util.Properties;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ObjectSpaceEstimatorTest extends TestCase {
    
    public ObjectSpaceEstimatorTest(String testName) {
        super(testName);
    }
    
    private class Vars {
        long a = 1;
        int b = 1;
        char c = 1;
        byte d = 1;
        float f = 1;
        short h = 1;
        boolean j = true;
        double g = 1;
        String s = "s";
        int[] arrayRef = new int[1];
    }
  
    public void test1() throws InterruptedException {
        
        /*
        test creation of a large object.
        
        object with array of 1000 primitive longs
               and  array of 1000 primitive ints
        */
        
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        
        final int nL = 1<<16;
        final int nI = 1<<16;
        
        ObjectSpaceEstimator est = new ObjectSpaceEstimator();
        est.setNArrayRefsFields(6);
        est.setNLongFields(1);
        est.setNIntFields(1);
        est.setNCharFields(1);
        est.setNByteFields(1);
        est.setNFloatFields(1);
        est.setNDoubleFields(1);
        est.setNShortFields(1);
        est.setNBooleanFields(1);
        // since strings are interned and re-used, will not include it except 
        //  as single addition
        //est.setNStrings(1, 1);
        
        long heapEstimate = nL * est.estimateSizeOnHeap() +
            ObjectSpaceEstimator.getArrayReferenceSize()
            + ObjectSpaceEstimator.estimateAStringSize(1);
       
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long nonHeapUsage = mbean.getNonHeapMemoryUsage().getUsed();

        Vars[] a = new Vars[nL];
        for (int i = 0; i < nL; ++i) {
            a[i] = new Vars();
        }
       
        long heapUsage2 = mbean.getHeapMemoryUsage().getUsed();
        long nonHeapUsage2 = mbean.getNonHeapMemoryUsage().getUsed();
        long used = heapUsage2 - heapUsage;
        long nhUsed = nonHeapUsage2 - nonHeapUsage;
    
        int code = a.hashCode();
        
        long eps = 8 * nL;
     
        System.out.println("heap used=" + used + " est=" + 
            heapEstimate);
        System.out.println("non-heap Used used=" + nhUsed);
        
        assertTrue(Math.abs(used - heapEstimate) <= eps);
        
        System.out.println("list size w/o items=" + 
            ObjectSpaceEstimator.estimateArrayList() 
            + " TLongObjectHashMap=" +
            ObjectSpaceEstimator.estimateTLongObjectHashMap());
    
    }
  
}
