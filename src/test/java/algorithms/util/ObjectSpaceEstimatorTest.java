package algorithms.util;

import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
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
    
    public void test0() throws InterruptedException {
        
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
    
        // seems to be off by 16 Bytes.
        // is the object overhead larger than 16 Bytes on 64 bit platforms?
        long eps = 16;
     
        /*
        long nonheapUsage = mbean.getNonHeapMemoryUsage().getUsed();
        
        long frameUsed = loadInMethod(mbean, nonheapUsage, nL, nI);
        
        long nonheapUsage2 = mbean.getNonHeapMemoryUsage().getUsed();
        long usedNH = nonheapUsage2 - nonheapUsage;
        
        System.out.println("heap used=" + used + " est=" + 
            heapEstimate);
                
        System.out.println("NON HEAP used=" + frameUsed + " then after frme=" + usedNH);
        
        long heapUsage3 = mbean.getHeapMemoryUsage().getUsed();
        long avail3 = totalMemory - heapUsage3;
        long used3 = avail - avail3;
        
        System.out.println("heap used 3 =" + used3);
        */
        
        assertTrue(Math.abs(used - heapEstimate) <= eps);
    }
    
    private long loadInMethod(final MemoryMXBean mbean,
        final long nonheapUsage, final int nL, final int nI) throws InterruptedException {
        
        /*
        Each frame contains:
            Local variable array
            Return value
            Operand stack
            function arguments
            Reference to runtime constant pool for class of the current method

        minimum reported non-heap from loading a method frame on 
            64 bit platform is 32 Bytes?  or that is from 2 methods, that is
            including the method invoking this.
        
        http://blog.jamesdbloom.com/images_2013_11_17_17_56/JVM_Internal_Architecture.png
        
        The memory is released after the method returns.
        */
     
        runThreads(1000, 100, mbean, nonheapUsage);
        
        long nonheapUsage2 = mbean.getNonHeapMemoryUsage().getUsed();
        long used = nonheapUsage2  - nonheapUsage;
        
        //the non-heap amount doesn't appear to scale with the number
        //   of threads created.
        //   appears to be staying under 256kB
        
        Thread.sleep(1000); 
        
        return used;
    }

    private void runThreads(int n, final int time, final MemoryMXBean mbean,
        final long nonheapUsage) {

        Thread[] rs = new Thread[n];
        for (int i = 0; i < n; ++i) {
            rs[i] = new Thread(
                new Runnable() {
                @Override
                public void run() {
                    try {
                        //long nonheapUsage2 = mbean.getNonHeapMemoryUsage().getUsed();
                        //long used = nonheapUsage2  - nonheapUsage;
                        //System.out.println("   nh usage=" + used);
                        Thread.sleep(time);
                    } catch (InterruptedException ex) {
                        Logger.getLogger(ObjectSpaceEstimatorTest.class.getName()).log(Level.SEVERE, null, ex);
                    }
                }
            });
        }
        for (int i = 0; i < n; ++i) {
            rs[i].start();
        }
    }
}
