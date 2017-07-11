package algorithms.util;

import gnu.trove.map.hash.TLongObjectHashMap;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.util.ArrayList;
import java.util.List;

/**
 * A class to estimate the amount of memory an object takes.
 * The user should know if it is to be on the heap or
 * stack (method local variables) when comparing the results to 
 * available memory.
 * 
 * Options for other data-types will be added as needed.
 * 
 * NOTE: haven't considered whether jvm option  UseCompressedOops is used.
 * 
 * NOTE: suggested other tools for determinging object size are
 *    java's Instrumentation.getObjectSize()
 *    and third party library JOL
 *        http://openjdk.java.net/projects/code-tools/jol/
 * 
 * @author nichole
 */
public class ObjectSpaceEstimator {
    
    /*
    memory usage of an object:
    overhead of each object + memory required for type padded to 8 byte values

    overhead = ref to class, gc state info, + synchronization info
    OVERHEAD = 16 bytes
    
    http://users.elis.ugent.be/~leeckhou/papers/SPE06.pdf

    ----------------
    Table III. Java types and their sizes in number of bits when used in the 
    heap (‘ﬁeld size’ column) and when used on the 
    stack (‘size on stack’ column).
                                         32-bit platform          64-bit platform
                                         Field   Size on          Field   Size on
    Java types                           size    stack            size    stack
    boolean                              32       32              32      64
    byte                                 32       32              32      64
    char                                 32       32              32      64
    short                                32       32              32      64
    int                                  32       32              32      64
    ﬂoat                                 32       32              32      64
    reference                            32       32              64      64
    array reference                      32       32              32      32
    returnAddress                        32       32              64      64
    long                                 64       64              64      128
    double                               64       64              64      128

    
    JVM HEAP: 
        -- shared among all virtual machine threads.  
           holds instance variables, static fields, array elements
        -- strings are stored on the heap
    JVM METHOD AREA: [THIS is logically part of the HEAP, 
        but small vms may not gc it.]
        -- shared among all virtual machine threads.  
           holds structures for compiled code:  
               runtime constant pool
               field data
               method data
               code for methods and constructors
    
    Thread Stack: 
        -- local variables, references to objects, and primitives
        Frame:  
            -- new one is created each time a method is invoked and destroyed 
               when it's completed.
               (holds data and partial results.)
            -- has a LIFO operand stack (max depth is determined at compile time)
               -- the operand stack is empty when frame is created 
               -- jvm loads constants, local variable values, or fields onto operand stack
        -- default size dpends on architecture and impl.  usually is 512K or 1024K
    */

    private int nBoolean = 0;
    private int nByte = 0;
    private int nChar = 0;
    private int nShort = 0;
    private int nInt = 0;
    private int nFloat = 0;
    private int nObjRefs = 0;
    private int nArrayRefs = 0;
    private int nLong = 0;
    private int nDouble = 0;
    private int nReturnAddress = 0;
    
    /**
     * Strings have special handling.  They are always on the heap and they
     * may be pooled (re-used) for some JVMs.
     */
    private int nStrings = 0;
    private int maxStringLength = 0;
     
    // sizes in bytes as [heap 32 bit, stack 32 bit, heap 54 bit, stack 64 bit]
    //private final static int objOverheadSz = 16;
    
    private final static int[] word3264 = new int[]{4, 4, 4, 8};
    
    private static String arch = System.getProperty("sun.arch.data.model");
    private static boolean is32Bit = ((arch != null) && arch.equals("64")) ? false : true;
    
    private final static int[] booleanSz = word3264;
    private final static int[] byteSz = word3264;
    private final static int[] charSz = word3264;
    private final static int[] shortSz = word3264;
    private final static int[] intSz  = word3264;
    private final static int[] floatSz = word3264;
    private final static int[] refSz = new int[]{4, 4, 8, 8};
    private final static int[] arrayRefSz = new int[]{4, 4, 8, 8};//= new int[]{4, 4, 4, 4};
    private final static int[] returnAddressSz = new int[]{4, 4, 8, 8};
    private final static int[] longSz = new int[]{8, 8, 8, 16};
    private final static int[] doubleSz = new int[]{8, 8, 8, 16};

    /**
     * @param nBoolean  the number of boolean primitives
     */
    public void setNBooleanFields(int nBoolean) {
        this.nBoolean = nBoolean;
    }
    
    public static long getWordSize() {
        if (is32Bit) {
            return 4;
        } else {
            return 8;
        }
    }
    
    public static long getObjectOverhead() {
        if (is32Bit) {
            return 8;
        } else {
            return 16;
        }
    }

    public static long getReferenceSize() {
        if (is32Bit) {
            return 4;
        } else {
            return 6;
        }
    }
    
    /**
     * @param nByte  the number of byte primitives
     */
    public void setNByteFields(int nByte) {
        this.nByte = nByte;
    }

    /**
     * @param nChar the number of char primitives to set
     */
    public void setNCharFields(int nChar) {
        this.nChar = nChar;
    }

    /**
     * @param nShort  the number of short primitives
     */
    public void setNShortFields(int nShort) {
        this.nShort = nShort;
    }

    /**
     * @param nInt  the number of int primitives
     */
    public void setNIntFields(int nInt) {
        this.nInt = nInt;
    }

    /**
     * @param nFloat  the number of float primitives
     */
    public void setNFloatFields(int nFloat) {
        this.nFloat = nFloat;
    }

    /**
     * @param nObjRefs  the number of object references
     */
    public void setNObjRefsFields(int nObjRefs) {
        this.nObjRefs = nObjRefs;
    }

    /**
     * @param nArrayRefs  the number of array references to set
     */
    public void setNArrayRefsFields(int nArrayRefs) {
        this.nArrayRefs = nArrayRefs;
    }

    /**
     * @param nLong  the number of long primitives
     */
    public void setNLongFields(int nLong) {
        this.nLong = nLong;
    }

    /**
     * @param nDouble the number of double primitives
     */
    public void setNDoubleFields(int nDouble) {
        this.nDouble = nDouble;
    }

    /**
     * @param nReturnAddress the nReturnAddress to set
     */
    public void setNReturnAddress(int nReturnAddress) {
        this.nReturnAddress = nReturnAddress;
    }
    
    public void setStrings(int numberOfStrings, int maximumNumberOfLetters) {
        nStrings = numberOfStrings;
        maxStringLength = maximumNumberOfLetters;
    }
    
    /**
     * A rough estimate of an instance of ArrayList without objects in it.
     * 
     * @return 
     */
    public static long estimateArrayList() {
        
        int index;
        long overhead;
        if (is32Bit) {
            index = 0;
            overhead = 8;
        } else {
            index = 2;
            overhead = 16;
        }
        long pad = getWordSize();
        
        long total = overhead;
        
        total += longSz[index];
        total += 5 * intSz[index];
        total += 3*arrayRefSz[index];
        
        return total;
        
        /*
        final int nLists = 10000;
        final int nObjects = 25;
         
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;

        List[] lists = new List[nLists];
        for (int i = 0; i < nLists/2; ++i) {
            lists[i] = new ArrayList();
            for (int j = 0; j < nObjects; ++j) {
                lists[i].add(new Object());
            }
        }
        
        long heapUsage1 = mbean.getHeapMemoryUsage().getUsed();
        long avail1 = totalMemory - heapUsage1;
        long used1 = heapUsage1 - heapUsage;
        
        for (int i = nLists/2; i < nLists; ++i) {
            for (int j = 0; j < nObjects; ++j) {
                lists[i].add(new Object());
            }
        }
       
        long heapUsage2 = mbean.getHeapMemoryUsage().getUsed();
        long avail2 = totalMemory - heapUsage2;
        long used2 = heapUsage2 - heapUsage;
    
        //System.out.format("%d lists with no objects size =%d\n", nLists, used1);
        //System.out.format("list with %d objects size =%d\n", nObjects, used2);
    
        long nObjectsSize = 
            nObjects * ObjectSpaceEstimator.getReferenceSize();
            //+ nObjects * ObjectSpaceEstimator.getObjectOverhead();    
        
        long diff = used2 - used1;
        diff /= (nLists/2);
        
        // assume that the List and the fields which are not the backing
        // arrays are small in size compared to the nObjects
        diff /= nObjectsSize;
        
        //System.out.println("list after subtracting nObjects from difference =" 
        //    + diff);
        
        return diff;
        */
    }
    
    /**
     * a rough estimate of an instance of TLongObjectHashMap without objects
     * in it.
    
     * @return 
     */
    public static long estimateTLongObjectHashMap() {
    
        int index;
        long overhead;
        if (is32Bit) {
            index = 0;
            overhead = 8;
        } else {
            index = 2;
            overhead = 16;
        }
        long pad = getWordSize();
        
        long total = overhead;
        
        total += 4 * longSz[index];
        total += booleanSz[index];
        total += refSz[index];
        total += 2*arrayRefSz[index];
        
        return total;
        
        /*
        final int nMaps = 10000;
        final int nObjects = 25;
         
        long totalMemory = Runtime.getRuntime().totalMemory();
        MemoryMXBean mbean = ManagementFactory.getMemoryMXBean();
        
        long heapUsage = mbean.getHeapMemoryUsage().getUsed();
        long avail = totalMemory - heapUsage;
        
        TLongObjectHashMap[] maps = new TLongObjectHashMap[nMaps];
        for (int i = 1; i < nMaps/2; ++i) {
            maps[i] = new TLongObjectHashMap();
            for (int j = 0; j < nObjects; ++j) {
                maps[i].put(j, new Object());
            }
        }
        
        long heapUsage1 = mbean.getHeapMemoryUsage().getUsed();
        long avail1 = totalMemory - heapUsage1;
        long used1 = heapUsage1 - heapUsage;
        
        for (int i = nMaps/2; i < nMaps; ++i) {
            maps[i] = new TLongObjectHashMap();
            for (int j = 0; j < nObjects; ++j) {
                maps[i].put(j, new Object());
            }
        }
       
        long heapUsage2 = mbean.getHeapMemoryUsage().getUsed();
        long avail2 = totalMemory - heapUsage2;
        long used2 = heapUsage2 - heapUsage;
    
        //System.out.format("%d maps with no objects size =%d\n", nMaps, used1);
        
        //System.out.format("map with %d objects size =%d\n", nObjects, used2);
        
        long nObjectsSize = 
            nObjects * ObjectSpaceEstimator.getReferenceSize();
        
        long diff = used2 - used1;
        diff /= (nMaps/2);
        
        // assume that the List and the fields which are not the backing
        // arrays are small in size compared to the nObjects
        diff /= (2*nObjects);
        
        //System.out.println("map after subtracting nObjects from difference =" 
        //    + diff);
        
        return diff;
        */
    }
   
    /**
     * estimate the size of an object in bytes for the given settings and for
     * placement on the heap.
     * 
     * Note that the string estimates are the maximum size in memory they could 
     * occupy and do not include consideration for pooling.
     * 
     * @return total size in bytes for the object placed on the heap.
     */
    public long estimateSizeOnHeap() {
        return estimateSize(true);
    }
    
    /**
     * estimate the size of an object in bytes for the given settings and for
     * placement on the stack (variables specific to a method frame, that is,
     * local variables).
     * 
     * Note that the string estimates will not be in the stack estimate as they
     * are purely heap objects.
     * 
     * @return total size in bytes for the object places on the stack.
     */
    public long estimateSizeOnStack() {
        return estimateSize(false);
    }
    
    private long estimateSize(boolean calcForHeap) {
        
        int index;
        long overhead;
        if (is32Bit) {
            if (calcForHeap) {
                index = 0;
            } else {
                index = 1;
            }
            overhead = 8;
        } else {
            if (calcForHeap) {
                index = 2;
            } else {
                index = 3;
            }
            overhead = 16;
        }
        long pad = getWordSize();
        
        long total = 0;

        // add fields
        total += nBoolean * booleanSz[index];
        total += nByte * byteSz[index];
        total += nChar * charSz[index];
        total += nShort * shortSz[index];
        total += nInt * intSz[index];
        total += nFloat * floatSz[index];
        total += nObjRefs * refSz[index];
        total += nArrayRefs * arrayRefSz[index];
        total += nLong * longSz[index];
        total += nDouble * doubleSz[index];
        total += nReturnAddress * returnAddressSz[index];
        
        if (calcForHeap) {
            total += nStrings * estimateAStringSize(maxStringLength);
        }
        
        // add object overhead
        total += overhead;

        // pad up to 8 byte boundary
        long padding = total % pad;
        total += padding;
        
        return total;
    }
   
    public static long estimateLongSize() {
        
        int idx;
        if (is32Bit) {
            idx = 0;
        } else {
            idx = 2;
        }
        
        return longSz[idx];
    }

    public static long estimateAStringSize(int maxNumberOfLetters) {
        
        // strings are always on the heap
        int idx;
        if (is32Bit) {
            idx = 0;
        } else {
            idx = 2;
        }
        
        long total = arrayRefSz[idx];
        total += maxNumberOfLetters * charSz[idx];
        total += 3*intSz[idx];
        
        long pad = total % getWordSize();
        
        total += pad;
        
        return total;
    }
    
}
