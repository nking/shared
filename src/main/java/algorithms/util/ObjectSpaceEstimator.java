package algorithms.util;

/**
 * A class to estimate the amount of memory an object takes.
 * The user should know if it is to be on the heap or
 * stack (method local variables) when comparing the results to 
 * available memory.
 * 
 * Options for other data-types will be added as needed.
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
    JVM METHOD AREA: [THIS is logically part of the HEAP, but in small vms may not gc it.]
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
    
    //TODO: handle string, which is a character array but is ALWAYS on the heap
    //   and for some jvm, strings are pooled objects, re-used.
    
    // sizes in bytes as [heap 32 bit, stack 32 bit, heap 54 bit, stack 64 bit]
    private final static int objOverheadSz = 16;
    
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
    private final static int[] arrayRefSz = new int[]{4, 4, 4, 4};
    private final static int[] returnAddressSz = new int[]{4, 4, 8, 8};
    private final static int[] longSz = new int[]{8, 8, 8, 16};
    private final static int[] doubleSz = new int[]{8, 8, 8, 16};

    /**
     * @param nBoolean  the number of boolean primitives
     */
    public void setNBooleanFields(int nBoolean) {
        this.nBoolean = nBoolean;
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
   
    /**
     * estimate the size of an object in bytes for the given settings and for
     * placement on the heap.
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
     * @return total size in bytes for the object places on the stack.
     */
    public long estimateSizeOnStack() {
        return estimateSize(false);
    }
    
    private long estimateSize(boolean calcForHeap) {
        
        int idx;
        if (is32Bit) {
            if (calcForHeap) {
                idx = 0;
            } else {
                idx = 1;
            }
        } else {
            if (calcForHeap) {
                idx = 2;
            } else {
                idx = 3;
            }
        }
        
        long total = 0;

        // add fields
        total += nBoolean * booleanSz[idx];
        total += nByte * byteSz[idx];
        total += nChar * charSz[idx];
        total += nShort * shortSz[idx];
        total += nInt * intSz[idx];
        total += nFloat * floatSz[idx];
        total += nObjRefs * refSz[idx];
        total += nArrayRefs * arrayRefSz[idx];
        total += nLong * longSz[idx];
        total += nDouble * doubleSz[idx];
        total += nReturnAddress * returnAddressSz[idx];

        // add object overhead
        total += 8;

        // pad up to 8 byte boundary
        long pad = total % 8;
        total += pad;
        
        return total;
    }
}
