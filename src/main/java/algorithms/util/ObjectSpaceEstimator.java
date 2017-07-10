package algorithms.util;

/**
 * A class to estimate the amount of memory an object takes.
 * The user should know if it is to be on the heap or
 * stack (method local variables) when comparing the results to 
 * available memory.
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
}
