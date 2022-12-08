package algorithms;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * a class to hold a fixed number of items of type T.  The latest capacity 
 * number of items are what is present in the vector.
 * 
 *<pre>
 * runtime complexity is:
 *     O(1) for inserts, and gets
 * Note that deletes are not allowed,
 *    but could be implemented in O(N) runtime.
 * 
 * where N is the fixed capacity.
 *
 *</pre>
 * 
 * @author nichole
 @param <T> class type to be held and sorted by this class.  It must implement
 * Comparable.
 */
@SuppressWarnings({"unchecked"})
public class FixedSizeVector<T> {

    /**
     *
     */
    protected T[] a = null;

    /**
     *
     */
    protected final int capacity;

    /**
     *
     */
    protected int n;

    /**
     *
     */
    protected int lastIdx;

    /**
     *
     @param fixedCapacity
     @param classTypeToHold
     */
    public FixedSizeVector(int fixedCapacity, Class<T> classTypeToHold) {

        if (fixedCapacity < 1) {
            throw new IllegalArgumentException(
            "fixedCapacity must be a positive non zero (arg was "
            + fixedCapacity + ")");
        }

        capacity = fixedCapacity;

        n = 0;

        lastIdx = -1;
       
        a = (T[]) Array.newInstance(classTypeToHold, capacity);
    }

    /**
     * add value to the fixed size sorted list using (T).compareTo to order
     * the items in the internal list.
     *
     * runtime complexity is O(log_2(capacity) + less than capacity).
     *
     @param value value to insert into vector
     */
    public void add(T value) {

        if (value == null) {
            return;
        }

        if (n < capacity) {
            n++;
        }
        
        lastIdx++;
        if (lastIdx >= capacity) {
            lastIdx = 0;
        }
        a[lastIdx] = value;
    }
    
    /**
     *
     @param itemIndex
     @return
     */
    public T get(int itemIndex) {
        
        if (itemIndex < 0 || itemIndex >= capacity) {
            throw new IllegalArgumentException("itemIndex is out of bounds: " +
                " capacity=" + capacity + " itemIndex=" + itemIndex);
        }
        if (n < capacity) {
            if (itemIndex > lastIdx) {
                throw new IllegalArgumentException("itemIndex is out of bounds: " +
                " capacity=" + capacity + " and last set index=" + lastIdx
                + " itemIndex=" + itemIndex);
            }
            return a[itemIndex];
        }
        
        /*
         0  1  2  3  
        
        c=4, last=1
        get i=2   idx=last-(c-i-1)   1-1=0
        
        c=4, last=3
        get i=2   idx=last-(c-i-1)   3-1=2
        
        c=4, last=1
        get i=1   idx=last-(c-i-1)   1-2=-1 to +c=3
        */
        
        int idx = lastIdx - (capacity - itemIndex - 1);
        if (idx < 0) {
            idx += capacity;
        }
        
        return a[idx];
    }

    /**
     * get a copy of the internal array.
     *
     * runtime complexity is O(N)
     *
     @return the internal array.  note that this is not a copy, intentionally.
     */
    T[] getArray() {

        return Arrays.copyOf(a, n);
    }

    /**
     * return the number of items in the internal array.  if the array is not
     * yet filled, the return will be less than the capacity, else will
     * be the same as the capacity.
     @return number of elements in vector
     */
    public int size() {
        return n;
    }
    
    /**
     * get the maximum size of the vector, given at instantiation.
     @return capacity
     */
    public int getFixedCapacity() {
        return capacity;
    }
}
