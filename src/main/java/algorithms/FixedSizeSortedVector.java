package algorithms;

import java.lang.reflect.Array;
import java.util.Arrays;

/**
 * a class to hold a fixed number of items of type T which are kept in a sorted
 * stated.
 *<pre>
 * runtime complexity is:
 *     O(N * (lg_2(k) + smaller than k))
 * where k is the fixed capacity and N is the number of times add is used.
 *
 * worse case runtime complexity is O(N * (k + lg_2(k)))
 * best case runtime complexity is O(N * (1 + lg_2(k)))
 *</pre>
 * 
 * first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
   * 
 * @author nichole
 @param <T> class type to be held and sorted by this class.  It must implement
 * Comparable.
 */
@SuppressWarnings({"unchecked"})
public class FixedSizeSortedVector<T extends Comparable<T>> {

    /**
     *
     */
    protected T[] a = null;

    /**
     *
     */
    protected final int size;

    /**
     *
     */
    protected int n;

    /**
     *
     */
    protected int availSlot;

    /**
     *
     @param fixedCapacity
     @param classTypeToHold
     */
    public FixedSizeSortedVector(int fixedCapacity, Class<T> classTypeToHold) {

        if (fixedCapacity < 1) {
            throw new IllegalArgumentException(
            "fixedCapacity must be a positive non zero (arg was "
            + fixedCapacity + ")");
        }

        size = fixedCapacity;

        n = 0;

        availSlot = -1;
       
        a = (T[]) Array.newInstance(classTypeToHold, size);

    }

    /**
     * add value to the fixed size sorted list using (T).compareTo to order
     * the items in the internal list.
     *
     * runtime complexity is O(log_2(capacity) + less than capacity).
     *
     @param value value to insert into vector
     @return true if added, else false
     */
    public boolean add(T value) {

        if (value == null) {
            return false;
        }

        if (n < size) {

            if (availSlot == -1) {
                availSlot = n;
            }

            insertIntoOpenSlot(value);

        } else {

            int compareIdx = n - 1;

            if ((n == 1) && (size == 1)) {
                compareIdx = 0;
            }

            int comp = value.compareTo(a[compareIdx]);

            if (comp >= 0) {
                return false;
            }

            // free up the last slot
            availSlot = compareIdx;

            n--;

            // insert value into array at position found by binarySearch
            insertIntoOpenSlot(value);

        }

        return true;
    }

    /**
     * Insert the value into the list while maintaining the sorted state
     * of the list.
     @param value value to insert into vector
     */
    private void insertIntoOpenSlot(T value) {

        int insIdx = Arrays.binarySearch(a, 0, n, value);
        if (insIdx < 0) {
            insIdx *= -1;
            insIdx--;
        }

        if ((availSlot > -1) && (insIdx > availSlot) && (a[availSlot].equals(value))) {
            // this depends upon logic of previous remove setting availSlot
            // to next value.
            boolean b = true;
            for (int i = insIdx; i > availSlot; --i) {
                if (!a[i].equals(value)) {
                    b = false;
                    break;
                }
            }

            if (b) {

                // no need to set value again
                n++;

                availSlot = -1;

                return;
            }
        }

        if (insIdx == availSlot) {

            a[availSlot] = value;

        } else if ((insIdx == (a.length - 1)) && (availSlot == (insIdx - 1))) {

            a[insIdx] = value;

        } else if (insIdx < availSlot) {

            // move all items from insIdx to availSlot down by 1
            for (int i = (availSlot - 1); i >= insIdx; i--) {
                a[i + 1] = a[i];
            }

            a[insIdx] = value;

        } else {

            int end = insIdx - 1;

            if ((availSlot > -1) && (a[insIdx] == value)) {
                while (((end + 1) <= n) && (a[end + 1] == value)) {
                    end++;
                }
            }

            // move items up from availSlot +1 to insIdx - 1
            // then insert value into insIdx - 1
            for (int i = availSlot; i < end; i++) {
                a[i] = a[i + 1];
            }

            a[insIdx] = value;
        }

        n++;

        availSlot = -1;
    }

    /**
     * get the internal array for the sorted list.  note this is not a copy in
     * order to keep the use small, so do not edit it and continue to use
     * the add method.  Note that the returned size may be smaller than
     * capacity if the vector was not completely filled.
     *
     * runtime complexity is O(1)
     *
     @return the internal array.  note that this is not a copy, intentionally.
     */
    public T[] getArray() {

        return a;
    }

    /**
     * return the number of items in the internal array.  if the array is not
     * yet filled, the return will be less than the capacity, else will
     * be the same as the capacity.
     @return number of items in the vector
     */
    public int getNumberOfItems() {
        return n;
    }
    
    /**
     * get the maximum size of the vector, given at instantiation.
     @return the capacity of the vector.
     */
    public int getFixedCapacity() {
        return size;
    }
}
