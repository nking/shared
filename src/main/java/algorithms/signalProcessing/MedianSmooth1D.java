package algorithms.signalProcessing;

import java.util.Arrays;

/**
  
 first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.
 
 * @author nichole
 */
public class MedianSmooth1D {
    
    /**
     * calculate a running median of a window of size xWindow, yWindow.
     * runtime complexity is
     *     n_rows * ((xWindow * yWindow) + ((n_cols)*lg2(xWindow * yWindow)))
     * so is roughly O(n_pixels * lg_2(window area)) where n_pixels = n_rows * n_cols
     *
     * NOTE: should only be used by a single thread.
     * 
     * NOTE: the border points outside of the window retain their 
     * initial values.
     *
     * @param input
     * @param window
     * @return
     */
    public float[] calculate(float[] input, int window) {

        if (input == null) {
            throw new IllegalArgumentException("input cannot be null");
        }
        if (input.length < window) {
            throw new IllegalArgumentException(
            "input.lenth must be equal to or greater than window");
        }

        int nW = window;

        int xh = window/2;

        //NOTE: to use zero-padding: output = input.createWithDimensions();
        float[] output = Arrays.copyOf(input, input.length);

        SortedVector sVec = new SortedVector(nW);

        // add the first nW to the sorted vector
        for (int i = 0; i < window; ++i) {
            sVec.append(input[i]);
        }

        assert(sVec.n == sVec.a.length);
        assert(sVec.sorted);

        //O(k) + (N)*lg2(k)
        float median;

        for (int i = (window - 1); i < input.length; ++i) {

            //O(1)
            median = sVec.getMedian();

            output[i - xh] = median;

            // remove each item from last column in window
            // and add each item in next column for window,

            if ((i + 1) < input.length) {

                assert(sVec.n == sVec.a.length);

                // remove : O(log_2(k))
                sVec.remove(input[i - window + 1]);

                assert(sVec.n == (sVec.a.length - 1));

                // add : O(log_2(k)) + < O(k)
                sVec.insertIntoOpenSlot(input[i + 1]);

                assert(sVec.n == sVec.a.length);
            }
        }

        return output;
    }
     
    /**
     * a fixed size list that keeps the contents sorted after the capacity is
     * reached.  points are added one at a time and removed one at a time
     * and there are rules to prevent removing when list is not full or
     * adding when list is full.
     */
    static class SortedVector {
        
        protected final float[] a;

        protected int n;

        protected int availSlot;

        protected boolean sorted;

        public SortedVector(int size) {

            a = new float[size];

            n = 0;

            availSlot = -1;

            sorted = false;
        }

        /**
         * append item value onto the end of the list.  Note that if the item
         * is added to the last slot, the list is immediately sorted into
         * ascending numerical order
         * afterwards as a side effect to keep the logic in the other
         * methods consistent.
         * runtime is usually O(1), but if append is used for the last item,
         * there is a sort adding O(N*log_2(N)).
         * For best use, append(v) the first size-1 items and thereafter use
         * insertIntoOpenSlot(v).
         *
         * @param value
         */
        public void append(float value) {

            if (n == (a.length)) {
                throw new IllegalArgumentException(
                    "0) there must be an empty slot in order to append." +
                    " remove and item then try insert again or construct larger list.");
            }

            a[n] = value;

            n++;

            if (n == a.length) {

                Arrays.sort(a);

                sorted = true;
            }
        }

        /**
         * Insert the value into the list while maintaining the sorted state
         * of the list.  Note that if there is not exactly one available slot
         * in the list, an IllegalArgumentException will be thrown.
         * runtime is usually O(log_2(N)) + less than O(N), but once per class lifetime
         * the sort may occur here adding O(N*log_2(N)).
         * @param value
         */
        public void insertIntoOpenSlot(float value) {

            if (n != (a.length - 1)) {
                String err = "1) the method is meant to be used only on a full list." 
                + " a.length=" + a.length + " n=" + n;
                throw new IllegalArgumentException(err);
            }

            if (!sorted) {
                // this can happen if the user used "size - 1" append()s followed
                // by insertIntoOpenSlot.  It's only needed once for lifetime
                // of object.

                if (availSlot != -1) {
                    throw new IllegalArgumentException(
                        "Error in the algorithm... should have been sorted already");
                }

                a[n] = value;

                n++;

                Arrays.sort(a);

                sorted = true;
                
                return;
            }

            int insIdx = Arrays.binarySearch(a, value);
            if (insIdx < 0) {
                insIdx *= -1;
                insIdx--;
            }

            if (insIdx == availSlot) {

                a[availSlot] = value;

            } else if (insIdx < availSlot) {

                // move all items from insIdx to availSlot down by 1
                for (int i = (availSlot - 1); i >= insIdx; i--) {
                    a[i + 1] = a[i];
                }

                a[insIdx] = value;

            } else {

                int end = insIdx - 1;

                // move items up from availSlot +1 to insIdx - 1
                // then insert value into insIdx - 1
                for (int i = availSlot; i < end; i++) {
                    a[i] = a[i + 1];
                }

                a[insIdx - 1] = value;
            }
            n++;
            availSlot = -1;            
        }

        /**
         * remove the item from the full list of items.
         * runtime is O(log_2(N)).
         * NOTE: this could be made O(1) runtime complexity 
         * at the expense
         * of 3 * space complexity.
         * @param value
         */
        public void remove(float value) {

            if (n != a.length) {
                throw new IllegalArgumentException(
                "2) the method is meant to be used only on a full list." 
                + " a.length=" + a.length + " n=" + n);
            }

            int rmIdx = Arrays.binarySearch(a, value);

            if (rmIdx < 0) {
                throw new IllegalArgumentException("could not find item in list");
            }

            availSlot = rmIdx;

            // to keep the list in a state where the next binary search works,
            // set the empty slot value to the proceeding value or max integer.
            if (availSlot == (a.length - 1)) {
                a[availSlot] = Float.POSITIVE_INFINITY;
            } else {
                a[availSlot] = a[availSlot + 1];
            }

            n--;            
        }

        /**
         * get median from the internal array.  Note that this will
         * throw an IllegalArgumentException if the list is not full.
         * runtime is O(1)
         * @return median
         */
         public float getMedian() {

            if (n != a.length) {
                // NOTE: in the use above, this is never invoked unless the
                // list a is full so this exception should never be thrown
                throw new IllegalArgumentException(
                    "3) the method is meant to be used only on a full list." 
                    + " a.length=" + a.length + " n=" + n);
            }

            int midIdx = ((n & 1) == 1) ? n/2 : (n - 1)/2;

            return a[midIdx];
        }
    }
    
}
