package algorithms;

import algorithms.misc.MiscMath0;

import java.math.BigInteger;
import java.util.concurrent.Semaphore;

/**
 a thread-safe class to permute numbers given in an array.
   <pre>
   The permute code is adapted from 
        from https://en.wikipedia.org/wiki/Heap%27s_algorithm
   which further references:
        Sedgewick, Robert. "a talk on Permutation Generation Algorithms
        http://www.cs.princeton.edu/~rs/talks/perms.pdf
   
   The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
   by Goetz et al.
   </pre>
   <pre>
   Example invocation and an outline of the internal use of the 2 semaphores and 1 AtomicBoolean:
       int[] seq = new int[]{1,2,3};
       int[] perm = new int[3];
       PermutationsWithAwait p = new PermutationsWithAwait(seq);
       p.getNext(perm);
   </pre>

 Note that the class has been made final because it starts a thread within the constructor and
 any subclass constructor would follow the started thread.
 */
public final class PermutationsWithAwait {

    private final Semaphore availableItem, computationLock;
    
    /**
     * the current permutation
     */
    private final int[] x;

    private final BigInteger nPermutations;

    /**
     * a count of the number of permutations calculated so far.
     * the variable is guarded between the semaphores availableItemLock and computationLock.
     */
    private BigInteger nCurrent;

    /**
     * permute the given set of numbers in a thread that waits for the getNext()
     * invocation to calculate the next permutation.

       <pre>
       The permute code is adapted from 
            from https://en.wikipedia.org/wiki/Heap%27s_algorithm
       which further references:
            Sedgewick, Robert. "a talk on Permutation Generation Algorithms
            http://www.cs.princeton.edu/~rs/talks/perms.pdf

       The semaphore model is adapted from Chap 12.1 of "Java Concurrency in Practice"
       by Goetz et al.
       </pre>
       
     * @param seq input sequence of numbers to permute
     * @throws java.lang.InterruptedException thrown if a thread is interrupted
     */
    public PermutationsWithAwait(int[] seq) throws InterruptedException {

        nPermutations = MiscMath0.factorialBigInteger(seq.length);
        nCurrent = BigInteger.ZERO;

        this.computationLock = new Semaphore(1);
        computationLock.acquire();

        this.availableItem = new Semaphore(0);

        int n = seq.length;
        this.x = new int[n];

        //output(A)
        System.arraycopy(seq, 0, x, 0, n);

        availableItem.release();

        Thread thread = new Thread(new Permuter(seq, x));
        thread.start();
    }
    
    /**
     * check whether all permutations have been returned in getNext() invocations.
     * @return true if all permutations have been returned in getNext() invocations,
     * else false if there are more permutations to be returned by the getNext() argument.
     */
    public boolean hasNext() {
        return (nCurrent.compareTo(nPermutations) == -1);
    }

    String[] getCounts() {
        return new String[]{nPermutations.toString(), nCurrent.toString()};
    }
    
    /**
     * get the next permutation of the original set
     * @param out an output array that will hold the results.  note that out.length
     * must be the same length as the set used during construction of this object instance.
     * @return returns true if a value was set into out, otherwise returns false.
     * @throws InterruptedException thrown when the semaphore acquire throws an InterruptedException
     */
    public boolean getNext(int[] out) throws InterruptedException {

        if (out.length != x.length) {
            throw new IllegalArgumentException("out.length must equal original set.length given to constructor");
        }

        if (nCurrent.compareTo(nPermutations) == -1) {

            //availableItem.tryAcquire(1, TimeUnit.SECONDS);
            availableItem.acquire();

            System.arraycopy(x, 0, out, 0, out.length);

            nCurrent = nCurrent.add(BigInteger.ONE);

            computationLock.release();

            return true;
        }

        return false;
    }
    
    private class Permuter implements Runnable {
        private final int[] in;
        private final int[] out;
        Permuter(final int[] in, final int[] out) {
           this.in = in;
           this.out = out;
        }

        @Override
        public void run() {
            final int n = in.length;
            int[] c = new int[n];
            int i = 0;
            int swap;
            while (i < n) {
                if (c[i] < i) {
                    if ((i & 1) != 1) {
                        // i is an even number
                        swap = in[0];
                        in[0] = in[i];
                        in[i] = swap;
                    } else {
                        swap = in[c[i]];
                        in[c[i]] = in[i];
                        in[i] = swap;
                    }

                    try {
                        //Acquires a permit from this semaphore, blocking until one is
                        //available, or the thread is interrupted.
                        computationLock.acquire();
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                        System.err.println("thread interruption: " + ex.getMessage());
                    }

                    // output permutation to instance member x
                    System.arraycopy(in, 0, out, 0, n);

                    //nCurrent = nCurrent.add(BigInteger.ONE);

                    //Releases a permit, increasing the number of available permits by
                    //one to the semaphore.  If any threads are trying to acquire a permit, then one is
                    //selected and given the permit that was just released.  That thread
                    //is (re)enabled for thread scheduling purposes.
                    availableItem.release();

                    //Swap has occurred ending the for-loop. Simulate the increment 
                    //of the for-loop counter
                    c[i] += 1;
                    //Simulate recursive call reaching the base case by bringing the 
                    //pointer to the base case analog in the array
                    i = 0;

                } else {
                    //Calling generate(i+1, A) has ended as the for-loop terminated. 
                    //Reset the state and simulate popping the stack by incrementing the pointer.
                    c[i] = 0;
                    i++;
                }
            }
        }
    }
}
