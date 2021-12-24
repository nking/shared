package algorithms;

import java.util.Arrays;
import java.util.concurrent.Semaphore;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicBoolean;

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
   <pre>
   Example invocation and an outline of the internal use of the 2 semaphores and 1 AtomicBoolean:
       int[] seq = new int[]{1,2,3}; int[] perm = new int[3];
       PermutationsWithAwait p = new PermutationsWithAwait(seq);
       p.getNext(perm);
      
       invokes these operations:
         constructor: availableItem = new Semaphore(0);
                      computationLock = new Semaphore(1);
                      computationLock.acquire(); //Acquires a permit, returns immed
                      availableItem.release(); //Releases a permit, incr nAvailPermits by +1
                      finished = new AtomicBoolean(false);
         PermThread:  computationLock.acquire(); // wait for computationLock.release() or thread interruption by another thread
         getNext():   if finished==true, return
                      availableItem.acquire(); // acquire or wait for availableItem.release() or thread interruption by another thread
                      copy data to out var
                      computationLock.release(); // Releases a permit, incr nAvailPermits by +1
         PermThread:  computationLock.acquire() // acquire or wait for computationLock.release() or thread interruption by another thread
                      data computation
                      availableItem.release(); // Releases a permit, incr nAvailPermits by +1
                      if permutations are done, sets finished = true
   </pre>
 */
public class PermutationsWithAwait {
    
    private final Semaphore availableItem, computationLock;
    
    /**
     * the current permutation
     */
    private final int[] x;
    
    /**
     * becomes true when the run-loop has ended for the permuter thread
     */
    private AtomicBoolean finished;
    
    //private final BigInteger nPermutations;
    
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
       
     * @param seq
     * @throws java.lang.InterruptedException
     */
    public PermutationsWithAwait(int[] seq) throws InterruptedException {
                
        int n = seq.length;
        this.x = new int[n];
        
        this.availableItem = new Semaphore(0);
        this.computationLock = new Semaphore(1);
        
        //nPermutations = MiscMath0.factorialBigInteger(n);
        
        computationLock.acquire();
        
        //output(A)
        System.arraycopy(seq, 0, x, 0, n);
        //nCurrent = BigInteger.ONE;
        
        availableItem.release();
        
        finished = new AtomicBoolean(false);
                
        Thread thread = new Thread(new Permuter(seq, finished));
        thread.start();
    }
    
    /**
     * check whether all permutations have been returned in getNext() invocations.
     * @return true if all permutations have been returned in getNext() invocations,
     * else false if there are more permutations to be returned by the getNext() argument.
     */
    public boolean hasNext() {
        return !finished.get();
    }
    
    /**
     * get the next permutation of the original set
     * @param out an output array that will hold the results.  note that out.length
     * must be the same length as the set used during construction of this object instance.
     * @throws InterruptedException 
     */
    public void getNext(int[] out) throws InterruptedException {
        if (out.length != x.length) {
            throw new IllegalArgumentException("out.length must equal original set.length given to constructor");
        }
        
        if (finished.get()) {
            return;
        }
        availableItem.tryAcquire(2, TimeUnit.SECONDS);
        
        System.arraycopy(x, 0, out, 0, out.length);
       
        computationLock.release();        
    }
    
    //TODO: consider using Callable so run can throw an exception
    private class Permuter implements Runnable {
        private final int[] s;
        final AtomicBoolean permDone;
        Permuter(int[] seq, AtomicBoolean permDone) {
           this.s = Arrays.copyOf(seq, seq.length);
           this.permDone = permDone;
        }

        @Override
        public void run() {
            final int n = s.length;
            int[] c = new int[n];
            int i = 0;
            int swap;
            while (i < n) {
                if (c[i] < i) {
                    if ((i & 1) != 1) {
                        // i is even number
                        swap = s[0];
                        s[0] = s[i];
                        s[i] = swap;
                    } else {
                        swap = s[c[i]];
                        s[c[i]] = s[i];
                        s[i] = swap;
                    }
                                        
                    try {
                        //Acquires a permit from this semaphore, blocking until one is
                        //available, or the thread is interrupted.
                        computationLock.acquire();
                    } catch (InterruptedException ex) {
                        Thread.currentThread().interrupt();
                    }
                    
                    // output permutation to instance member x
                    System.arraycopy(s, 0, x, 0, n);
                    
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
            permDone.set(true);
        }
    }
}
