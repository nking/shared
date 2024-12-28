package algorithms.concurrency.diningPhilosophers;

import junit.framework.TestCase;

public class DinersTest extends TestCase {

    /*
    implementations of a few ways to solve the dining philosophers
    who are sitting at a round table and have one fork in between them
    to share.

    the dining philosophers problem asks how to avoid deadlock while
    enabling concurrency.

    a hierarchical method is implemented.

    and a message passing distributed algorithm is present .

    and I made an even and odd method so that all even numbered chairs eat then all odd,
    etc. with a correction for case when n is odd to alternate first and last
    chairs eating.  It's a toy model that is close to sequential execution, mostly present for runtime comparison.
    The method only uses 2 threads.  It could be adapted to use only 1 thread.
    On a platform that will run the process on only 1 CPU, using fewer threads
    is closer to sequential non-threaded code.  The runtime is dominated by
    the construction of the threads.  And that thread construction cost stays the
    same as n increases.
    The method is faster than hierarchical when n >= 20.

    Hierarchical has the problem of using a thread per fork, and so quickly needs
    alot of memory if a threadpool is not used.

     */

    public void test0() throws InterruptedException {

        int n = 5;
        int bites = 3;
        int thinkTimeMSec = 10;
        int timeoutSec = 1;
        AbstractDiningPhilosophers dp;
        int nRuns = 10;

        // for my  platform,
        // at n~20, Hierarchical begins to take longer than EvenOdd
        // see summary above w/ caveat that garbage collection times are in there too

        for (int i = 0; i < nRuns; ++i) {
            dp = new EvenOdd(n, thinkTimeMSec);
            dp.dine(bites, timeoutSec);

            dp = new Hierarchical(n, thinkTimeMSec);
            dp.dine(bites, timeoutSec);

            dp = new HierarchicalPooled(n, thinkTimeMSec);
            dp.dine(bites, timeoutSec);

            n *= 2;
        }

        System.out.printf("end test0\n");
        System.out.flush();
    }

}
