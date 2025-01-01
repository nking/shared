package algorithms.concurrency.diningPhilosophers;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.ReentrantLock;

public class Hierarchical extends AbstractDiningPhilosophers {
    protected ReentrantLock[] forkLocks;

    public Hierarchical(int n, int thinkTimeMsec) {
        super(n, thinkTimeMsec, Hierarchical.class.getSimpleName());

        this.forkLocks = new ReentrantLock[n];
        for (int i = 0; i < n; ++i) {
            forkLocks[i] = new ReentrantLock();
        }
    }

    public Hierarchical(int n, int thinkTimeMsec, String label) {
        super(n, thinkTimeMsec, label);

        this.forkLocks = new ReentrantLock[n];
        for (int i = 0; i < n; ++i) {
            forkLocks[i] = new ReentrantLock(true);
        }
    }

    //   0   1   2
    // 2  0    1   2
    protected ReentrantLock getLowerFork(int i) {
        if (i == 0) {
            return forkLocks[0];
        }
        return forkLocks[i - 1];
    }

    protected ReentrantLock getHigherFork(int i) {
        if (i == 0) {
            return forkLocks[n - 1];
        }
        return forkLocks[i];
    }

    public void dine(int bitesPer, int timeoutSec) throws InterruptedException {

        initDine(bitesPer);

        // start all threads.
        // each tries to obtain low number fork, high number fork, eat, release forks. think
        List<Thread> threadList = new ArrayList<>();
        for (int i = 0; i < n; ++i) {
            Thread thr = new Thread(eat(i, timeoutSec));
            threadList.add(thr);
            thr.start();
        }
        for (int i = 0; i < n; ++i) {
            threadList.get(i).join();
        }

        closeDine();
    }//end dine

    protected Runnable eat(int i, int timeoutSec) {
        return new Runnable() {
            @Override
            public void run() {
                while (bites.containsKey(i)) {
                    try {
                        boolean s = getLowerFork(i).tryLock(timeoutSec, TimeUnit.MILLISECONDS);
                        if (!s) continue;
                        s = getHigherFork(i).tryLock(timeoutSec, TimeUnit.MILLISECONDS);
                        if (!s) continue;
                        if (bites.containsKey(i)) {
                            if (bites.get(i) == 1) {
                                bites.remove(i);
                            } else {
                                bites.put(i, bites.get(i) - 1);
                            }
                            //System.out.printf("%d eats\n", i);
                        }
                    } catch (InterruptedException e) {
                        Thread.currentThread().interrupt();
                    } finally {
                        if (getLowerFork(i).isHeldByCurrentThread()) {
                            getLowerFork(i).unlock();
                        }
                        if (getHigherFork(i).isHeldByCurrentThread()) {
                            getHigherFork(i).unlock();
                            try {
                                Thread.sleep(thinkTimeMilliSec);
                            } catch (InterruptedException e) {
                                Thread.currentThread().interrupt();
                            }
                        }
                    }
                }
            }
        };
    }
}
