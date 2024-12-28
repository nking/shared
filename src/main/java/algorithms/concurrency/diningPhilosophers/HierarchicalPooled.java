package algorithms.concurrency.diningPhilosophers;

import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class HierarchicalPooled extends Hierarchical {

    public HierarchicalPooled(int n, int thinkTimeMsec) {
        super(n, thinkTimeMsec, HierarchicalPooled.class.getSimpleName());
    }

    @Override
    public void dine(int bitesPer, int timeoutSec) throws InterruptedException {
        initDine(bitesPer);

        ExecutorService exec = null;
        try {
            exec = Executors.newFixedThreadPool(10);
            for (int i = 0; i < n; ++i) {
                exec.submit(eat(i, timeoutSec));
                //exec.execute(eat(i, timeoutSec));
            }
        } finally {
            if (exec != null && !exec.isShutdown()) {
                //exec. awaitTermination((n/2)*timeoutSec, TimeUnit. SECONDS);
                exec.shutdownNow();
            }
        }

        closeDine();
    }
}
