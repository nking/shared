package algorithms.concurrency.diningPhilosophers;

public interface IDiningPhilosophers {
    public void dine(int bitesPer, int timeoutSec) throws InterruptedException;
}
