package algorithms.concurrency.diningPhilosophers;

public class WithMessagePassing extends AbstractDiningPhilosophers {

    public WithMessagePassing(int n, int thinkTimeMsec, String label) {
        super(n, thinkTimeMsec, label);
        throw new UnsupportedOperationException("not yet implemented");
    }

    @Override
    public void dine(int bites, int timeoutSec) {
        throw new UnsupportedOperationException("not yet implemented");
    }

    /*

    NOTE: when implementing, if not using blocking queues,
    consider Thread.onSpinWait() for condition...

        Chandy-Misra message passing algorithm:

        init:
           each process gets a unique int ID.
           the fork between 2 people is given to the lower ID diner
           and marked as dirty.

                  diner   0      1    2
              fork      d  d       d      and marked as 'dirty'

         thinking:
            if receives fork req from neighbor: clean it and gives to neighbor
         preparing to eat:
            any chopsticks that diner doesn't have are requested from neighbor.
            if neighbor asks diner for a fork they have:
                if dirty, clean it and give to neighbor.
                else if clean, puts request in a queue, and keeps fork
         has both forks:
            can begin eating.
            any requests are deferred (stored in queue).
            forks are now dirty.
         immed after eating
            diner cleans forks and gives forks to requesters.
             the priority is lowered here using the state "dirty" for the fork.

         Wait for graph (from wikipedia):
            -- node is a process P_i.
            -- edge from process P_i to P_j implies P_j holds a resource needed by P_i
               hence P_i is waiting on P_j to release its lock on resource.
            -- if P_i is waiting on more than 1 lock, there may be relationships
               between the locks such as disjunctive (OR) or conjunctive (AND).
               In the conjunctive case, a possible deadlock is a cycle(s).
               In the disjunctive case, a possible deadlock is a knot(s) where a knot is
               an inescapable section of a directed graph.

          will refer to fork with number and 'd' for dirty or 'c' for clean.
          f0d is fork 0 dirty.

          for 3 diners and a fair queueing system:
          ----------------------
          init: P0 has f0d, f1d
                P1 has f2d
          ---------------------
          P1 -> req f1c -> P0
          P2 -> req f0c -> P0
          P2 -> req f2c -> P1
          P0 eats
          --------------------
          P0 cleans forks and gives (releases mutex) to P1 and P2.
          P1 has fc1, fc2
          P2 has fc0
          P0 is thinking
          -------------------
          P0  -> req f1c -> P1
          P0 -> req f2c -> P2
          P1 eats
          -------------------
          P1 cleans forks and gives (releases mutex) to P0 and P2.
          P2 has has fc0, f2c
          P0 has ...

          no cycles

          TODO: add mutex details and use P language
         */
}
