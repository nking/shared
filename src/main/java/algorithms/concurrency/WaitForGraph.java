package algorithms.concurrency;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * NOTE: this class is not thread-safe.  Its a toy model.
 */
public class WaitForGraph {
    // key = process ID, value = set of resource IDs
    Map<Integer, Set<Integer>> procResWaitMap = new HashMap<>();
    // key = resource ID, value = process ID... is a 1-1 mapping reverse of procResWaitMap
    Map<Integer, Integer> resProcAllocMap = new HashMap<>();

    // for reverse lookup

    // a graph path is procResWaitMap process ID -> resource ID -> resProcAllocMap resourceID -> processID
    //    etc.

    /**
     * add a resource allocation to the graph (internally, also removes it from process wait map).
     * @param processID
     * @param resourceID
     */
    public void addAlloc(int processID, int resourceID) {
        //int currProcID = resProcAllocMap.get(resourceID);
        resProcAllocMap.put(resourceID, processID);
        if (procResWaitMap.containsKey(processID)) {
            procResWaitMap.get(processID).remove(resourceID);
        }
    }
    protected void addRequest(int processID, int resourceID) {
        procResWaitMap.putIfAbsent(processID, new HashSet<>());
        procResWaitMap.get(processID).add(resourceID);
    }

    public boolean addRequestIfCannotDeadlock(int processID, int resourceID) {
        if (cannotDeadlock(processID, resourceID)) {
            addRequest(processID, resourceID);
            return true;
        }
        return false;
    }

    /**
     * if requesting resourceID for process processID cannot create a deadlock, return true.
     * Note that if the resource is already held by processID, this method will return false.
     * @param processID
     * @param resourceID
     * @return
     */
    public boolean cannotDeadlock(int processID, int resourceID) {

        Set<Integer> visitedP = new HashSet<>();
        visitedP.add(processID);

        Set<Integer> cycle = new HashSet<Integer>();
        cycle.add(processID);

        boolean hasCycle = hasCycle(resourceID, visitedP, cycle);
        return !hasCycle;
    }

    protected boolean hasCycle(int resourceID, Set<Integer> visitedP, Set<Integer> cycle) {

        if (!resProcAllocMap.containsKey(resourceID)) {
            return false;
        }
        int processID = resProcAllocMap.get(resourceID);

        if (cycle.contains(processID)) {
            return true;
        }

        if (visitedP.contains(processID) || !procResWaitMap.containsKey(processID)) {
            return false;
        }

        for (int resID : procResWaitMap.get(processID)) {
            Set<Integer> cycle2 = new HashSet<>(cycle);
            if (hasCycle(resID, visitedP, cycle2)) {
                return true;
            }
        }

        return false;
    }
}
