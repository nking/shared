package algorithms.concurrency;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * class for checking a wait for graph of processes waiting on resources.
 *
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
     * @param processId the process id
     * @param resourceId the resource id
     */
    public void addAlloc(int processId, int resourceId) {
        //int currProcID = resProcAllocMap.get(resourceId);
        resProcAllocMap.put(resourceId, processId);
        if (procResWaitMap.containsKey(processId)) {
            procResWaitMap.get(processId).remove(resourceId);
        }
    }
    
    protected void addRequest(int processID, int resourceID) {
        procResWaitMap.putIfAbsent(processID, new HashSet<>());
        procResWaitMap.get(processID).add(resourceID);
    }

    /**
     * add a request to the wait-for-graph if the request cannot deadlock, else return false without
     * adding the request.
     * @param processId
     * @param resourceId
     * @return true if the request was added.
     */
    public boolean addRequestIfCannotDeadlock(int processId, int resourceId) {
        if (cannotDeadlock(processId, resourceId)) {
            addRequest(processId, resourceId);
            return true;
        }
        return false;
    }

    /**
     * determine whether a request for a resource can deadlock.
     * @param processId
     * @param resourceId
     * @return whether the request would lead to a deadlock in the wait-for-graph.
     */
    public boolean cannotDeadlock(int processId, int resourceId) {

        Set<Integer> visitedP = new HashSet<>();
        visitedP.add(processId);

        boolean hasCycle = hasCycle(resourceId, visitedP);
        return !hasCycle;
    }

    /**
     * determine if there is a cycle in the wait for graph.
     * @param resourceId resource if
     * @param visitedP visited map of process ids
     * @return true if a cycle is found.
     */
    protected boolean hasCycle(int resourceId, Set<Integer> visitedP) {

        if (!resProcAllocMap.containsKey(resourceId)) {
            return false;
        }
        int processID = resProcAllocMap.get(resourceId);

        if (visitedP.contains(processID)) {
            return true;
        }

        if (!procResWaitMap.containsKey(processID)) {
            return false;
        }

        for (int resID : procResWaitMap.get(processID)) {
            if (hasCycle(resID, visitedP)) {
                return true;
            }
        }

        return false;
    }
}
