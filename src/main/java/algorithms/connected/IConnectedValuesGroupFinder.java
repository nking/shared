package algorithms.connected;

import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import java.util.List;

/**
 * given a set of points, finds the connected among them and places them into
 * groups. note that connected here means adjacent to one another and adjacent
 * is defined by the default "4 neighbor" offsets, but can be overridden to use
 * all 8 neighbors.
 *
 *
 * @author nichole
 */
public interface IConnectedValuesGroupFinder {

    public void setDebug(boolean setDebugToTrue);

    public void setMinimumNumberInCluster(int n);

    public void setToUse8Neighbors();
    
    public void setValuesToExclude(TIntSet values);

    /**
     * find the groups of connected points in pixIdxs where connected means is
     * adjacent to another point in the group and having this pixelValue. The
     * adjacency by default is using the 4 neighbor pattern search unless the
     * user has set that to 8 neighbors. The runtime complexity is essentially
     * O(pixIdxs.size()).
     *
     * @param data
     */
    public List<TLongSet> findGroups(int[][] data);

}