package algorithms.connected;

import algorithms.disjointSets.UnionFind;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * given a set of points, finds the connected among them and places them into
 * groups. note that connected here means adjacent to one another and adjacent
 * is defined by the default "4 neighbor" offsets, but can be overridden to use
 * all 8 neighbors.
 *
 * The runtime complexity is O(N_points).  The resulting pixel coordinates use
 * pixelIdx = (row * data[0].length) + col.
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinder implements IConnectedValuesGroupFinder {

    protected UnionFind uf;

    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;

    protected int minimumNumberInCluster = 3;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    /**
     *
     */
    protected boolean debug = false;
    
    /**
     *
     */
    protected TIntSet excludeValues = null;

    /**
     *
     */
    public ConnectedValuesGroupFinder() {

        this.log = Logger.getLogger(this.getClass().getName());
    }

    /**
     *
     @param setDebugToTrue
     */
    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    /**
     *
     @param n
     */
    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }

    /**
     *
     */
    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    /**
     *
     @param values
     */
    public void setValuesToExclude(TIntSet values) {
        this.excludeValues = values;
    }

    /**
     * find the groups of connected points in pixIdxs where connected means is
     * adjacent to another point in the group and having this pixelValue. The
     * adjacency by default is using the 4 neighbor pattern search unless the
     * user has set that to 8 neighbors. 
     * 
     *The runtime complexity is O(N_points).  The resulting pixel coordinates use
     * pixelIdx = (row * data[0].length) + col.
     *
     @param data array of data.  all rows must be same length.
     @return 
     */
    public List<TLongSet> findGroups(int[][] data) {

        uf = new UnionFind(data.length * data[0].length);

        findClustersInData(data);

        List<TLongSet> groupList = extractClusters();
        
        return groupList;
    }

    /**
     *
     @param data
     */
    protected void findClustersInData(int[][] data) {

        //int w = data.length;
        //int h = data[0].length;
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = new int[]{-1,  0, 1, 0};
            dys = new int[]{ 0, -1, 0, 1};
        } else {
            /*
            for 8 neighbor, can use 4 offsets instead of 8 if visiting all pix
            to avoid repeating calculations
            
             2  *  *  *       (2,0) 1:2,1:1,2:1
             1  *  *  +       (1,1) 0:1,0:0,1:0,2:0,2:1,2:2,1:2,0:2
             0  *  *  *       (2,1) 1:1,1:0,2:0,3:0,3:1,3:2,2:2,1:2
                0  1  2             X: 1:1,1:0,2:0, 1:2
                                    use: +1,-1  +1,0  +1,+1  0:1
             */
            dxs = new int[]{1, 1, 1, 0};
            dys = new int[]{-1, 0, 1, 1};
        }

        int width = data[0].length;
        
        for (int uX = 0; uX < data.length; ++uX) {
            for (int uY = 0; uY < data[uX].length; ++uY) {
                
                int uPixIdx = (uX * width) + uY;
                
                int uValue = data[uX][uY];
                
                if (excludeValues != null && excludeValues.contains(uValue)) {
                    continue;
                }
                
                for (int k = 0; k < dxs.length; ++k) {
                    
                    int vX = uX + dxs[k];
                    int vY = uY + dys[k];
                    
                    if (vX < 0 || vY < 0 || vX >= data.length || vY >= data[0].length) {
                        continue;
                    }
                    
                    int vValue = data[vX][vY];
                    
                    if (excludeValues != null && excludeValues.contains(vValue)) {
                        continue;
                    }
                    
                    if (uValue != vValue) {
                        continue;
                    }
                    
                    int vPixIdx =(vX * width) + vY;

                    if (uf.find(uPixIdx) != uf.find(vPixIdx)) {
                        uf.union(uPixIdx, vPixIdx);
                    }
                }
            }
        }
    }

    protected List<TLongSet> extractClusters() {

        // key = repr node index, value = set of pixels w/ repr
        Map<Integer, Set<Integer>> reprMap = uf.getComponents();

        // rewrite the above into a list
        List<TLongSet> groups = new ArrayList<TLongSet>();

        for (Map.Entry<Integer, Set<Integer>> entry : reprMap.entrySet()) {
            if (entry.getValue().size() >= minimumNumberInCluster) {
                TLongSet set = new TLongHashSet();
                groups.add(set);
                for (int idx : entry.getValue()) {
                    set.add(idx);
                }
            }
        }

        log.finest("number of groups =" + groups.size());
        
        return groups;
    }

}