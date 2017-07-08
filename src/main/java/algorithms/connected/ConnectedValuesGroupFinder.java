package algorithms.connected;

import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TIntObjectIterator;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * given a set of points, finds the connected among them and places them into
 * groups. note that connected here means adjacent to one another and adjacent
 * is defined by the default "4 neighbor" offsets, but can be overridden to use
 * all 8 neighbors.
 *
 * The runtime complexity is essentially O(N_points).
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinder implements IConnectedValuesGroupFinder {

    private final DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();

    // key = pixIdx, value = disjoint set node with key pixIdx
    private TIntObjectMap<DisjointSet2Node<Integer>> pixNodes = null;

    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;

    protected int minimumNumberInCluster = 3;

    protected Logger log = Logger.getLogger(this.getClass().getName());

    protected boolean debug = false;
    
    protected TIntSet excludeValues = null;

    public ConnectedValuesGroupFinder() {

        this.log = Logger.getLogger(this.getClass().getName());
    }

    public void setDebug(boolean setDebugToTrue) {
        this.debug = setDebugToTrue;
    }

    public void setMinimumNumberInCluster(int n) {
        this.minimumNumberInCluster = n;
    }

    public void setToUse8Neighbors() {
        use4Neighbors = false;
    }
    
    public void setValuesToExclude(TIntSet values) {
        this.excludeValues = values;
    }

    /**
     * find the groups of connected points in pixIdxs where connected means is
     * adjacent to another point in the group and having this pixelValue. The
     * adjacency by default is using the 4 neighbor pattern search unless the
     * user has set that to 8 neighbors. The runtime complexity is essentially
     * O(pixIdxs.size()).
     *
     * @param data
     */
    public List<TIntSet> findGroups(int[][] data) {

        initMap(data);

        findClustersIterative(data);

        List<TIntSet> groupList = prune();
        
        return groupList;
    }

    protected void findClustersIterative(int[][] data) {

        int w = data.length;
        int h = data[0].length;
        
        int[] dxs;
        int[] dys;
        if (use4Neighbors) {
            dxs = new int[]{-1,  0, 1, 0};
            dys = new int[]{ 0, -1, 0, 1};
        } else {
            /*
            for 8 neighbor, can use 4 offsets instead of 8 if visiting all pix
            
             2  *  *  *       (2,0) 1:2,1:1,2:1
             1  *  *  +       (1,1) 0:1,0:0,1:0,2:0,2:1,2:2,1:2,0:2
             0  *  *  *       (2,1) 1:1,1:0,2:0,3:0,3:1,3:2,2:2,1:2
                0  1  2             X: 1:1,1:0,2:0, 1:2
                                    use: +1,-1  +1,0  +1,+1  0:1
             */
            dxs = new int[]{1, 1, 1, 0};
            dys = new int[]{-1, 0, 1, 1};
        }

        PixelHelper ph = new PixelHelper();
        
        for (int uX = 0; uX < data.length; ++uX) {
            for (int uY = 0; uY < data[uX].length; ++uY) {
                
                int uPixIdx = ph.toPixelIndex(uX, uY, w);
                
                int uValue = data[uX][uY];
                
                if (excludeValues != null && excludeValues.contains(uValue)) {
                    continue;
                }
                
                for (int k = 0; k < dxs.length; ++k) {
                    
                    int vX = uX + dxs[k];
                    int vY = uY + dys[k];
                    
                    if (vX < 0 || vY < 0 || vX >= w || vY >= h) {
                        continue;
                    }
                    
                    int vValue = data[vX][vY];
                    
                    if (excludeValues != null && excludeValues.contains(vValue)) {
                        continue;
                    }
                    
                    if (uValue != vValue) {
                        continue;
                    }
                    
                    int vPixIdx = ph.toPixelIndex(vX, vY, w);
                    
                    processPair(uPixIdx, vPixIdx);
                }
            }
        }
    }

    protected void processPair(int uPoint, int vPoint) {

        DisjointSet2Node<Integer> uNode = pixNodes.get(uPoint);
        DisjointSet2Node<Integer> uParentNode
            = disjointSetHelper.findSet(uNode);
        assert(uParentNode != null);

        int uGroupId = uParentNode.getMember().intValue();

        DisjointSet2Node<Integer> vNode = pixNodes.get(vPoint);
        DisjointSet2Node<Integer> vParentNode
            = disjointSetHelper.findSet(vNode);
        assert(vParentNode != null);

        int vGroupId = vParentNode.getMember().intValue();

        DisjointSet2Node<Integer> merged =
            disjointSetHelper.union(uParentNode, vParentNode);

        pixNodes.put(uGroupId, merged);

        pixNodes.put(vGroupId, merged);
    }

    private List<TIntSet> prune() {

        // key = repr node index, value = set of pixels w/ repr
        TIntObjectMap<TIntSet> map = new TIntObjectHashMap<TIntSet>();

        TIntObjectIterator<DisjointSet2Node<Integer>> iter =
            pixNodes.iterator();
        for (int i = 0; i < pixNodes.size(); ++i) {

            iter.advance();

            int pixIdx = iter.key();
            DisjointSet2Node<Integer> node = iter.value();

            DisjointSet2Node<Integer> repr = disjointSetHelper.findSet(node);

            int reprIdx = repr.getMember().intValue();

            TIntSet set = map.get(reprIdx);
            if (set == null) {
                set = new TIntHashSet();
                map.put(reprIdx, set);
            }
            set.add(pixIdx);
        }

        log.finest("number of groups before prune=" + map.size());

        // rewrite the above into a list
        List<TIntSet> groups = new ArrayList<TIntSet>();

        TIntObjectIterator<TIntSet> iter2 = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            iter2.advance();

            TIntSet idxs = iter2.value();

            if (idxs.size() >= minimumNumberInCluster) {
                groups.add(idxs);
            }
        }

        log.finest("number of groups after prune=" + groups.size());
        
        return groups;
    }

    private void initMap(int[][] data) {

        int w = data.length;
             
        pixNodes = new TIntObjectHashMap<DisjointSet2Node<Integer>>();
        
        PixelHelper ph = new PixelHelper();
        
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                
                int pixIdx = ph.toPixelIndex(i, j, w);
            
                DisjointSet2Node<Integer> pNode =
                    disjointSetHelper.makeSet(
                    new DisjointSet2Node<Integer>(Integer.valueOf(pixIdx)));

                pixNodes.put(pixIdx, pNode);
            }
        }
    }
}