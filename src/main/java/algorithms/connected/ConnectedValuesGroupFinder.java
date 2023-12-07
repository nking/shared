package algorithms.connected;

import algorithms.disjointSets.DisjointSet2Helper;
import algorithms.disjointSets.DisjointSet2Node;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongObjectIterator;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
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

    protected final DisjointSet2Helper disjointSetHelper = new DisjointSet2Helper();

    // key = pixIdx, value = disjoint set node with key pixIdx
    protected TLongObjectMap<DisjointSet2Node<Long>> pixNodes = null;

    /**
     * uses the 4 neighbor region if true, else the 8-neighbor region
     */
    protected boolean use4Neighbors = true;

    /**
     *
     */
    protected int minimumNumberInCluster = 3;

    /**
     *
     */
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
     * The runtime complexity is essentially O(pixIdxs.size()).
     *
     @param data
     @return 
     */
    public List<TLongSet> findGroups(int[][] data) {

        initMap(data);

        findClustersIterative(data);

        List<TLongSet> groupList = prune();
        
        return groupList;
    }

    /**
     *
     @param data
     */
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

        PixelHelper ph = new PixelHelper();
        
        for (int uX = 0; uX < data.length; ++uX) {
            for (int uY = 0; uY < data[uX].length; ++uY) {
                
                long uPixIdx = ph.toPixelIndex(uX, uY, w);
                
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
                    
                    long vPixIdx = ph.toPixelIndex(vX, vY, w);
                    
                    processPair(uPixIdx, vPixIdx);
                }
            }
        }
    }

    /**
     *
     @param uPoint
     @param vPoint
     */
    protected void processPair(long uPoint, long vPoint) {

        DisjointSet2Node<Long> uReprNode = disjointSetHelper.findSet(pixNodes.get(uPoint));
        assert(uReprNode != null);

        DisjointSet2Node<Long> vReprNode = disjointSetHelper.findSet(pixNodes.get(vPoint));
        assert(vReprNode != null);

        DisjointSet2Node<Long> mergedNode = disjointSetHelper.union(uReprNode, vReprNode);

        pixNodes.put(uPoint, mergedNode);
        pixNodes.put(vPoint, mergedNode);
        pixNodes.put(uReprNode.getMember().intValue(), mergedNode);
        pixNodes.put(vReprNode.getMember().longValue(), mergedNode);
    }

    private List<TLongSet> prune() {

        // key = repr node index, value = set of pixels w/ repr
        TLongObjectMap<TLongSet> map = new TLongObjectHashMap<TLongSet>();

        TLongObjectIterator<DisjointSet2Node<Long>> iter = pixNodes.iterator();
        for (int i = 0; i < pixNodes.size(); ++i) {

            iter.advance();

            long pixIdx = iter.key();
            DisjointSet2Node<Long> node = iter.value();

            DisjointSet2Node<Long> repr = disjointSetHelper.findSet(node);

            long reprIdx = repr.getMember().longValue();

            TLongSet set = map.get(reprIdx);
            if (set == null) {
                set = new TLongHashSet();
                map.put(reprIdx, set);
            }
            set.add(pixIdx);
        }

        log.finest("number of groups before prune=" + map.size());

        // rewrite the above into a list
        List<TLongSet> groups = new ArrayList<TLongSet>();

        TLongObjectIterator<TLongSet> iter2 = map.iterator();
        for (int i = 0; i < map.size(); ++i) {
            iter2.advance();

            TLongSet idxs = iter2.value();

            if (idxs.size() >= minimumNumberInCluster) {
                groups.add(idxs);
            }
        }

        log.finest("number of groups after prune=" + groups.size());
        
        return groups;
    }

    private void initMap(int[][] data) {

        int w = data.length;
             
        pixNodes = new TLongObjectHashMap<DisjointSet2Node<Long>>();
        
        PixelHelper ph = new PixelHelper();
        
        for (int i = 0; i < data.length; ++i) {
            for (int j = 0; j < data[i].length; ++j) {
                
                long pixIdx = ph.toPixelIndex(i, j, w);
            
                DisjointSet2Node<Long> pNode = disjointSetHelper.makeSet(new DisjointSet2Node<Long>(Long.valueOf(pixIdx)));

                pixNodes.put(pixIdx, pNode);
            }
        }
    }
}