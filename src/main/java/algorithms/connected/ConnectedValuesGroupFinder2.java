package algorithms.connected;

import algorithms.StackLongLarge;
import algorithms.util.PixelHelper;
import gnu.trove.iterator.TLongObjectIterator;
import gnu.trove.map.TLongLongMap;
import gnu.trove.map.TLongObjectMap;
import gnu.trove.map.hash.TLongLongHashMap;
import gnu.trove.map.hash.TLongObjectHashMap;
import gnu.trove.set.TIntSet;
import gnu.trove.set.TLongSet;
import gnu.trove.set.hash.TLongHashSet;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * This version is similar to ConnectedValuesGroupFinder.java
 * except that uses smaller amount of space at the
 * expense of a few more O(1) iterations when merging groups.
 * 
 * given a set of points, finds the connected among them and places them into
 * groups. note that connected here means adjacent to one another and adjacent
 * is defined by the default "4 neighbor" offsets, but can be overridden to use
 * all 8 neighbors.
 *
 * The runtime complexity is essentially O(N_points).
 *
 * @author nichole
 */
public class ConnectedValuesGroupFinder2 implements IConnectedValuesGroupFinder {

    // key = pixel index, value = key of pixNodes
    private TLongLongMap pixKeyMap = new TLongLongHashMap();
    
    private TLongObjectMap<TLongSet> keySetMap = new TLongObjectHashMap<TLongSet>();
    
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
    public ConnectedValuesGroupFinder2() {

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
     * user has set that to 8 neighbors. The runtime complexity is essentially
     * O(pixIdxs.size()).
     *
     @param data
     @return 
     */
    public List<TLongSet> findGroups(int[][] data) {

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
        
        int n2 = w * h;
        n2 += (int)(Math.log(n2)/Math.log(2));
            
        StackLongLarge stack = new StackLongLarge(n2);
        
        for (int uX = 0; uX < data.length; ++uX) {
            for (int uY = 0; uY < data[uX].length; ++uY) {
                
                long uPixIdx = ph.toPixelIndex(uX, uY, w);
                
                int uValue = data[uX][uY];
                
                if (excludeValues != null && excludeValues.contains(uValue)) {
                    continue;
                }
                
                stack.push(uPixIdx);
            }
        }
        
        TLongSet visited = new TLongHashSet();
        int[] xy = new int[2];
        
        while (!stack.isEmpty()) {
            
            long uIdx = stack.pop();
            if (visited.contains(uIdx)) {
                continue;
            }
            
            ph.toPixelCoords(uIdx, w, xy);
            
            int uX = xy[0];
            int uY = xy[1];
            int uValue = data[uX][uY];
            
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

                long vIdx = ph.toPixelIndex(vX, vY, w);

                processPair(uIdx, vIdx);
                
                if (!visited.contains(vIdx)) {
                    stack.push(vIdx);
                }
            }
            
            visited.add(uIdx);
        }
    }
    
    private long getKeyRecipr(long pixIdx) {
        
        if (!pixKeyMap.containsKey(pixIdx)) {
            return -1;
        }
        long prevKey = pixIdx;
        long key = pixKeyMap.get(pixIdx);
        long tmp;
        while (prevKey != key) {
            tmp = key;
            key = pixKeyMap.get(key);
            pixKeyMap.put(prevKey, key);
            prevKey = tmp;
        }
        
        return key;
    }

    /**
     *
     @param uPoint
     @param vPoint
     */
    protected void processPair(long uPoint, long vPoint) {
        
        if (pixKeyMap.containsKey(uPoint) && pixKeyMap.containsKey(vPoint)) {
            // put all in u
            long uKey = pixKeyMap.get(uPoint);
            long vKey = pixKeyMap.get(vPoint);
            
            //System.out.println("*u=" + uPoint + " *v=" + vPoint
            //    + " orig keys=" + uKey + " " + vKey);
            
            if (uKey == vKey) {
                return;
            }
            uKey = getKeyRecipr(uPoint);
            vKey = getKeyRecipr(vPoint);
            if (uKey == vKey) {
                return;
            }
            //System.out.println("  rec keys=" + uKey + " " + vKey);
            TLongSet uSet = keySetMap.get(uKey);
            TLongSet vSet = keySetMap.remove(vKey);
            uSet.addAll(vSet);
            pixKeyMap.put(vPoint, uKey);
            pixKeyMap.put(vKey, uKey);
            
            //System.out.println("  key=>" + uKey);
            //System.out.println("  getKeyRecipr(uPoint)=" + getKeyRecipr(uPoint));
            //System.out.println("  getKeyRecipr(vPoint)=" + getKeyRecipr(vPoint));
            //System.out.println("  getKeyRecipr(uKey)=" + getKeyRecipr(uKey));
            //System.out.println("  getKeyRecipr(vKey)=" + getKeyRecipr(vKey));
            
        } else if (pixKeyMap.containsKey(uPoint)) {
            long uKey = getKeyRecipr(uPoint);
            TLongSet uSet = keySetMap.get(uKey);
            uSet.add(vPoint);
            pixKeyMap.put(vPoint, uKey);
            
            //System.out.println("*u=" + uPoint + " v=" + vPoint +
            //    " key=" + uKey);
            
        } else if (pixKeyMap.containsKey(vPoint)) {
            long vKey = getKeyRecipr(vPoint);
            TLongSet vSet = keySetMap.get(vKey);
            vSet.add(uPoint);
            pixKeyMap.put(uPoint, vKey);
            
            //System.out.println("u=" + uPoint + " *v=" + vPoint +
            //    " key=" + vKey);
            
        } else {
            TLongSet pixSet = new TLongHashSet();
            pixSet.add(uPoint);
            pixSet.add(vPoint);
            pixKeyMap.put(uPoint, uPoint);
            pixKeyMap.put(vPoint, uPoint);
            keySetMap.put(uPoint, pixSet);
            
            //System.out.println("u=" + uPoint + " v=" + vPoint +
            //    " key=" + uPoint);
        }
        assert(keySetMap.get(pixKeyMap.get(uPoint)) != null);
        assert(keySetMap.get(pixKeyMap.get(vPoint)) != null);
    
    }

    private List<TLongSet> prune() {

        // remove from keySetMap, sets smaller than min limit
        
        List<TLongSet> groups = new ArrayList<TLongSet>();
            
        TLongObjectIterator<TLongSet> iter = keySetMap.iterator();
        
        for (int i = 0; i < keySetMap.size(); ++i) {

            iter.advance();
            
            TLongSet pixSet = iter.value();
            
            if (pixSet.size() < minimumNumberInCluster) {
                continue;
            }

            groups.add(pixSet);
        }
        
        keySetMap = null;
        pixKeyMap = null;

        log.finest("number of groups after prune=" + groups.size());
        
        return groups;
    }

}