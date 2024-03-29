package algorithms.bipartite;

import algorithms.bipartite.MinCostUnbalancedAssignment.LeftNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.RightNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SinkNode;
import algorithms.bipartite.MinCostUnbalancedAssignment.SourceNode;
import gnu.trove.map.TIntObjectMap;
import gnu.trove.map.hash.TIntObjectHashMap;

/**
 *
 * @author nichole
 */
public class PathNodes {
    
    private final TIntObjectMap<RightNode> rightNodes;
    
    private final TIntObjectMap<LeftNode> leftNodes;
    
    private final SourceNode sourceNode;
        
    private final SinkNode sinkNode;
    
    /**
     *
     @param nLeft
     @param nRight
     */
    public PathNodes(int nLeft, int nRight) {
        
        rightNodes = new TIntObjectHashMap<RightNode>();
        leftNodes = new TIntObjectHashMap<LeftNode>();
        sourceNode = new SourceNode();
        sinkNode = new SinkNode();
        
        init(nLeft, nRight);
    }
    
    private void init(int nLeft, int nRight) {
     
        for (int i = 0; i < nLeft; ++i) {
            LeftNode node = new LeftNode();
            node.setKey(Long.MAX_VALUE);
            node.index = i;
            leftNodes.put(i, node);
        }
        for (int i = 0; i < nRight; ++i) {
            RightNode node = new RightNode();
            node.setKey(Long.MAX_VALUE);
            node.index = i;
            rightNodes.put(i, node);
        }
        
        sourceNode.setKey(Long.MAX_VALUE);
        sourceNode.index = nLeft;
        
        sinkNode.setKey(Long.MAX_VALUE);
        sinkNode.index = nRight;       
    }
    
    /**
     *
     @param key
     */
    public void resetNodeExceptData(long key) {
     
        for (int i = 0; i < leftNodes.size(); ++i) {
            LeftNode node = leftNodes.get(i);
            
            node.resetExceptData();
            
            node.setKey(key);
        }
        for (int i = 0; i < rightNodes.size(); ++i) {
            RightNode node = rightNodes.get(i);
            node.resetExceptData();
            node.setKey(key);
        }
        
        sourceNode.setKey(key);
        
        sinkNode.setKey(key);
    }

    /**
     @return the rightNodes
     */
    public TIntObjectMap<RightNode> getRightNodes() {
        return rightNodes;
    }

    /**
     @return the leftNodes
     */
    public TIntObjectMap<LeftNode> getLeftNodes() {
        return leftNodes;
    }

    /**
     @return the sourceNode
     */
    public SourceNode getSourceNode() {
        return sourceNode;
    }

    /**
     @return the sinkNode
     */
    public SinkNode getSinkNode() {
        return sinkNode;
    }
}
