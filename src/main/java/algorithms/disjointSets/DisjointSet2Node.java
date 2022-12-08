package algorithms.disjointSets;

/**
 * a node for the forest implementation of the disjoint set.

   first implemented in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * 
 * @author nichole
 @param <T> parameter type for the container to hold
 */
public class DisjointSet2Node<T> {

    /**
     * the set member data to be held in a disjoint set
     */
    protected T member = null;
    
    /**
     * pointer to the set representative
     */
    protected DisjointSet2Node<T> parent = null;

    /**
    upper bound of edges between this node and the longest path to it's descendants
    */
    protected int rank = 0;
    
    /**
     *
     */
    protected Object data = null;

    /**
     *
     @param member
     */
    public DisjointSet2Node(T member) {
        this.member = member;
    }
    
    /**
     *
     */
    public DisjointSet2Node() {
    }
   
    /**
     * get the member data
     @return 
     */
    public T getMember() {
        return member;
    }

    /**
     *
     @param member
     */
    public void setMember(T member) {
        this.member = member;
    }
    
    /**
     *
     @param data
     */
    public void setDeta(Object data) {
        this.data = data;
    }
    
    /**
     *
     @return
     */
    public Object getObject() {
        return data;
    }

    /**
     * get the set representative
     @return 
     */
    public DisjointSet2Node<T> getParent() {
        return parent;
    }

    /**
     *
     @param theParent
     */
    public void setParent(DisjointSet2Node<T> theParent) {
        this.parent = theParent;
    }

    /**
     *
     @return
     */
    public int getRank() {
        return rank;
    }

    /**
     *
     @param rank
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    @Override
    public String toString() {
        
        StringBuilder sb = new StringBuilder();
        sb.append("    [parent=");
        if (parent != null) {
            if (parent.equals(this)) {
                sb.append("self");
            } else if (parent.member != null) {
                sb.append(parent.member);
            } else {
                sb.append(parent.hashCode());
            }
        }
        sb.append(", ");
        sb.append("member=");
        if (member != null) {
            sb.append(member.toString());
        }
        sb.append(", ");
        sb.append("rank=").append(Integer.toString(rank)).append("; ");
        
        sb.append("data=");
        if (data != null) {
            sb.append(data.toString());
        }
        sb.append(", ");
        if (parent != null) {
            sb.append("hashcode=").append(parent.hashCode());
            sb.append(", ");
        }
        sb.append("] ");
        
        return sb.toString();
    }
  
}
