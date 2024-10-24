package algorithms;

import java.util.Comparator;
import java.util.TreeSet;

/**
 * class with methods to find the k largest or smallest values in an array.
 * The runtime complexities are O(N*log(k)) where N is the length of the array and the space
 * complexities are O(k).
 <br/>
 * To make a simple parallel version for very datasets larger than Integer.MAX_VALUE,
 * one could partition the data in j arrays. solve for top k from each j,
 * then solve for the top k among the j result arrays.
 * The total runtime complexity for the simple parallel version would be O((N + j) * log)k)).
 * If have access to all j arrays at the same time and can write to same heap
 * (no large communication costs), one can use a single heap and the same approach to
 * result in runtime complexity O(N*log(k)).
 */
public class KSelect {

    //NOTE: can be made faster by using a YFastTrie and a Map to count the multiplicity of numbers.

    protected static class Node<T> implements Comparable<T>{
        double time;
        T val;
        public Node(T value) {
            if (!(value instanceof Integer) && !(value instanceof Double)) {
                throw new IllegalArgumentException("value must be an interger or double");
            }
            this.time = System.nanoTime();
            this.val = value;
        }

        @Override
        public boolean equals(Object obj) {
            if (!(obj instanceof Node)) {
                return false;
            }
            Node other = (Node)obj;
            if (time != other.time) return false;
            return (val == other.val);
        }

        @Override
        public int hashCode() {
            int hash = Double.hashCode(time);
            hash ^= val.hashCode();
            hash *= Integer.MAX_VALUE;
            return hash;
        }

        /** default is ascending.
         * only implemented for integer and double
         * */
        @Override
        public int compareTo(T o) {
            Node<T> other = (Node)o;
            int c;
            if (val instanceof Integer) {
                c = Integer.compare((int)val, (int) other.val);
            } else {
                c = Double.compare((double)val, (double) other.val);
            }
            if (c == 0) {
                return Double.compare(time, other.time);
            }
            return c;
        }
    }
    protected static class MaxNode<T> extends Node<T>{
        public MaxNode(T value) {
            super(value);
        }
        @Override
        public int compareTo(T o) {
            return -1*super.compareTo(o);
        }
    }

    /**
     * find the top k elements of k, that is, the k largest values in a.
     * The runtime complexity is O(n * log(k)) where n is a.length
     * and the space cmplexity is O(k).
     * @param a
     * @param k
     * @return the k largest elements of k in descending order.
     */
    public static double[] topK(double[] a, int k) {
        TreeSet<MaxNode<Double>> heap = new TreeSet<>();

        for (double v : a) {
            heap.add(new MaxNode<>(v));
            if (heap.size() > k) {
                heap.pollLast();
            }
        }
        double[] out = new double[Math.min(k, heap.size())];
        int iOut = 0;
        for (MaxNode<Double> v : heap) {
            out[iOut++] = v.val;
        }
        return out;
    }

    /**
     * find the bottom k elements of k, that is, the k smallest values in a.
     * The runtime complexity is O(n * log(k)) where n is a.length
     * and the space cmplexity is O(k).
     * @param a
     * @param k
     * @return the k smallest elements of k in ascending order.
     */
    public static double[] bottomK(double[] a, int k) {
        TreeSet<Node<Double>> heap = new TreeSet<>();

        for (double v : a) {
            heap.add(new Node<>(v));
            if (heap.size() > k) {
                heap.pollLast();
            }
        }
        double[] out = new double[Math.min(k, heap.size())];
        int iOut = 0;
        for (Node<Double> v : heap) {
            out[iOut++] = v.val;
        }
        return out;
    }

    /**
     * find the top k elements of k, that is, the k largest values in a.
     The runtime complexity is O(n * log(k)) where n is a.length
     * and the space cmplexity is O(k).
     * @param a
     * @param k
     * @return the k largest elements of k in descending order.
     */
    public static int[] topK(int[] a, int k) {
        // heap needs add, remove and max or min.
        // heap has size k, all values <= current min.
        // need a multiset.

        // could use fibonacci heap or YFastTrie wrapper or a TreeSet for values.
        // the TreeSet would need to create an object that is different for having same
        // value from a (object implements equals and hashcode and uses time, e.g.).
        // The comparator can be in the TreeSet rather than the object in order
        // to allow other methods to use the object

        TreeSet<MaxNode<Integer>> heap = new TreeSet<>();

        for (int v : a) {
            heap.add(new MaxNode<>(v));
            if (heap.size() > k) {
                heap.pollLast();
            }
        }
        int[] out = new int[Math.min(k, heap.size())];
        int iOut = 0;
        for (MaxNode<Integer> v : heap) {
            out[iOut++] = v.val;
        }
        return out;
    }

    /**
     * find the bottom k elements of k, that is, the k smallest values in a.
     The runtime complexity is O(n * log(k)) where n is a.length
     * and the space cmplexity is O(k).
     * @param a
     * @param k
     * @return the k smallest elements of k in ascending order.
     */
    public static int[] bottomK(int[] a, int k) {
        TreeSet<Node<Integer>> heap = new TreeSet<>();

        for (int v : a) {
            heap.add(new Node<>(v));
            if (heap.size() > k) {
                heap.pollLast();
            }
        }
        int[] out = new int[Math.min(k, heap.size())];
        int iOut = 0;
        for (Node<Integer> v : heap) {
            out[iOut++] = v.val;
        }
        return out;
    }

    /**
     * O(n)
     * @param a
     * @param k
     * @return
     */
    //public static int kthSmallest(double[] a, int k) {

    //}
}