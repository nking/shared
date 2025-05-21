package algorithms;

/**
 * data structure which is relatively easy to implement and has O(log(n))
 * r.t.c on average for all operations, though worse case r.t.c. is O(n).
 *
 * the implementation follows the paper by William Pugh, 1989
 <pre>
 reference:
 "A SkipList Cookbook", by Willianm Pugh, 1989
 UMIACS-TR-89-72.1
 CS-TR-2286.1
 </pre>
 */
public class SkipList {

    protected final double P;
    protected final int MAX_LEVEL;
    protected int level;
    protected Node header;

    protected static class Node {
        Object data;
        final int key;
        final Node[] forward;
        public Node(int key, Object data, int nLevels) {
            this.key = key;
            this.data = data;
            forward = new Node[nLevels];
        }
    }

    /**
     * ideally, maxLevel should be log_{1/p}(n) where n is the number of expected elements
     * @param p
     * @param maxLevel
     */
    public SkipList(double p, int maxLevel) {
        this.P = p;
        this.MAX_LEVEL = maxLevel;
        this.header = new Node(-1, null, MAX_LEVEL+1);
    }

    /**
     * given expected number of elements in list and a probability threshold
     * the method internally calculates the ideal max level os
     * log_{1/p}(n) where n is the number of expected elements
     *
     * @param expectedNumberOfElements
     * @param p
     */
    public SkipList(int expectedNumberOfElements, double p) {
        this.P = p;
        this.MAX_LEVEL = (int)Math.ceil(Math.log(
                expectedNumberOfElements)/Math.log(1./P));
        this.header = new Node(-1, null, MAX_LEVEL+1);
    }

    /**
     * given expected number of elements in list, the method
     * internally sets a probability threshold to 0.5 and
     * calculates the ideal max level os
     * log_{1/p}(n) where n is the number of expected elements
     *
     * @param expectedNumberOfElements
     * @param p
     */
    public SkipList(int expectedNumberOfElements) {
        this.P = 0.5;
        this.MAX_LEVEL = (int)Math.ceil(Math.log(
                expectedNumberOfElements)/Math.log(1./P));
        this.header = new Node(-1, null, MAX_LEVEL+1);
    }

    /**
     *
     * @param key
     * @param update can be null.  if not null, should be length MAX_LEVEL
     * @return
     */
    protected Node findNode(int key, Node[] update) {
        Node x = header;
        for (int i = level; i >= 0; --i) {
            while (x.forward[i] != null && x.forward[i].key < key) {
                x = x.forward[i];
            }
            if (update != null) {
                update[i] = x;
            }
        }
        x = x.forward[0];
        return x;
    }

    public Object search(int key) {
        Node x = findNode(key, null);

        if (x != null && x.key == key) {
            return x.data;
        }
        return null;
    }

    /**
     * insert key into skiplist.  NOTE that if key exists, the new data overwrites
     * the old data and the old data is returned.
     * @param key
     * @param data
     * @return after inserting the key into the list,
     * returns former data for given key if key was already in list,
     * else returns null if key was not in list.
     */
    public Object insert(int key, Object data) {
        Node[] update = new Node[MAX_LEVEL];
        Node x = findNode(key, update);

        if (x != null && x.key == key) {
            Object old = x.data;
            x.data = data;
            return old;
        }

        int newLevel = randomLevel();
        if (newLevel > level) {
            for (int i = level + 1; i <= newLevel; ++i) {
                update[i] = header;
            }
            level = newLevel;
        }
        x = new Node(key, data, newLevel);
        for (int i = 0; i <= newLevel; ++i) {
            x.forward[i] = update[i].forward[i];
            update[i].forward[i] = x;
        }
        return null;
    }

    private int randomLevel() {
        int lvl = 0;
        while (Math.random() < P && lvl < MAX_LEVEL) {
            ++lvl;
        }
        return lvl;
    }

    /**
     * insert key into skiplist.  NOTE that if key exists, the new data overwrites
     * the old data and the old data is returned.
     * @param key
     * @param data
     * @return after inserting the key into the list,
     * returns former data for given key if key was already in list,
     * else returns null if key was not in list.
     */
    public Object delete(int key) {
        Node[] update = new Node[MAX_LEVEL];
        Node x = findNode(key, update);
        Object data = null;
        if (x != null && x.key == key) {
            data = x.data;
            for (int i = 0; i <= level; ++i) {
                if (!update[i].forward[i].equals(x)) {
                    break;
                }
                update[i].forward[i] = x.forward[i];
            }
            while (level > 0 && header.forward[level] == null) {
                --level;
            }
        }

        int newLevel = randomLevel();
        if (newLevel > level) {
            for (int i = level + 1; i <= newLevel; ++i) {
                update[i] = header;
            }
            level = newLevel;
        }
        x = new Node(key, data, newLevel);
        for (int i = 0; i < newLevel; ++i) {
            x.forward[i] = update[i].forward[i];
            update[i].forward[i] = x;
        }
        return data;
    }
}
