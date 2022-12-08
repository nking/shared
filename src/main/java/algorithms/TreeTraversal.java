package algorithms;

import algorithms.heapsAndPQs.HeapNode;
import java.util.ArrayDeque;
import java.util.Queue;
import java.util.Set;
import java.util.Stack;

/**
 *
 * @author nichole
 */
public class TreeTraversal {

    /**
     * root, left subtree, right subtree.
     <pre>
       e.g.
               0
            1         2
           3  4      5  6
         7  
        
      visits: 0, 1, 3, 7, 4, 2, 5, 6
     </pre>
     @param root root of tree
     */
    public void preorderRecursive(BinaryTreeNode root) {
        if (root != null) {
            System.out.println(root.getData());
            preorderRecursive(root.getLeft());
            preorderRecursive(root.getRight());
        }
    }

    /**
     * left subtree, root, right subtree
     <pre>
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                        
     visits: 7, 3, 1, 4, 10, 0, 5, 8, 9, 2, 6, 12, 11, 13
     </pre>
     @param root tree root
     */
    public void inorderRecursive(BinaryTreeNode root) {
        if (root != null) {
            inorderRecursive(root.getLeft());
            System.out.printf("%d, ", root.getData());
            inorderRecursive(root.getRight());
        }
    }

    /**
     * left subtree, right subtree, root
     <pre>
       e.g.
                  0
            1         2
           3  4      5  6
         7  
        
     visits: 7, 3, 4, 1, 5, 6, 2, 0
     </pre>
     @param root tree root
     */
    public void postorderRecursive(BinaryTreeNode root) {
        if (root != null) {
            postorderRecursive(root.getLeft());
            postorderRecursive(root.getRight());
            System.out.println(root.getData());
        }
    }

    /**
     * root, left subtree, right subtree
     <pre>
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                        
     visits: 0, 1, 3, 7, 4, 10, 2, 5, 8, 9, 6, 11, 12, 13
     </pre>
     @param node
     */
    public void preorderIterative(BinaryTreeNode node) {
        Stack<BinaryTreeNode> stack = new Stack<>();
        while (!stack.isEmpty() || node != null) {
            if (node != null) {
                System.out.printf("%d, ", node.getData());
                stack.push(node);
                node = node.getLeft();
            } else {
                node = stack.pop();
                node = node.getRight();
            }
        }
        System.out.println();
    }
    
    /**
     * left subtree, root, right subtree
     <pre>
     e.g.
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                        
     recursive visits: 7, 3, 1, 4, 10, 0, 5, 8, 9, 2, 6, 12, 11, 13
     iterative visits: 7, 3, 1, 4, 10, 0, 5, 8, 9, 2, 6, 12, 11, 13
     </pre>
     @param node
     */
    public void inorderIterative(BinaryTreeNode node) {
        Stack<BinaryTreeNode> stack = new Stack<>();
        int c = 0;
        while (!stack.isEmpty() || node != null) {
            c++;
            if (node != null) {
                stack.push(node);
                node = node.getLeft();
            } else {
                node = stack.pop();
                System.out.printf("%d, ", node.getData());
                node = node.getRight();
            }
        }
        System.out.printf("  nCalls=%d", c);
        System.out.println();
    }

    static class Snapshot {
        final BinaryTreeNode node;
        int stage;
        public Snapshot(BinaryTreeNode node, int stage) {
            this.node = node;
            this.stage = stage;
        }
    }
    /* left subtree, root, right subtree
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
    7, 3, 1, 4, 10, 0, 5, 8, 9, 2, 6, 12, 11, 13
         */

    /**
     *
     @param node
     */

    public void inorderIterative2(BinaryTreeNode node) {        
        // recursion to iteration: 
        //    https://www.codeproject.com/Articles/418776/How-to-replace-recursive-functions-using-stack-and
        Stack<Snapshot> s = new Stack<Snapshot>();
        s.push(new Snapshot(node, 0));
        int c = 0;
        Snapshot currentSnapshot;
        while (!s.isEmpty()) {
            c++;
            currentSnapshot = s.pop();
            // 2 recursive function calls, so 3 stages
            if (currentSnapshot.node != null) {
                switch(currentSnapshot.stage) {
                    case 0: {
                        currentSnapshot.stage++;
                        s.push(currentSnapshot);
                        s.push(new Snapshot(currentSnapshot.node.getLeft(), 0));                        
                        break;
                    }
                    case 1: {
                       System.out.printf("%d, ", currentSnapshot.node.getData());                        
                        currentSnapshot.stage++;
                        s.push(currentSnapshot);
                        break; 
                    }
                    case 2:{
                        s.push(new Snapshot(currentSnapshot.node.getRight(), 0));
                        break;
                    }
                }
            }
        }
        System.out.printf("  nCalls=%d", c);
        System.out.println();
    }
    
    /**
     * left subtree, right subtree, root
     <pre>
       e.g.
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                       
      7, 3, 10, 4, 1, 9, 8, 5, 12, 13, 11, 6, 2, 0
     </pre>
     this is actually bottom-up post-order iterative while recursive
     is top-down.
     TODO:  read this and convert here:
         https://www.cs.odu.edu/~zeil/cs361/latest/Public/recursionConversion/index.html
     @param node
     */
    public void postorderIterative(BinaryTreeNode node) {
        Stack<BinaryTreeNode> stack = new Stack<>();
        Stack<BinaryTreeNode> stack2 = new Stack<>();
        stack.push(node);
        while (!stack.isEmpty()) {
            node = stack.pop();
            stack2.push(node);
            if (node.getLeft() != null) {
                stack.push(node.getLeft());
            }
            if (node.getRight() != null) {
                stack.push(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.getData());
        }
        System.out.println();
    }

    /**
     * a.k.a. breadth first traversal
     <pre>
       e.g.
               0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                       
     0, 1, 2, 3, 4, 5, 6, 7, 10, 8, 11, 9, 12, 13           
     </pre>
     @param node tree root
     */
    public void levelOrderIterative(BinaryTreeNode node) {
        Queue<BinaryTreeNode> queue = new ArrayDeque<>();
        while (node != null) {
            System.out.printf("%d, ", node.getData());
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
            node = queue.poll(); // returns null if empty
        }
        System.out.println();
    }
    
    /**
     *
     @param node
     */
    public static void printLevelOrder(HeapNode node) {
        
        Queue<HeapNode> queue = new ArrayDeque<>();
        Queue<Integer> nodeLevel = new ArrayDeque<>();
        Queue<Long> nodeParent = new ArrayDeque<>();
        Queue<Character> leftOrRight = new ArrayDeque<>();
        
        int level = 0;
        Character lOrR = '-';
        Long parent = -1L;
        
        while (node != null) {
            
            System.out.printf("L=%d, [%d, %s], %s of key %d\n", 
                level, node.getKey(), node.getData() != null ?
                ((Integer)node.getData()).toString() : "-",
                lOrR, parent);
            
            ++level;
            
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
                nodeLevel.add(level);
                nodeParent.add(node.getKey());
                leftOrRight.add('L');
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
                nodeLevel.add(level);
                nodeParent.add(node.getKey());
                leftOrRight.add('R');
            }
            //TODO: rearrange to check empty just once instead of at while loop too
            if (queue.isEmpty() || nodeLevel.isEmpty()) {
                break;
            }
            node = queue.poll(); // returns null if empty
            level = nodeLevel.poll();
            lOrR = leftOrRight.poll();
            parent = nodeParent.poll();
        }
    }


    /**
     * get the reverse level-order traversal of tree node.
     * implemented as post-order traversal but using a queue for the first
     * stack:
     * adapted from https://www.geeksforgeeks.org/inorder-tree-traversal-without-recursion/?ref=gcse 
     <pre>
       e.g.
               0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                        * 
      13, 12, 9, 11, 8, 10, 7, 6, 5, 4, 3, 2, 1, 0
     </pre>
     @param node
     */
    public void reverseLevelOrderIterative(BinaryTreeNode node) {
        Queue<BinaryTreeNode> queue = new ArrayDeque<>();//FIFO
        Stack<BinaryTreeNode> stack2 = new Stack<>();
        queue.add(node);
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.getData());
        }
        System.out.println();
    }
    
    /**
     * given a tree represented by node, return a doubly-linked list of nodes
     * visited in a reverse level-order traversal.
     * 
     @param node the tree to be traversed from bottom up to this node using
     * reverse level-order traversal.
     * NOTE that any next and prev links in the tree are overwritten by the
     * DoublyLinkedList, so if those need to be preserved, give this method
     * node.copyTree() instead.
     @return a double-linked list of nodes in reverse level order traversal.
     */
    public DoublyLinkedList<NAryTreeNode> getReverseLevelOrderIterative2(NAryTreeNode node) {
        
        DoublyLinkedList<NAryTreeNode> out = new DoublyLinkedList<NAryTreeNode>();
        
        Queue<NAryTreeNode> queue = new ArrayDeque<>();//FIFO
        Stack<NAryTreeNode> stack2 = new Stack<>();
        queue.add(node);
        
        Set<NAryTreeNode> children;
        
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            children = node.getChildren();
            if (children == null) {
                continue;
            }
            queue.addAll(children);
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            out.add(node);
        }
        return out;
    }
    
    /**
     * a.k.a. breadth first traversal
     <pre>
       e.g.
               0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
                       
     0, 1, 2, 3, 4, 5, 6, 7, 10, 8, 11, 9, 12, 13           
     </pre>
     @param node n-ary tree root
     @return 
     */
    public DoublyLinkedList<NAryTreeNode> getLevelOrderIterative(NAryTreeNode node) {
        Queue<NAryTreeNode> queue = new ArrayDeque<>();
        DoublyLinkedList<NAryTreeNode> out = new DoublyLinkedList<NAryTreeNode>();
        Set<NAryTreeNode> children;
        while (node != null) {
            out.add(node);
            children = node.getChildren();
            if (children == null) {
                continue;
            }
            queue.addAll(children);

            node = queue.poll(); // returns null if empty
        }
        return out;
    }
    
    /**
     * 
     @param node tree root
     @return a double-linked list of nodes in reverse level order traversal.
     */
    public DoublyLinkedList<BinaryTreeNode> getReverseLevelOrderIterative(BinaryTreeNode node) {
        DoublyLinkedList<BinaryTreeNode> out = new DoublyLinkedList<BinaryTreeNode>();
        Queue<BinaryTreeNode> queue = new ArrayDeque<>();//FIFO
        Stack<BinaryTreeNode> stack2 = new Stack<>();
        queue.add(node);
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            if (node.getLeft() != null) {
                queue.add(node.getLeft());
            }
            if (node.getRight() != null) {
                queue.add(node.getRight());
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            out.add(node);
        }
        return out;
    }

}
