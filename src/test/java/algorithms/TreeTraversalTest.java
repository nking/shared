package algorithms;

import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayDeque;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Stack;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class TreeTraversalTest extends TestCase {

    static class Node {
        int data;
        Node left, right,parent;
        Node(int data) {
            this.data = data;
        }
        public String toString() {
            return Integer.toString(data);
        }
    }

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
     * @param root
     */
    public void preorderRecursive(Node root) {
        if (root != null) {
            System.out.println(root.data);
            preorderRecursive(root.left);
            preorderRecursive(root.right);
        }
    }

    /**
     * left subtree, root, right subtree
     <pre>
       e.g.
                 0
            1         2
           3  4      5  6
         7  
        
      visits: 1, 3, 7, 4, 0, 2, 5, 6
     </pre>
     * @param root
     */
    public void inorderRecursive(Node root) {
        if (root != null) {
            inorderRecursive(root.left);
            System.out.printf("%d, ", root.data);
            inorderRecursive(root.right);
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
     * @param root
     */
    public void postorderRecursive(Node root) {
        if (root != null) {
            postorderRecursive(root.left);
            postorderRecursive(root.right);
            System.out.println(root.data);
        }
    }

    /**
     * root, left subtree, right subtree
     <pre>
       e.g.
               0
            1         2
           3  4      5  6
         7  
        
     visits: 0, 1, 3, 7, 4, 2, 5, 6
     </pre>
     */
    public void preorderIterative(Node node) {
        Stack<Node> stack = new Stack<>();
        while (!stack.isEmpty() || node != null) {
            if (node != null) {
                System.out.printf("%d, ", node.data);
                stack.push(node);
                node = node.left;
            } else {
                node = stack.pop();
                node = node.right;
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
     */
    public void inorderIterative(Node node) {
        Stack<Node> stack = new Stack<>();
        int c = 0;
        while (!stack.isEmpty() || node != null) {
            c++;
            if (node != null) {
                stack.push(node);
                node = node.left;
            } else {
                node = stack.pop();
                System.out.printf("%d, ", node.data);
                node = node.right;
            }
        }
        System.out.printf("  nCalls=%d", c);
        System.out.println();
    }

    static class Snapshot {
        final Node node;
        int stage;
        public Snapshot(Node node, int stage) {
            this.node = node;
            this.stage = stage;
        }
    }
    /*
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
         */
    public void inorderIterative2(Node node) {        
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
                        s.push(new Snapshot(currentSnapshot.node.left, 0));                        
                        break;
                    }
                    case 1: {
                       System.out.printf("%d, ", currentSnapshot.node.data);                        
                        currentSnapshot.stage++;
                        s.push(currentSnapshot);
                        break; 
                    }
                    case 2:{
                        s.push(new Snapshot(currentSnapshot.node.right, 0));
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
            1         2
           3  4      5  6
         7  
     </pre>
     this is actually bottom-up post-order iterative while recursive
     is top-down.
     TODO:  read this and convert here:
         https://www.cs.odu.edu/~zeil/cs361/latest/Public/recursionConversion/index.html
     */
    public void postorderIterative(Node node) {
        Stack<Node> stack = new Stack<>();
        Stack<Node> stack2 = new Stack<>();
        stack.push(node);
        while (!stack.isEmpty()) {
            node = stack.pop();
            stack2.push(node);
            if (node.left != null) {
                stack.push(node.left);
            }
            if (node.right != null) {
                stack.push(node.right);
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.data);
        }
        System.out.println();
    }

    /**
     * a.k.a. breadth first traversal
     <pre>
       e.g.
               0
            1         2
           3  4      5  6
         7  
     </pre>
     * @param node
     */
    public void levelOrderIterative(Node node) {
        Queue<Node> queue = new ArrayDeque<>();
        while (node != null) {
            System.out.printf("%d, ", node.data);
            if (node.left != null) {
                queue.add(node.left);
            }
            if (node.right != null) {
                queue.add(node.right);
            }
            node = queue.poll(); // returns null if empty
        }
        System.out.println();
    }

    /**
     * implemented as post-order traversal but using a queue for the first
     * stack.
     * adapted from https://www.geeksforgeeks.org/inorder-tree-traversal-without-recursion/?ref=gcse 
     <pre>
       e.g.
               0
            1         2
           3  4      5  6
         7  
     </pre>
     */
    public void reverseLevelOrderIterative(Node node) {
        Queue<Node> queue = new ArrayDeque<>();//FIFO
        Stack<Node> stack2 = new Stack<>();
        queue.add(node);
        while (!queue.isEmpty()) {
            node = queue.remove();//retrieves first in list
            stack2.push(node);
            if (node.left != null) {
                queue.add(node.left);
            }
            if (node.right != null) {
                queue.add(node.right);
            }
        }
        while (!stack2.isEmpty()) {
            node = stack2.pop();
            System.out.printf("%d, ", node.data);
        }
        System.out.println();
    }

    public void est0() {

        /*
                0
            1         2
           3  4      5  6
         7                
         */
        Node[] nodes = new Node[8];
        int i;
        for (i = 0; i < nodes.length; ++i) {
            nodes[i] = new Node(i);
        }
        nodes[0].left = nodes[1];
        nodes[0].right = nodes[2];
        nodes[1].left = nodes[3];
        nodes[1].right = nodes[4];
        nodes[2].left = nodes[5];
        nodes[2].right = nodes[6];
        nodes[3].left = nodes[7];
        
        nodes[1].parent = nodes[0];
        nodes[2].parent = nodes[0];
        nodes[3].parent = nodes[1];
        nodes[4].parent = nodes[1];
        nodes[7].parent = nodes[3];
        nodes[5].parent = nodes[2];
        nodes[6].parent = nodes[2];

        System.out.println("pre-order:");
        preorderRecursive(nodes[0]);
        preorderIterative(nodes[0]);
        
        System.out.println("in-order:");
        inorderRecursive(nodes[0]);
        inorderIterative(nodes[0]);

        System.out.println("post-order:");
        postorderRecursive(nodes[0]);
        postorderIterative(nodes[0]);
        
        System.out.println("level-order:");
        levelOrderIterative(nodes[0]);

        System.out.println("reverse-level-order:");
        reverseLevelOrderIterative(nodes[0]);
    }
    
    public void est1() {

        /*
                0
            1           2
           3  4      5     6
         7            8
                        9
         */
        Node[] nodes = new Node[10];
        int i;
        for (i = 0; i < nodes.length; ++i) {
            nodes[i] = new Node(i);
        }
        nodes[0].left = nodes[1];
        nodes[0].right = nodes[2];
        nodes[1].left = nodes[3];
        nodes[1].right = nodes[4];
        nodes[2].left = nodes[5];
        nodes[2].right = nodes[6];
        nodes[3].left = nodes[7];
        nodes[5].right = nodes[8];
        nodes[8].right = nodes[9];
        
        nodes[1].parent = nodes[0];
        nodes[2].parent = nodes[0];
        nodes[3].parent = nodes[1];
        nodes[4].parent = nodes[1];
        nodes[7].parent = nodes[3];
        nodes[5].parent = nodes[2];
        nodes[6].parent = nodes[2];
        nodes[8].parent = nodes[5];
        nodes[9].parent = nodes[8];

        System.out.println("pre-order:");
        preorderRecursive(nodes[0]);
        preorderIterative(nodes[0]);
        
        System.out.println("in-order:");
        inorderRecursive(nodes[0]);
        inorderIterative(nodes[0]);

        System.out.println("post-order:");
        postorderRecursive(nodes[0]);
        postorderIterative(nodes[0]);
        
        System.out.println("level-order:");
        levelOrderIterative(nodes[0]);

        System.out.println("reverse-level-order:");
        reverseLevelOrderIterative(nodes[0]);
    }
    
    public void test2() {

        /*
                0
            1           2
           3  4      5     6
         7     10     8      11
                        9   12 13
         */
        Node[] nodes = new Node[14];
        int i;
        for (i = 0; i < nodes.length; ++i) {
            nodes[i] = new Node(i);
        }
        nodes[0].left = nodes[1];
        nodes[0].right = nodes[2];
        nodes[1].left = nodes[3];
        nodes[1].right = nodes[4];
        nodes[2].left = nodes[5];
        nodes[2].right = nodes[6];
        nodes[3].left = nodes[7];
        nodes[4].right = nodes[10];
        nodes[5].right = nodes[8];
        nodes[6].right = nodes[11];
        nodes[8].right = nodes[9];
        nodes[11].left = nodes[12];
        nodes[11].right = nodes[13];
        
        nodes[1].parent = nodes[0];
        nodes[2].parent = nodes[0];
        nodes[3].parent = nodes[1];
        nodes[4].parent = nodes[1];
        nodes[7].parent = nodes[3];
        nodes[5].parent = nodes[2];
        nodes[6].parent = nodes[2];
        nodes[8].parent = nodes[5];
        nodes[9].parent = nodes[8];
        nodes[10].parent = nodes[4];
        nodes[11].parent = nodes[6];
        nodes[12].parent = nodes[11];
        nodes[13].parent = nodes[11];

        System.out.println("pre-order:");
        preorderRecursive(nodes[0]);
        preorderIterative(nodes[0]);
        
        System.out.println("in-order:");
        inorderRecursive(nodes[0]); System.out.println();
        inorderIterative(nodes[0]);
        inorderIterative2(nodes[0]);

        System.out.println("post-order:");
        postorderRecursive(nodes[0]);
        postorderIterative(nodes[0]);
        
        System.out.println("level-order:");
        levelOrderIterative(nodes[0]);

        System.out.println("reverse-level-order:");
        reverseLevelOrderIterative(nodes[0]);

    }

}
