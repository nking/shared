package algorithms;

import java.util.ArrayDeque;
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
            preorderRecursive(root.left);
            System.out.println(root.data);
            preorderRecursive(root.right);
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
            preorderRecursive(root.left);
            preorderRecursive(root.right);
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
            1         2
           3  4      5  6
         7  
         
     recursive visits: 1, 3, 7, 4, 0, 2, 5, 6
     iterative visits: 7, 3, 1, 4, 0, 5, 2, 6
     </pre>
     */
    public void inorderIterative(Node node) {
        /*Stack<Node> stack = new Stack<>();
        while (!stack.isEmpty() || node != null) {
            if (node != null) {
                stack.push(node);
                node = node.left;
            } else {
                node = stack.pop();
                System.out.printf("%d, ", node.data);
                node = node.right;
            }
        }
        System.out.println();
        */
//TODO: visiting right subtree is not correct yet     
        Stack<Node> stack = new Stack<>();
        Stack<Node> stackL = new Stack<>();
        Node p;
        while (!stack.isEmpty() || node != null) {
            if (node != null) {
                stack.add(node);
                node = node.left;
            } else {
                node = stack.pop(); 
                while (node.right == null && 
                    node.parent != null && !stack.isEmpty() && (stack.peek().data == node.parent.data)) {
                    stackL.add(node);
                    node = stack.pop();
                    if (node.right != null) {
                        break;
                    }
                }
                System.out.printf("%d, ", node.data);
                node = node.right;
                while (!stackL.isEmpty()) {
                    Node n2 = stackL.pop();
                    System.out.printf("%d, ", n2.data);
                }
            }
        }
        System.out.println();
        /*
                0
            1           2
           3  4      5     6
         7            8
                        9
        
        */
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
    
    public void test1() {

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

}
