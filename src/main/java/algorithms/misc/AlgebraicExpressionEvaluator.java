package algorithms.misc;

import java.util.*;

/**
 *
 * class to evaluate an algebraic expression using a reduce-shift parser pattern to build parse tree
 * and then uses reverse level order traversal to evaluate the parse tree.

 @author Nichole
 */
public class AlgebraicExpressionEvaluator {

    protected static Map<Character, Integer> opAlgPriorityMap = new HashMap<>();
    static {
        opAlgPriorityMap.put('^', 3);
        opAlgPriorityMap.put('*', 2);
        opAlgPriorityMap.put('/', 2);
        opAlgPriorityMap.put('+', 1);
        opAlgPriorityMap.put('-', 1);
        opAlgPriorityMap.put('(', 0);
        opAlgPriorityMap.put(')', 0);
    }

    /**
     * evaluates the algebraic expression.
     *
     * Recognized operators are +, -, *, /,  and ^ as power operator.
     * The method does not evaluate binary operators (e.g. |, &&, or ^ as XOR).
     * The method handles parenthesization in context of use in ordering and prioritizing operations,
     * but does not handle square brackets or other separators.
     * The meothd can handle spaces.
     *
     * The runtime exception IllegalArgumentException is thrown if expression isn't a valid format,
     * such as having imbalanced parentheses or other syntax errors.
     *
     * NOTE that algebraic order of operations give highest priority to power operator,
     * then to multiplicaiton and division equally, then to addition and subtraction equally.
     *
     * Adding parentheses can be used to change priority and order of operations.
     *
     <pre>
     The expression tree method was built following the article
     https://www.geeksforgeeks.org/program-to-convert-infix-notation-to-expression-tree/
     authored by Julkar9, and following the java code contributed by aashish1995.

     The building of the expression tree combines a conversion of the input infix string to a postfix format
     using a stack.  The postfix format is itself a tree built with a stack.

     Here I adapted the abstract syntax tree builder to lex and parse numbers and I added a reverse level
     order evaluator.
     </pre>

     * @param string algebraic expression that can include parentheses.  the string must be in infix format.
     *               infix format places the operator in between the two operands (numbers only for this method).
     * In contrast, prefix format places the operator before each pair of operands and
     * postfix places the operator after the operands.
     */
    public double evaluateAlgebraicExpression(String expression) {

        char[] chars = expression.toCharArray();

        // trim off any trailing semi-colons from end of expression.
        int j = chars.length - 1;
        while (j > -1 && (Character.isSpaceChar(chars[j]) || (chars[j] == ';'))) --j;

        if (j == 0) {
            if (Character.isSpaceChar(chars[0])) return 0;
            else return Double.parseDouble(String.copyValueOf(chars,0, 1));
        };

        int i = 0;
        while (i < j && Character.isSpaceChar(chars[i])) ++i;

        OpPair root = buildExprTree(chars, i, j);

        return reverseLevelOrderEval(root, chars);
    }

    /**
     * build an algebraic expression tree from the input string that has infix format.
     *
     <pre>
     This method is adapted from code available at
     https://www.geeksforgeeks.org/program-to-convert-infix-notation-to-expression-tree/
     authored by Julkar9, and aashish1995
     </pre>
     * @param chars
     * @param i
     * @param j
     * @return
     */
    protected OpPair buildExprTree(char[] chars, int i, int j) {

        //numberStack is the equiv of the postfix format string
        Stack<OpPair> numberStack = new Stack<>();
        Stack<CRange> charStack = new Stack<>();

        OpPair t;
        int[] i01 = new int[2];
        /*
        while (i <= j) {
            if c == '('
                cStack.push(c)
            else if c == ')' {
                while (!cStack.isEmpty() && cStack.peek != '(') {
                    t = new OpPair(cStack.pop())
                    t.right = nStack.pop();
                    t.left = nStack.pop();
     =>              nStack.push(t);
                }
                cStack.pop();//this is the opening parenthesis
            } else if c is an operator {
                while (!cStack.isEmpty() && cStack.peek != '(' && obeys precedence) {
                    t = new OpPair(cStack.pop());
                    t.right = nStack.pop();
                    t.left = nStack.pop();
      -->            nStack.push(t);
                }
                cStack.push(new CRange(i, i));
            } else {
                c is a number
      =>          nStack.push(new OPair(i,i))
            }

        }
         */

        boolean prevWasNumber = false;
        while (i <= j) {
            if (chars[i] == '(') {

                charStack.add(new CRange(i, i));
                prevWasNumber = false;

            } else if (chars[i] == ')') {
                prevWasNumber = true;//result is number

                // a Reduce stage of shift-reduce parsing

                while (!charStack.isEmpty() && isNotLeftParenthesis(charStack.peek(), chars)) {

                    t = new OpPair(charStack.pop());
                    t.right = numberStack.pop();
                    t.left = numberStack.pop();
                    numberStack.push(t);
                }
                // pop the opening parenthesis:
                charStack.pop();

            } else if (prevWasNumber && opAlgPriorityMap.containsKey(chars[i])) {
                prevWasNumber = false;

                // a Reduce stage of shift-reduce parsing

                while (!charStack.isEmpty() && isNotLeftParenthesis(charStack.peek(), chars)
                    && ((chars[i] != '^' && priorityisGEQ(charStack.peek(), i, chars))
                    || (chars[i] == '^' && priorityisGT(charStack.peek(), i, chars)))) {

                    t = new OpPair(charStack.pop());
                    t.right = numberStack.pop();
                    t.left = numberStack.pop();

                    numberStack.push(t);
                }
                charStack.push(new CRange(i, i));
            } else {

                // a Shift stage of shift-reduce parsing

                // we have positive or negative numbers that
                // may begin with a . if the user dropped the 0 as an integer in a decimal number.
                prevWasNumber = true;
                parseForNumber(i, j, chars, i01);
                i = i01[1];

                t = new OpPair(new CRange(i01[0], i01[1]));
                numberStack.push(t);
            }

            ++i;
            while (i < j && Character.isSpaceChar(chars[i])) ++i;
        }

        // the Accept stage of shift-reduce parsing

        assert(charStack.isEmpty());
        assert(numberStack.size() == 1);
        return numberStack.pop();
    }

    protected void parseForNumber(int i, int j, char[] chars, int[] outputIdxs) {
        // parse number starting at i, and do not parse beyond j
        // and put starting and ending indexes in outputIdxs

        while (i < j && Character.isSpaceChar(chars[i])) ++i;
        outputIdxs[0] = i;
        if (chars[i] == '-') ++i;
        if (chars[i] == '+') ++i;
        if (chars[i] == '.') ++i;
        while (((i+1) < j) && (Character.isDigit(chars[i+1]) || chars[i+1] == '.' )) ++i;

        outputIdxs[1] = i;
    }

    protected double reverseLevelOrderEval(OpPair node, char[] chars) {
        Queue<OpPair> q = new ArrayDeque<>();
        Stack<OpPair> s = new Stack<>();
        q.offer(node);

        while (!q.isEmpty()) {
            node = q.poll();
            s.push(node);
            if (node.left != null) q.offer(node.left);
            if (node.right != null) q.offer(node.right);
        }

        // key = node hashcode, value=result
        Map<Integer, Double> resultsMap = new HashMap<>();

        double result=0;
        char op;
        while (!s.isEmpty()) {
            node = s.pop();

            if (node.left == null && node.right==null) {
                // leaf nodes have literals stored in data
                resultsMap.put(node.hashCode(), parseNumber(chars, node.data));
            } else {
                assert(node.data.i0 == node.data.i1);
                op = chars[node.data.i0];
                double v1;

                if (resultsMap.containsKey(node.left.hashCode())) {
                    v1 = resultsMap.get(node.left.hashCode());
                } else {
                    v1 = parseNumber(chars, node.left.data);
                }
                double v2;
                if (resultsMap.containsKey(node.right.hashCode())) {
                    v2 = resultsMap.get(node.right.hashCode());
                } else {
                    v2 = parseNumber(chars, node.right.data);
                }

                switch(op) {
                    case '^': {
                        result = Math.pow(v1, v2);
                        resultsMap.put(node.hashCode(), result);
                        break;
                    }
                    case '*': {
                        result = v1 * v2;
                        resultsMap.put(node.hashCode(), result);
                        break;
                    }
                    case '/': {
                        result = v1 / v2;
                        resultsMap.put(node.hashCode(), result);
                        break;
                    }
                    case '+': {
                        result = v1 + v2;
                        resultsMap.put(node.hashCode(), result);
                        break;
                    }
                    case '-': {
                        result = v1 - v2;
                        resultsMap.put(node.hashCode(), result);
                        break;
                    }
                    default:
                        throw new IllegalArgumentException("ERROR: operator not recognized: " + node.data);
                }// end switch
            }// end else

        }//end while
        return result;
    }
    public String toString(char[] chars, CRange cRange) {
        return String.copyValueOf(chars, cRange.i0, cRange.i1 - cRange.i0 + 1);
    }

    protected double parseNumber(char[] chars, OpPair node) {
        return parseNumber(chars, node.data);
    }

    protected double parseNumber(char[] chars, CRange nodeData) {
        String s = String.copyValueOf(chars, nodeData.i0, nodeData.i1 - nodeData.i0 + 1);
        return Double.parseDouble(s);
    }

    protected boolean priorityisGT(CRange c0, int i1, char[] chars) {
        assert(c0.i0 == c0.i1);
        int p0 = opAlgPriorityMap.get(chars[c0.i0]);
        int p1 = opAlgPriorityMap.get(chars[i1]);
        return p0 > p1;
    }
    protected boolean priorityisGEQ(CRange c0, int i1, char[] chars) {
        assert(c0.i0 == c0.i1);
        int p0 = opAlgPriorityMap.get(chars[c0.i0]);
        int p1 = opAlgPriorityMap.get(chars[i1]);
        return p0 >= p1;
    }

    protected boolean isNotLeftParenthesis(CRange c, char[] chars) {
        return !(c.i0 == c.i1 && chars[c.i0] == '(');
    }

    protected static class OpPair {
        OpPair left = null;
        OpPair right = null;
        CRange data;
        public OpPair(CRange cr) {
            this.data = cr;
        }
    }

    protected static class CRange {
        final int i0; // beginning index of operator or number
        final int i1; // ending index of operator or number, inclusive
        public CRange(int idx0, int idx1) {
            this.i0 = idx0;
            this.i1 = idx1;
        }
    }

}
