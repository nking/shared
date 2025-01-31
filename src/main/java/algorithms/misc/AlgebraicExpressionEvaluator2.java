package algorithms.misc;

import java.util.*;

/**
 class to evaluate an algebraic expression using a shift-reduce parse pattern to build
 a parse tree. The tree nodes are evaluated in the reduce stages.

 @author Nichole
 */
public class AlgebraicExpressionEvaluator2 {

    /**
     * sentinel
     */
    protected static double sentinel = Double.POSITIVE_INFINITY;

    /**
     * map of operation priorities
     */
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
     * The method does not evaluate binary operators (e.g. |, &amp;&amp;, or ^ as XOR).
     * The method handles parenthesization in context of use in ordering and prioritizing operations,
     * but not handle square brackets or other separators.
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

     * @param expression string algebraic expression that can include parentheses.  the string must be in infix format.
     *               infix format places the operator in between the two operands (numbers only for this method).
     * In contrast, prefix format places the operator before each pair of operands and
     * postfix places the operator after the operands.
     * @return evaluation
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

        return parseAndEvaluate(chars, i, j);
    }

    /**
     * build an algebraic expression tree from the input string that has infix format and evaluate
     * the expression in the reduce stages.

     <pre>
     This method is adapted from code available at
     https://www.geeksforgeeks.org/program-to-convert-infix-notation-to-expression-tree/
     authored by Julkar9, and aashish1995
     but edited to add evaluation.

     </pre>
     * @param chars character array
     * @param i start index of chars for use in method
     * @param j stop index, inclusive of chars  for use in method
     * @return evaluation
     */
    protected double parseAndEvaluate(char[] chars, int i, int j) {

        //numberStack is the equiv of the postfix format string
        Stack<Double> numberStack = new Stack<>();
        Stack<CRange> charStack = new Stack<>();

        double result = 0.;

        int[] i01 = new int[2];

        boolean prevWasNumber = false;
        while (i <= j) {
            if (chars[i] == '(') {

                charStack.add(new CRange(i, i));
                prevWasNumber = false;

            } else if (chars[i] == ')') {
                prevWasNumber = true;//result is number

                // a Reduce stage of shift-reduce parsing

                while (!charStack.isEmpty() && isNotLeftParenthesis(charStack.peek(), chars)) {
                    result = evaluate(chars, charStack.pop(), numberStack.pop(), numberStack.pop());
                    numberStack.push(result);
                }
                // pop the opening parenthesis:
                charStack.pop();

            } else if (prevWasNumber && opAlgPriorityMap.containsKey(chars[i])) {
                prevWasNumber = false;

                // a Reduce stage of shift-reduce parsing

                while (!charStack.isEmpty() && isNotLeftParenthesis(charStack.peek(), chars)
                    && ((chars[i] != '^' && priorityisGEQ(charStack.peek(), i, chars))
                    || (chars[i] == '^' && priorityisGT(charStack.peek(), i, chars)))) {

                    result = evaluate(chars, charStack.pop(), numberStack.pop(), numberStack.pop());
                    numberStack.push(result);
                }
                charStack.push(new CRange(i, i));

            } else {

                // a Shift stage of shift-reduce parsing

                // we have positive or negative numbers that
                // may begin with a . if the user dropped the 0 as an integer.
                prevWasNumber = true;
                parseForNumber(i, j, chars, i01);
                i = i01[1];

                numberStack.push( parseNumber(chars, new CRange(i01[0], i01[1])) );
            }

            ++i;
            while (i < j && Character.isSpaceChar(chars[i])) ++i;
        }

        // the Accept stage of shift-reduce parsing

        assert(charStack.isEmpty());
        assert(numberStack.size() == 1);
        return result;
    }

    /**
     * evaluates vLeft operator vRight.
     * @param chars array of characters
     * @param op operator
     * @param vRight right value
     * @param vLeft left value
     * @return result of vLeft operator vRight
     */
    private double evaluate(char[] chars, CRange op, double vRight, double vLeft) {
        double result = 0;

        switch(chars[op.i0]) {
            case '^': {
                result = Math.pow(vLeft, vRight);
                break;
            }
            case '*': {
                result = vLeft * vRight;
                break;
            }
            case '/': {
                result = vLeft / vRight;
                break;
            }
            case '+': {
                result = vLeft + vRight;
                break;
            }
            case '-': {
                result = vLeft - vRight;
                break;
            }
            default:
                throw new IllegalArgumentException("ERROR: operator not recognized: " + chars[op.i0]);
        }

        return result;
    }

    /**
     *
     * @param chars character array
     * @param i start index of chars for use in method
     * @param j stop index, inclusive of chars  for use in method
     * @param outputIdxs the starting and ending indexes, inclusive
     *                   of the number just parsed from chars
     */
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

    /**
     * string of chars bounded by cRange
     * @param chars character array
     * @param cRange node olding index range
     * @return string of cars bounded by cRange
     */
    public String toString(char[] chars, CRange cRange) {
        return String.copyValueOf(chars, cRange.i0, cRange.i1 - cRange.i0 + 1);
    }

    /**
     * number from chars[nodeData.i0] to chars[nodeData.i1] inclusive
     * @param chars character array
     * @param nodeData daa defining character range of interest
     * @return number from chars[nodeData.i0] to chars[nodeData.i1] inclusive
     */
    protected double parseNumber(char[] chars, CRange nodeData) {
        String s = String.copyValueOf(chars, nodeData.i0, nodeData.i1 - nodeData.i0 + 1);
        return Double.parseDouble(s);
    }

    /**
     * determine if the operation priority for chars[c0.i0] is greater than
     * priority for chars[i1]
     * @param c0 1st index in comparison is co.i0
     * @param i1 2nd index in comparison
     * @param chars character array
     * @return true if operation priority for chars[c0.i0] is greater than
     *   priority for chars[i1]
     */
    protected boolean priorityisGT(CRange c0, int i1, char[] chars) {
        assert(c0.i0 == c0.i1);
        int p0 = opAlgPriorityMap.get(chars[c0.i0]);
        int p1 = opAlgPriorityMap.get(chars[i1]);
        return p0 > p1;
    }

    /**
     * determine if the operation priority for chars[c0.i0] is greater than
     * or equal to priority for chars[i1]
     * @param c0 1st index in comparison is co.i0
     * @param i1 2nd index in comparison
     * @param chars character array
     * @return true if operation priority for chars[c0.i0] is greater than
     *   or equal to priority for chars[i1]
     */
    protected boolean priorityisGEQ(CRange c0, int i1, char[] chars) {
        assert(c0.i0 == c0.i1);
        int p0 = opAlgPriorityMap.get(chars[c0.i0]);
        int p1 = opAlgPriorityMap.get(chars[i1]);
        return p0 >= p1;
    }

    /**
     * check c.i0 == '('
     * @param c range
     * @param chars character array
     * @return true if c.i0 == '('
     */
    protected boolean isNotLeftParenthesis(CRange c, char[] chars) {
        return !(c.i0 == c.i1 && chars[c.i0] == '(');
    }

    /**
     * holds an index range
     */
    protected static class CRange {
        /**
         * beginning index of operator or number
         */
        final int i0;
        /**
         * ending index, inclusive, or operator or number
         */
        final int i1;

        /**
         * range node
         * @param idx0 start index
         * @param idx1 stop index, inclusive
         */
        public CRange(int idx0, int idx1) {
            this.i0 = idx0;
            this.i1 = idx1;
        }
    }

}
