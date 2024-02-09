package algorithms.misc;

import junit.framework.TestCase;

public class AlgebraicExpressionEvaluatorTest extends TestCase {

    public void test0() {
        //          012345678901234
        String s = "((-12+3)*2-5*2)";
        AlgebraicExpressionEvaluator eval = new AlgebraicExpressionEvaluator();
        double r = eval.evaluateAlgebraicExpression(s);
        System.out.printf("expression=%s, algebraic evaluation=%.3f\n", s, r);
        assertEquals(-28.0, r);

        s = "((-12+3)*2-5*2^4)";
        r = eval.evaluateAlgebraicExpression(s);
        System.out.printf("expression=%s, algebraic evaluation=%.3f\n", s, r);
        assertEquals(-98.0, r);

        s = "1";
        r = eval.evaluateAlgebraicExpression(s);
        System.out.printf("expression=%s, algebraic evaluation=%.3f\n", s, r);
        assertEquals(1.0, r);
    }
}
