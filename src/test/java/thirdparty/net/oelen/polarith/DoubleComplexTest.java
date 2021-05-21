package thirdparty.net.oelen.polarith;

// A test class, which tests the DoubleComplex class.

import java.util.Optional;
import java.util.function.BiFunction;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import thirdparty.net.oelen.polarith.DoubleComplex;
import thirdparty.net.oelen.polarith.DoubleDouble;


public class DoubleComplexTest {
    
    static final String[][] TESTS_3 = {
        {"ABS", "3", "4", "5"},
        {"ABS", "4", "3", "5"},
        {"ABS1", "3", "4", "7"},
        {"ABS1", "4", "3", "7"},
        {"ABSINF", "3", "4", "4"},
        {"ABSINF", "4", "3", "4"},
        {"SQRABS", "4", "3", "25"},
        {"REAL", "4", "3", "4"},
        {"IMAG", "4", "3", "3"},
        
    };
    
    static final String[][] TESTS_4 = {
        {"ADD1", "1", "2", "2", "2"},
        {"SUB1", "1", "2", "0", "2"},
        {"RECIP", "1", "2", "0.2", "-0.4"},
        {"RECIP", "2", "1", "0.4", "-0.2"},
        {"NEG", "1", "2", "-1", "-2"},
        {"CONJ", "1", "2", "1", "-2"},
        
    };
    
    static final String[][] TESTS_5 = {
        {"ADD", "1", "2", "-33",     "-32", "2"},
        {"ADD", "1e20", "2e20", "-33",     "99999999999999999967", "200000000000000000000"},
        
        {"SUB", "1", "2", "-33",    "34", "2"},
        
        {"MUL", "1", "2", "3",    "3", "6"},
        {"DIV", "1", "0", "4",    "0.250000000000000000000", "0"},
        {"DIV", "10e10", "2", "-4e0",   "-25.0e+9", "-50e-2"},
    };
    
    static final String[][] TESTS_6 = {
        {"ADD", "1", "2", "-33", "34", "-32", "36"},
        {"ADD", "1", "2", null, "34", "1", "36"},
        {"ADD", "1", "2", "-33", null, "-32", "2"},
        {"ADD", "1e20", "2e20", "-33", "34", "99999999999999999967", "200000000000000000034"},
        
        {"SUB", "1", "2", "-33", "34", "34", "-32"},
        {"SUB", "1", "2", null, "34", "1", "-32"},
        {"SUB", "1", "2", "-33", null, "34", "2"},
        
        {"MUL", "1", "2", "3", "4", "-5", "10"},
        {"MUL", "1", "2", null, "4", "-8", "4"},
        {"MUL", "1", "2", "3", null, "3", "6"},
        {"DIV", "1", "0", "3", "4", "0.12", "-0.16"},
        {"DIV", "1", "0", "4", "3", "0.16", "-0.12"},
        {"DIV", "1", "2", "3", "4", "0.44", "0.08"},
        {"DIV", "1", "2", "4", "3", "0.4", "0.2"},
        {"DIV", "1", "-2", "-4", null, "-0.25", "0.5"},
        {"DIV", "1", "-2", null, "-4", "0.5", "0.25"},
    };
    
    public static void main(String[] args) {
        for (String[] test : TESTS_3) {
            Optional<Function<DoubleComplex,DoubleDouble>> op = getDoubleFunction(test[0]);
            DoubleComplex a = new DoubleComplex(test[1], test[2]);
            DoubleDouble expected = new DoubleDouble(test[3]);
            op.map(f->f.apply(a))
              .filter(r->!expected.isNear(r, 1e-30))
              .ifPresent(r->System.err.println("Error " + test[0] + " " + print(a) + " ---> " + r + "//" + expected));
        }
        
        for (String[] test : TESTS_4) {
            Optional<Function<DoubleComplex,DoubleComplex>> op = getComplexFunction(test[0]);
            DoubleComplex a = new DoubleComplex(test[1], test[2]);
            DoubleComplex expected = new DoubleComplex(test[3], test[4]);
            op.map(f->f.apply(a))
              .filter(r->!expected.isNear(r, 1e-30))
              .ifPresent(r->System.err.println("Error " + test[0] + " " + print(a) + " ---> " + print(r) + "//" + print(expected)));
        }
        
        for (String[] test : TESTS_5) {
            Optional<BiFunction<DoubleComplex,DoubleDouble,DoubleComplex>> op = getComplexDoubleFunction(test[0]);
            DoubleComplex a = new DoubleComplex(test[1], test[2]);
            DoubleDouble b = new DoubleDouble(test[3]);
            DoubleComplex expected = new DoubleComplex(test[4], test[5]);
            op.map(f->f.apply(a, b))
              .filter(r->!expected.isNear(r, 1e-30))
              .ifPresent(r->System.err.println("Error " + test[0] + " " + print(a) + "   " + b + " ---> " + print(r) + "//" + print(expected)));
        }
        
        for (String[] test : TESTS_6) {
            if (test[3] == null || test[4] == null) {
                continue;
            }
            Optional<BinaryOperator<DoubleComplex>> op = getComplexOperator(test[0]);
            DoubleComplex a = new DoubleComplex(test[1], test[2]);
            DoubleComplex b = new DoubleComplex(test[3], test[4]);
            DoubleComplex expected = new DoubleComplex(test[5], test[6]);
            op.map(f->f.apply(a, b))
              .filter(r->!expected.isNear(r, 1e-30))
              .ifPresent(r->System.err.println("Error " + test[0] + " " + print(a) + "   " + print(b) + " ---> " + print(r) + "//" + print(expected)));
        }
        
        for (String[] test : TESTS_6) {
            Optional<TriFunction<DoubleComplex,DoubleDouble,DoubleDouble,DoubleComplex>> op = getTriFunction(test[0]);
            DoubleComplex a = new DoubleComplex(test[1], test[2]);
            DoubleDouble bre = test[3] != null ? new DoubleDouble(test[3]) : null;
            DoubleDouble bim = test[4] != null ? new DoubleDouble(test[4]) : null;
            DoubleComplex expected = new DoubleComplex(test[5], test[6]);
            op.map(f->f.apply(a, bre, bim))
              .filter(r->!expected.isNear(r, 1e-30))
              .ifPresent(r->System.err.println("Error " + test[0] + " " + print(a) + "   " + bre + "||" + bim  + " ---> " + print(r) + "//" + print(expected)));
        }
        
        for (int i=-50; i<=50; i++) {
            if (i == 0) continue;
            DoubleDouble x = new DoubleDouble(i);
            DoubleDouble y = x.root(121);
            DoubleDouble diff = x.sub(y.pow(121));
            System.out.println("  " + i + "   " + diff.div(x));
        }
        System.out.println();
        for (int i=1; i<=100; i++) {
            DoubleDouble x = new DoubleDouble(i);
            DoubleDouble diff = x.log().exp().sub(x);
            System.out.println("  " + i + "   " + diff.div(x));
        }
    }
    
    
    
    
    private static Optional<Function<DoubleComplex, DoubleDouble>> getDoubleFunction(String f) {
        switch (f) {
            case "ABS":
                return Optional.of(DoubleComplex::abs);
            case "ABS1":
                return Optional.of(DoubleComplex::abs1);
            case "ABSINF":
                return Optional.of(DoubleComplex::absinf);
            case "SQRABS":
                return Optional.of(DoubleComplex::sqrabs);
            case "REAL":
                return Optional.of(DoubleComplex::real);
            case "IMAG":
                return Optional.of(DoubleComplex::imag);
            default:
                return Optional.empty();
        }
    }
    
    
    
    
    private static Optional<Function<DoubleComplex, DoubleComplex>> getComplexFunction(String f) {
        switch (f) {
            case "ADD1":
                return Optional.of(DoubleComplex::add1);
            case "SUB1":
                return Optional.of(DoubleComplex::sub1);
            case "RECIP":
                return Optional.of(DoubleComplex::recip);
            case "NEG":
                return Optional.of(DoubleComplex::neg);
            case "CONJ":
                return Optional.of(DoubleComplex::conj);
            default:
                return Optional.empty();
        }
    }
    
    
    
    
    private static Optional<BiFunction<DoubleComplex,DoubleDouble,DoubleComplex>> getComplexDoubleFunction(String f) {
        switch (f) {
            case "ADD":
                return Optional.of(DoubleComplex::add);
            case "SUB":
                return Optional.of(DoubleComplex::sub);
            case "MUL":
                return Optional.of(DoubleComplex::mul);
            case "DIV":
                return Optional.of(DoubleComplex::div);
            default:
                return Optional.empty();
        }
    }
    
    
    
    
    private static Optional<BinaryOperator<DoubleComplex>> getComplexOperator(String f) {
        switch (f) {
            case "ADD":
                return Optional.of(DoubleComplex::add);
            case "SUB":
                return Optional.of(DoubleComplex::sub);
            case "MUL":
                return Optional.of(DoubleComplex::mul);
            case "DIV":
                return Optional.of(DoubleComplex::div);
            default:
                return Optional.empty();
        }
    }
    
    
    
    
    private static Optional<TriFunction<DoubleComplex,DoubleDouble,DoubleDouble,DoubleComplex>> getTriFunction(String f) {
        switch (f) {
            case "ADD":
                return Optional.of(DoubleComplex::add);
            case "SUB":
                return Optional.of(DoubleComplex::sub);
            case "MUL":
                return Optional.of(DoubleComplex::mul);
            case "DIV":
                return Optional.of(DoubleComplex::div);
            default:
                return Optional.empty();
        }
    }
    
    
    
    private static String print(DoubleComplex a) {
        return "" + a.real().toString() + "||" + a.imag().toString();
    }    
}



interface TriFunction<T,U,V,R> {
    R apply(T t, U u, V v);
}
