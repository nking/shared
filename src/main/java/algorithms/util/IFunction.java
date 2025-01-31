package algorithms.util;

/**
 * interface for classes providing an implementation
 * of a function and its derivative.
 * 
 * @author nichole
 */
public interface IFunction {
    
    /** evaluates the objective at this set of coefficients.  NOTE that any data
     * specific to a function, such as input multi-dimensional observations,
     * can be part of the concrete implementation  instance variables.
     @param coeffs
     @return evaluation 
    */
    double f (double[] coeffs);
    
    /** estimates the gradient
     @param coeffs
     @return gradient estimate
    */
    double[] der(double[] coeffs);
}
