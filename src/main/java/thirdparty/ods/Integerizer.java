package thirdparty.ods;

/*
The class is from the open datastructures source code
http://opendatastructures.org/ods-java.pdf

"The source code available there is released under a Creative Commons
Attribution license, meaning that anyone is free to share: to copy, distribute
and transmit the work; and to remix: to adapt the work, including
the right to make commercial use of the work. The only condition on
these rights is attribution: you must acknowledge that the derived work
contains code and/or text from opendatastructures.org.
https://github.com/patmorin/ods

Edits were made to the code in project
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet, The MIT License (MIT)
and then moved to this project

*/

/**
 *
 * @author nichole
 @param <T>
 */


public interface Integerizer<T> {

    /**
     *
     @param x
     @return
     */
    public int intValue(T x);
}
