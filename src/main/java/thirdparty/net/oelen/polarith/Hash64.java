/*
The author of this code is Wilco Oelen and he offers it
freely without copyright, but asks that his pages are referenced
as the source if used.
He has a webpage with information on the polynomial software
he ported and more modern versions which require jini bindings:
https://woelen.homescience.net/science/math/exps/polynomials/
https://woelen.homescience.net/science/math/exps/polynomials/software.html
The code here is from the Java port of RPoly, CPoly and MPSolve 1996 algorithms:
https://woelen.homescience.net/science/math/exps/polynomials/software/polsolve.tgz
*/
package thirdparty.net.oelen.polarith;

import java.io.UnsupportedEncodingException;

// A strong very fast 64 byte hash generator.
// This can be used on any datatype which can
// be represented as a byte array. A helper
// function for strings is supplied.

class Hash64 {
    private static final long[] byteTable;
    static {
        byteTable = new long[256];
        long h = 0x544B2FBACAAF1684L;
        for (int i = 0; i < 256; i++) {
            for (int j = 0; j < 31; j++) {
                h = (h >>> 7) ^ h;
                h = (h << 11) ^ h;
                h = (h >>> 10) ^ h;
            }
            byteTable[i] = h;
        }
    }
    
    private static final long HSTART = 0xBB40E64DA205B064L;
    private static final long HMULT = 7664345821815920749L;
    
    
    public static long hash(String msg) {
        long hsh;
        byte[] data = null;
        try {
            data = msg.getBytes("UTF-8");
            hsh = hash(data);
        } 
        catch (UnsupportedEncodingException ign) {
            // This never occurs, the UTF-8 encoding is a standard
            // encoding and this exception hence never is thrown.
            hsh = msg.hashCode();
        }
        return hsh;
    }
    
    
    
    public static long hash(byte[] data) {
        long h = HSTART;
        final long hmult = HMULT;
        final long[] ht = byteTable;
        for (int len = data.length, i = 0; i < len; i++) {
            h = (h * hmult) ^ ht[data[i] & 0xff];
        }
        return h;
    }
}
