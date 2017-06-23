package algorithms.util;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.Enumeration;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.logging.Level;
import java.util.logging.Logger;
import junit.framework.TestCase;

/**
 *
 * @author nichole
 */
public class ResourceFinderTest extends TestCase {
    
    public ResourceFinderTest(String testName) {
        super(testName);
    }
    
    public void testJarResource() throws IOException {
        
        String sep = System.getProperty("file.separator");
        String cwd = System.getProperty("user.dir");
          
        String jarFileName = "trove4j-3.0.3.jar";
        jarFileName = cwd + sep + "lib" + sep + jarFileName; 
        
        InputStream inStream = ResourceFinder.findJarEntry(
            jarFileName, 
            "META-INF/MANIFEST.MF");
        
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        int c;
        while ((c = inStream.read()) != -1) {
            out.write(c);
        }
        inStream.close();
        String contents = out.toString();
        out.close();

        assertNotNull(contents);
        
        assertTrue(contents.length() > 0);
    }

    
}
