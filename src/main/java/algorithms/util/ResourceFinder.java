package algorithms.util;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.net.URL;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
   first implemented in projects
     https://github.com/nking/two-point-correlation
     w/ Copyright (c) 2013-2015 Nichole King
     http://nking.github.io/two-point-correlation/
     using The MIT License (MIT)
     and
     https://github.com/nking/curvature-scale-space-corners-and-transformations
     w/ Copyright (c) 2014 Climb With Your Feet
     and using The MIT License (MIT)

   then moved to this shared library project which has the same copyright
   and license.

 * @author nichole
 */
public class ResourceFinder {

    /**
     *
     */
    public final static String sep = System.getProperty("file.separator");

    /**
     * find the jar file then the entry for the filePath within it.
     * 
     * NOTE: to read the file into a string, can use:
     * <pre>
            inStream = findJarEntry(...);
            ByteArrayOutputStream out = new ByteArrayOutputStream();
            int c;
            while ((c = inStream.read()) != -1) {
                out.write(c);
            }
            inStream.close();
            String contents = out.toString();
            out.close();
     * </pre>
     * It is the invoker's responsibility to close the stream.
     @param jarPath
     @param filePath
     @return
     * @throws IOException 
     */
    public static InputStream findJarEntry(String jarPath, String filePath) throws IOException {
        
        //String cwd = System.getProperty("user.dir");
        //String jarFileName = "trove4j-3.0.3.jar";
        //jarFileName = cwd + sep + "lib" + sep + jarFileName; 
        
        File f = new File(jarPath);
        if (!f.exists()) {
            throw new IOException("could not find jar at " + jarPath);
        }
        
        JarFile jarFile = new JarFile(jarPath);
        
        /*for (Enumeration<JarEntry> em = jarFile.entries(); 
            em.hasMoreElements();) {
            String s= em.nextElement().toString();
            System.out.println("FILE=" + s);
        }*/
        
        JarEntry je = jarFile.getJarEntry(filePath);

        return jarFile.getInputStream(je);
    }
    
    /**
     *
     @param fileName
     @return
     * @throws IOException
     */
    public static String findFileInResources(String fileName) throws IOException {

        String dirPath = findResourcesDirectory();

        String filePath = dirPath + sep + fileName;

        File f = new File(filePath);
        if (!f.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return filePath;
    }

    /**
     *
     @return
     * @throws IOException
     */
    public static String findOutputTestDirectory() throws IOException {

        String binDir = findDirectory("bin");
        
        String testDir = binDir + sep + "test-classes";
        
        File f = new File(testDir);
        if (!f.exists()) {
            throw new IOException("could not directory bin/test-classes");
        }
        return testDir;
    }
    
    /**
     *
     @return
     * @throws IOException
     */
    public static String findResourcesDirectory() throws IOException {
        try {
            return findDirectory("resources");
        } catch (IOException e) {
            return findDirectory(new String[]{"src", "main"}, "resources");
        }
    }
    
    /**
     *
     @return
     * @throws IOException
     */
    public static String findTestResourcesDirectory() throws IOException {
        try {
            return findDirectory("testresources");
        } catch (IOException e) {
            return findDirectory(new String[]{"src", "test"}, "resources");
        }
    }

    /**
     *
     @param dirName
     @return
     * @throws IOException
     */
    public static String findDirectory(String dirName) throws IOException {

        String cwd = getBaseDir();
        
        assert(cwd != null && !cwd.equals(""));
        
        String filePath = cwd + sep + dirName;

        File f = new File(filePath);
        
        //System.out.println("looking for " + f.toString() + " exists=" + f.exists());
        
        if (!f.exists()) {
            ClassLoader cls = ResourceFinder.class.getClassLoader();
            URL url = cls.getResource(dirName);
            if (url == null) {
                throw new IOException("could not find directory named " + dirName);
            }
            f = new File(url.getPath());
            if (!f.exists()) {
                throw new IOException("could not find directory named " + dirName);
            }
            //System.out.println("   found in classloader");
        }
        return filePath;
    }
    
    /**
     *
     @param path
     @param dirName
     @return
     * @throws IOException
     */
    public static String findDirectory(String[] path, String dirName) throws IOException {

        String dir = getBaseDir();
        
        assert(dir != null && !dir.equals(""));

        for (String p : path) {
            dir = dir + sep + p;
        }

        String filePath = dir + sep + dirName;

        File f = new File(filePath);
        if (!f.exists()) {
            ClassLoader cls = ResourceFinder.class.getClassLoader();
            URL url = cls.getResource(dirName);
            if (url == null) {
                throw new IOException("could not find directory named " + dirName);
            }
            f = new File(url.getPath());
            if (!f.exists()) {
                throw new IOException("could not find directory named " + dirName);
            }
        }
        return filePath;
    }

    /**
     *
     @param fileName
     @return
     * @throws IOException
     */
    public static String findFileInTestResources(String fileName) throws IOException {

        try {

            String dirPath = findTestResourcesDirectory();//findDirectory("testresources");
            String filePath = dirPath + sep + fileName;

            File f = new File(filePath);
            
            //System.out.println("looking for " + f.toString() + " exists=" + f.exists());

            if (!f.exists()) {
                throw new IOException("could not find file at " + filePath);
            }
            return filePath;

        } catch (IOException e) {

            String dirPath = 
                findDirectory(new String[]{"src", "test"}, "resources");
            
            String filePath = dirPath + sep + fileName;

            File f = new File(filePath);
            if (!f.exists()) {
                throw new IOException("could not find file at " + filePath);
            }
            return filePath;
        }
    }

    /**
     *
     @param serializationFileName
     @return
     * @throws IOException
     */
    public static String findFileInCWD(String serializationFileName) throws IOException {

        String filePath = getCWD() + sep + serializationFileName;

        return filePath;
    }

    /**
     * get the base directory of the project.  note that the string does not end with
     * a file separator "/".
     * @return
     * @throws IOException
     */
    public static String getCWD() throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL basedir = cls.getResource(".");
        if (basedir == null) {
            throw new IOException("base path not found");
        }

        String filePath = basedir.getPath();

        return filePath;
    }

    /**
     * get the base directory of the project.  note that the string does not end with
     * a file separator "/".
     * @return
     * @throws IOException
     */
    public static String getBaseDir() throws IOException {

        ClassLoader cls = ResourceFinder.class.getClassLoader();

        URL basedir = cls.getResource(".");
        if (basedir == null) {
            throw new IOException("base path not found");
        }

        String filePath = basedir.getPath();

        String[] rep = new String[]{"/bin/classes/", "/bin/test-classes/"};
        for (String r : rep) {
            if (filePath.endsWith(r)) {
                filePath = filePath.replace(r, "");
            }
        }

        return filePath;
    }

    /**
     *
     @param fileName
     @return
     * @throws IOException
     */
    public static String getAFilePathInCWD(String fileName) throws IOException {

        String filePath = getCWD() + sep + fileName;

        return filePath;
    }

    /**
     *
     @return
     * @throws IOException
     */
    public static String findTmpDataDirectory() throws IOException {

        return findDirectory("tmpdata");
    }

    /**
     *
     @param fileName
     @return
     * @throws IOException
     */
    public static File findFileInTmpData(String fileName) throws IOException {

        String filePath = getAFilePathInTmpData(fileName);

        File fl = new File(filePath);
        if (!fl.exists()) {
            throw new IOException("could not find file at " + filePath);
        }
        return fl;
    }

    /**
     *
     @param fileName
     @return
     * @throws IOException
     */
    public static String getAFilePathInTmpData(String fileName) throws IOException {

        String baseDir = findTmpDataDirectory();

        String filePath = baseDir + sep + fileName;

        return filePath;
    }

    /**
     *
     @param fileContent
     @param fileName
     @return
     * @throws IOException
     */
    public static String writeToCWD(String fileContent, String fileName) throws IOException {

        String filePath = getAFilePathInCWD(fileName);

        return writeDataToDirectory(fileContent, filePath);
    }

    /**
     *
     @param fileContent
     @param filePath
     @return
     * @throws IOException
     */
    protected static String writeDataToDirectory(String fileContent, String filePath) throws IOException {

        FileWriter fw = null;
        BufferedWriter writer = null;
        try {
            File file = new File(filePath);
            if (file.exists()) {
                boolean deleted = file.delete();
                if (!deleted) {
                    throw new IOException("could not delete file " + file.getPath());
                }
            }
            boolean dreated = file.createNewFile();
            if (!dreated) {
                throw new IOException("could not create new file " + file.getPath());
            }

            fw = new FileWriter(file);
            writer = new BufferedWriter(fw);
            writer.write(fileContent);

            writer.flush();

        } finally {

            if (writer != null) {
                writer.close();
            }
            if (fw != null) {
                fw.close();
            }
        }
        return filePath;
    }

}
