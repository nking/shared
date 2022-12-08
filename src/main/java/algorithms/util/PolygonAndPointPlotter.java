package algorithms.util;

import algorithms.misc.Misc0;
import algorithms.misc.MiscMath0;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;

/**
 * writes an html file using the d3.js library.
 * https://github.com/d3/
 * which has a BSD 3-clause "New" or "Revised" License
 * https://github.com/d3/d3/blob/master/LICENSE
 * 
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
public class PolygonAndPointPlotter {

    /**
     *
     */
    protected final StringBuffer plotContent;

    /**
     *
     */
    protected int plotNumber = 0;

    /**
     *
     */
    protected boolean dataMinMaxAreSet = false;
    
    /**
     *
     */
    protected float minX,

    /**
     *
     */
    maxX, minY,

    /**
     *
     */
    maxY;

    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     * @throws IOException
     */
    public PolygonAndPointPlotter(float minX, float maxX, float minY, 
        float maxY) throws IOException {

        plotContent = getTemplateHtmlPlot();

        setDataMinMax(plotContent, minX, maxX, minY, maxY);
    }

    /**
     *
     * @throws IOException
     */
    public PolygonAndPointPlotter() throws IOException {

        plotContent = getTemplateHtmlPlot();
    }

    /**
     *
     @param plotContent
     @param minX
     @param maxX
     @param minY
     @param maxY
     */
    protected void setDataMinMax(StringBuffer plotContent, 
        float minX, float maxX, float minY, float maxY) {

        this.minX = minX;
        this.maxX = maxX;
        this.minY = minY;
        this.maxY = maxY;
        
        StringBuilder dataSB = new StringBuilder();

        dataSB.append("\nvar xmin = ").append(minX).append(";\n");
        dataSB.append("var xmax = ").append(maxX).append(";\n");
        dataSB.append("var ymin = ").append(minY).append(";\n");
        dataSB.append("var ymax = ").append(maxY).append(";\n");

        String srchFor = "/* === DO NOT REMOVE THIS == START DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        insertOffset += srchFor.length();
        plotContent.insert(insertOffset, dataSB.toString());

        dataMinMaxAreSet = true;
    }

    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float minX, float maxX, float minY, float maxY,
        float[] xPoints, int[] yPoints, float[] xPolygon, int[] yPolygon, 
        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float minX, float maxX, float minY, float maxY,
        int[] xPoints, float[] yPoints, int[] xPolygon, float[] yPolygon, 
        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float minX, float maxX, float minY, float maxY,
        int[] xPoints, int[] yPoints, int[] xPolygon, int[] yPolygon, 
        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float minX, float maxX, float minY, float maxY,
        float[] xPoints, float[] yPoints, float[] xPolygon, float[] yPolygon, 
        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }

    /**
     *
     @param minX
     @param maxX
     @param minY
     @param maxY
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(double minX, double maxX, double minY, double maxY,
                        double[] xPoints, double[] yPoints, double[] xPolygon, double[] yPolygon,
                        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }

    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float[] xPoints, int[] yPoints, float[] xPolygon, 
        int[] yPolygon, String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }

    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(int[] xPoints, float[] yPoints, int[] xPolygon, 
        float[] yPolygon, String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }

    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float[] xPoints, float[] yPoints, float[] xPolygon, 
        float[] yPolygon, String plotLabel) {
        
        addPlot(xPoints, yPoints, null, null, xPolygon, yPolygon, plotLabel);
    }
    
    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(int[] xPoints, int[] yPoints, int[] xPolygon, 
        int[] yPolygon, String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float[] xPoints, float[] yPoints, int[] xPolygon, 
        int[] yPolygon, String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param xPoints
     @param yPoints
     @param xLinePairs
     @param yLinePairs
     @param plotLabel
     */
    public void addPlotWithLines(float[] xPoints, float[] yPoints, 
        float[] xLinePairs, float[] yLinePairs, 
        String plotLabel) {
        
        float[] xErrPoints = null; 
        float[] yErrPoints = null;

        addPlot(xPoints, yPoints, xErrPoints, yErrPoints,
            xLinePairs, yLinePairs, plotLabel);
        
        String srch = "renderPlot('";
        int srchIdx = this.plotContent.lastIndexOf(srch);
        plotContent.insert(srchIdx + srch.length() - 2, 
            "WithLines");
    }

    /**
     *
     @param xPoints
     @param yPoints
     @param xErrPoints
     @param yErrPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(double[] xPoints, double[] yPoints,
                        double[] xErrPoints, double[] yErrPoints,
                        double[] xPolygon, double[] yPolygon,
                        String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                Misc0.convertToNumberArray(xErrPoints), Misc0.convertToNumberArray(yErrPoints),
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
    }
    
    /**
     *
     @param xPoints
     @param yPoints
     @param xErrPoints
     @param yErrPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(float[] xPoints, float[] yPoints, 
        float[] xErrPoints, float[] yErrPoints, 
        float[] xPolygon, float[] yPolygon, 
        String plotLabel) {

        if (!dataMinMaxAreSet) {

            float minX0 = MiscMath0.findMin(xPoints);
            float maxX0 = MiscMath0.findMax(xPoints);
            float minY0 = MiscMath0.findMin(yPoints);
            float maxY0 = MiscMath0.findMax(yPoints);

            addPlot(minX0, maxX0, minY0, maxY0,
                    Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                    Misc0.convertToNumberArray(xErrPoints), Misc0.convertToNumberArray(yErrPoints),
                    Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);

        } else {

            addPlot(minX, maxX, minY, maxY,
                    Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                    Misc0.convertToNumberArray(xErrPoints), Misc0.convertToNumberArray(yErrPoints),
                    Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);
        }
    }

    /**
     *
     @param xmn
     @param xmx
     @param ymn
     @param ymx
     @param xPoints
     @param yPoints
     @param xErrPoints
     @param yErrPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(Number xmn, Number xmx, Number ymn, Number ymx,
                        Number[] xPoints, Number[] yPoints,
                        Number[] xErrPoints, Number[] yErrPoints,
                        Number[] xPolygon, Number[] yPolygon, String plotLabel) {

        if (!dataMinMaxAreSet) {
            xmn = MiscMath0.findMin(xPoints);
            xmx = MiscMath0.findMax(xPoints);
            ymn = MiscMath0.findMin(yPoints);
            ymx = MiscMath0.findMax(yPoints);
            // do not set dataMinMaxAreSet to true to let the next plot scale as needed

            if (xmn.doubleValue() < 0) {
                xmn = 1.1 * xmn.doubleValue();
            } else {
                xmn = 0.9 * xmn.doubleValue();
            }
            if (xmx.doubleValue() < 0) {
                xmx = 0.9 * xmx.doubleValue();
            } else {
                xmx = 1.1 * xmx.doubleValue();
            }

            if (ymn.doubleValue() < 0) {
                ymn = 1.1 * ymn.doubleValue();
            } else {
                ymn = 0.9 * ymn.doubleValue();
            }
            if (ymx.doubleValue() < 0) {
                ymx = 0.9 * ymx.doubleValue();
            } else {
                ymx = 1.1 * ymx.doubleValue();
            }
        }

        StringBuffer dataSB = new StringBuffer("\n");

        //  ===== add plotLabel data =====
        dataSB.append("var data_plot_label_").append(plotNumber)
                .append(" = '").append(plotLabel).append("';\n");

        //  ===== add points data =====
        if (xPoints == null) {
            dataSB.append("var data_points_").append(plotNumber)
                    .append(" = undefined;\n");
        } else {
            dataSB.append("var data_points_").append(plotNumber).append(" = [\n");
            for (int i = 0; i < xPoints.length; i++) {
                if (i > 0) {
                    dataSB.append(",\n");
                }
                dataSB.append("    {x:").append(xPoints[i]).append(", y:").append(yPoints[i]);
                if (xErrPoints != null && xErrPoints.length> 0) {
                    dataSB.append(", dx:").append(xErrPoints[i]).append(", dy:").append(yErrPoints[i]);
                }
                dataSB.append("}");
            }
            dataSB.append("\n];\n");
        }

        if (xPolygon == null) {
            dataSB.append("var data_polygon_").append(plotNumber)
                    .append(" = undefined;\n");
        } else {
            //  ===== add polygon =====
            dataSB.append("var data_polygon_").append(plotNumber).append(" = [\n");
            dataSB.append("    [");
            for (int ii = 0; ii < xPolygon.length; ii++) {
                String xStr = String.format("%.7f", xPolygon[ii]);
                String yStr = String.format("%.7f", yPolygon[ii]);
                if (ii > 0) {
                    dataSB.append(", ");
                }
                dataSB.append("    {x:").append(xStr)
                        .append(", y:").append(yStr).append("}");
            }
            dataSB.append("],\n ");
            dataSB.append("];\n");
        }

        dataSB.append("var xmin_").append(plotNumber).append("=")
                .append(xmn).append(";\n");
        dataSB.append("var xmax_").append(plotNumber).append("=")
                .append(xmx).append(";\n");
        dataSB.append("var ymin_").append(plotNumber).append("=")
                .append(ymn).append(";\n");
        dataSB.append("var ymax_").append(plotNumber).append("=")
                .append(ymx).append(";\n");

        // ======= add RENDER statement ==========
        dataSB.append("\nrenderPlot('plot").append(plotNumber)
                .append("', data_points_").append(plotNumber)
                .append(", data_polygon_").append(plotNumber)
                .append(", data_plot_label_").append(plotNumber)
                .append(", ")
                .append(" xmin_").append(plotNumber).append(", ")
                .append(" xmax_").append(plotNumber).append(", ")
                .append(" ymin_").append(plotNumber).append(", ")
                .append(" ymax_").append(plotNumber)
                .append( ");\n\n");

        String srchFor = "/* === DO NOT REMOVE THIS == END DATA */";
        int insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, dataSB.toString());
        dataSB = null;


        // ========== add the PLOT DIVS ==============
        StringBuffer plotDivs = new StringBuffer();
        plotDivs.append("<div id='plot").append(plotNumber).append("' class='plot'></div>\n");


        srchFor = "<!-- === DO NOT REMOVE THIS == END PLOT DIVS -->";
        insertOffset = plotContent.indexOf(srchFor);
        if (insertOffset == -1) {
            throw new IllegalStateException("cannot find END DATA marker");
        }
        plotContent.insert(insertOffset, plotDivs.toString());
        plotDivs = null;

        plotNumber++;
    }

    /**
     *
     @param xPoints
     @param yPoints
     @param xPolygon
     @param yPolygon
     @param plotLabel
     */
    public void addPlot(double[] xPoints, double[] yPoints, double[] xPolygon, 
        double[] yPolygon, String plotLabel) {

        addPlot(minX, maxX, minY, maxY,
                Misc0.convertToNumberArray(xPoints), Misc0.convertToNumberArray(yPoints),
                null, null,
                Misc0.convertToNumberArray(xPolygon), Misc0.convertToNumberArray(yPolygon), plotLabel);

    }
    
    /**
     *
     @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected final StringBuffer getTemplateHtmlPlot() throws FileNotFoundException, IOException {
        return getTemplateHtmlPlot("plot_points_and_polygon.html");
    }

    /**
     *
     @param fileName
     @return
     * @throws FileNotFoundException
     * @throws IOException
     */
    protected final StringBuffer getTemplateHtmlPlot(String fileName)
        throws FileNotFoundException, IOException 
    {
            
        try {
            
            String path = ResourceFinder.findFileInResources(fileName);
        
            StringBuffer sb = new StringBuffer();
            BufferedReader in = null;
            try {
                in = new BufferedReader(new FileReader(new File(path)));
                String line = in.readLine();
                while (line != null) {
                    sb.append(line).append("\n");
                    line = in.readLine();
                }
            } finally {
                if (in != null) {
                    in.close();
                }
            }
            return sb;
            
        } catch (IOException ex) {
            
            // this class and resources might be in a jar file, so look
            // for that
            
            String sep = System.getProperty("file.separator");
            String cwd = System.getProperty("user.dir");
                        
            String jarFilePath = "com.climbwithyourfeet.shared.jar";
            jarFilePath = cwd + sep + "lib" + sep + jarFilePath; 
          
            InputStream inStream = null;
            ByteArrayOutputStream out = null;
            
            try {
                inStream = ResourceFinder.findJarEntry(jarFilePath, fileName);
                out = new ByteArrayOutputStream();
                int c;
                while ((c = inStream.read()) != -1) {
                    out.write(c);
                }
                StringBuffer contents = new StringBuffer(out.toString());

                return contents;
                
            } finally {
                if (inStream != null) {
                    inStream.close();
                }
                if (out != null) {
                    out.close();
                }
            }            
        }

        
    }

    /**
     *
     @return
     * @throws IOException
     */
    public String writeFile() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon.html");
    }

    /**
     *
     @return
     * @throws IOException
     */
    public String writeFile2() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon2.html");
    }

    /**
     *
     @return
     * @throws IOException
     */
    public String writeFile3() throws IOException {
        return writeToFile(this.plotContent.toString(), "points_and_polygon3.html");
    }
    
    /**
     *
     @param num
     @return
     * @throws IOException
     */
    public String writeFile(Integer num) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon" + num.toString() +".html");
    }
    
    /**
     *
     @param fileSuffix
     @return
     * @throws IOException
     */
    public String writeFile(String fileSuffix) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon_" + fileSuffix + ".html");
    }
    
    /**
     *
     @param num
     @return
     * @throws IOException
     */
    public String writeFile(long num) throws IOException {
        return writeToFile(this.plotContent.toString(), 
            "points_and_polygon" + Long.toString(num) +".html");
    }

    /**
     *
     @param fileContent
     @param fileName
     @return
     * @throws IOException
     */
    protected String writeToFile(String fileContent, String fileName) throws IOException {

        return ResourceFinder.writeToCWD(fileContent, fileName);
    }

}
