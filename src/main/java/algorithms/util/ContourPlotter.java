package algorithms.util;

import gnu.trove.iterator.TLongFloatIterator;
import gnu.trove.map.TLongFloatMap;
import java.io.IOException;

/**
 * uses the d3.js library and contour example to make an
 * html contour plot.
 * 
 * https://github.com/d3/d3-contour
 * which has a BSD 3-clause "New" or "Revised" License
 * https://github.com/d3/d3/blob/master/LICENSE
 * 
 * https://bl.ocks.org/mbostock/7f5f22524bd1d824dd53c535eda0187f
 * is Released under the GNU General Public License, version 3.
 * 
 * @author nichole
 */
public class ContourPlotter {
        
    private StringBuilder init() {
        StringBuilder content = new StringBuilder();
        content.append("<!DOCTYPE html>%n");
        content.append("<style>");
        content.append("%n");
        content.append(".domain { display: none;}");
        content.append("%n");
        content.append(".axis line { fill: none; stroke: black; stroke-width: 1; shape-rendering: crispEdges; }");
        content.append("%n");
        content.append("text { dominant-baseline: central; font-size: 14pt; font-weight: bold; }");
        content.append("%n");
        content.append("</style>");
        content.append("%n");
        content.append("<canvas width='960' height='960'></canvas>");
        content.append("%n");
        content.append("<svg width='1060' height='1060' style='position: absolute; left: 0; top: 0;'></svg>");
        content.append("%n");
        content.append("<script src='https://d3js.org/d3.v4.min.js'></script>");
        content.append("%n");
        content.append("<script src='https://d3js.org/d3-contour.v1.min.js'></script>");
        content.append("%n");
        content.append("<script src='https://d3js.org/d3-scale-chromatic.v1.min.js'></script>");
        content.append("%n");
        content.append("<script>");
        content.append("%n");
        return content;
    }
    
    public String writeFile(TLongFloatMap pixValueMap, int width, int height,
        String fileName) throws IOException {
        
        StringBuilder content = init();
        content.append("var values=[");
        
        TLongFloatIterator iter2 = pixValueMap.iterator();
        for (int ii = 0; ii < pixValueMap.size(); ++ii) {
            iter2.advance();
            long pixIdx = iter2.key();
            float p = iter2.value();
        
            content.append(String.format("%.3f", p));
            content.append(", ");
            if (ii > 0 && ((ii % 10) == 0)) {
                content.append("%n");
            }
        }

        content.append("];%n");
        content.append("var n=").append(width)
            .append(", m=").append(height).append("%n");
    
        addClose(content);
    
        return writeToFile(content.toString(), fileName);
    }
    
    private void addClose(StringBuilder content) {
        
        content.append("var w=n, h=m;");
        content.append("%n");
        content.append("var svg = d3.select('svg'),");
        content.append("%n");
        content.append("  width = +svg.attr('width'),");
        content.append("%n");
        content.append("  height = +svg.attr('height'),");
        content.append("%n");
        content.append("  margin = {top: 0, right: 100, bottom: 100, left: 0};");
        content.append("%n");
        content.append("var x = d3.scaleLinear().domain([0, n]).range([0, width]);");
        content.append("%n");
        content.append("y = d3.scaleLinear().domain([m, 0]).range([height - margin.bottom, margin.top]);");
        content.append("%n");
        content.append("var color = d3.scaleSequential(d3.interpolateYlGnBu)");
        content.append("%n");
        content.append("   .domain([0, 1.8]); // Points per square pixel.");
        content.append("%n");
        content.append("svg.append('g')");
        content.append("%n");
        content.append("   .attr('transform', 'translate(0,' + (height - margin.bottom) + ')')");
        content.append("%n");
        content.append("   .call(d3.axisBottom(x).ticks(null, '.1f'))");
        content.append("%n");
        content.append("   .select('.tick:last-of-type text')");
        content.append("%n");
        content.append("   .select(function() { return this.parentNode.appendChild(this.cloneNode()); })");
        content.append("%n");
        content.append("   .attr('y', -3).attr('dy', null)");
        content.append("%n");
        content.append("   .attr('font-weight', 'bold')");
        content.append("%n");
        content.append("   .text('X');%n");
        content.append("svg.append('g')%n");
        content.append("   .attr('transform', 'translate(' + 960 + ',0)')%n");
        content.append("   .call(d3.axisRight(y).ticks(null, '.1s'))%n");
        content.append("   .select('.tick:last-of-type text')%n");
        content.append("   .select(function() { return this.parentNode.appendChild(this.cloneNode()); })%n");
        content.append("   .attr('x', 3)%n");
        content.append("   .attr('text-anchor', 'start')%n");
        content.append("   .attr('font-weight', 'bold')%n");
        content.append("   .text('Y');%n");
        content.append("var canvas = document.querySelector('canvas'),%n");
        content.append("   context = canvas.getContext('2d'),%n");
        content.append("   color = d3.scaleSequential(d3.interpolateRdBu).domain([-1, 1]),%n");
        content.append("   path = d3.geoPath(null, context),%n");
        content.append("   thresholds = d3.range(-1.2, 1, 0.2),%n");
        content.append("   contours = d3.contours().size([n, m]);%n");
        content.append("context.scale((canvas.width/n), (canvas.height/m));%n");
        content.append("var dv = 0;%n");
        content.append("contours%n");
        content.append("   .thresholds(thresholds.map(function(v) { return v + dv; }))%n");
        content.append("   (values)%n");
        content.append("   .forEach(fill);%n");
        content.append("function fill(geometry) {%n");
        content.append("   context.beginPath();%n");
        content.append("   path(geometry);%n");
        content.append("   context.fillStyle = color(geometry.value);%n");
        content.append("   context.fill();%n");
        content.append("}%n");
        content.append("</script>%n");
    }
    
    protected String writeToFile(String fileContent, String fileName) throws IOException {

        return ResourceFinder.writeToCWD(fileContent, fileName);
    }

}
