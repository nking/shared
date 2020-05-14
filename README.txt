code containing algorithms used by more than one project,
specifically by 
http://nking.github.io/two-point-correlation/
and
http://nking.github.io/curvature-scale-space-corners-and-transformations
================================================================

To build and run the project, you need to have installed java,
ant, and jacoco.

http://java.com
https://ant.apache.org/bindownload.cgi
http://www.jacoco.org/jacoco/

The project requires java 1.7 or greater.
The ant version should be 1.9.6 or greater.
The jacoco version should be jacoco-0.7.5.201505241946 or greater.
    Then set an environment variable called JACOCO_HOME
    to the path of the base direcotry of jacoco.
    The build scripts looks for $JACOCO_HOME/lib/jacocoant.jar

-----
Build
-----
To list the targets:
  ant

To compile and package :
  ant package

To compile the main source code:
  ant compile

To run all tests:
  ant runTests

To run a specific test:
  ant runTest -Dtest=package.TestName

To generate coverage reports:
  ant runCoverage
  
------------
Dependencies
------------
If you build the jar file to use as an API
NOTE that it depends upon lib/trove4j-3.0.3.jar
lib/netlib-java-all-1.1.2.jar and lib/mtj-1.0.4.jar
