code containing algorithms used by more than one project,
specifically by 
http://nking.github.io/two-point-correlation/
and
http://nking.github.io/curvature-scale-space-corners-and-transformations
================================================================

To build and run the project, you need to have installed java,
ant, and jacoco (which requires eclemma and provides download options).

http://java.com
https://ant.apache.org/bindownload.cgi
http://www.jacoco.org/jacoco/
and https://www.jacoco.org/installation.html#manual

The project requires java 1.7 or greater.
The ant version should be 1.9.6 or greater.
The jacoco version should be jacoco-0.7.5.201505241946 or greater.
    Then set an environment variable called JACOCO_HOME
    to the path of the base directory of jacoco.
    The build script looks for $JACOCO_HOME/lib/jacocoant.jar

Note: if you choose to download jacoco and eclemma in a directory independent
of an IDE, you can simply download the latest jacoco file from
http://www.jacoco.org/jacoco/ and unpack it into a directory that you've
named JACOCO_HOME in your environment variables (e.g. defined in your path via your shell profile file, etc).
Then download the latest eclemma file from https://www.jacoco.org/download.html
into the $JACOCO_HOME directory and unpack it.
Remember to refresh your environment variables once you've updated them.

Also useful is spotbugs, the successor of findbugs:
   https://spotbugs.readthedocs.io/en/stable/installing.html

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

-----------
The project also includes in various directories, python, rust and c
code.
   The build.xml has targets for the c code.
      (The target for ISPC depends upon a few APIs that I haven't 
       written up here yet.)

   The python and rust code do not have targets in build.xml,
   but can be built from their source code directories.
    I'll add notes to building and running them here at some point.
  
------------
Dependencies
------------
If you build the jar file to use as an API
NOTE that it depends upon lib/trove4j-3.0.3.jar
lib/netlib-java-all-1.1.2.jar and lib/mtj-1.0.4.jar
which are included in the lib directory (the build
file uses that location).
