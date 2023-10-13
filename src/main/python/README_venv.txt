to create a venv:
   python3 -m venv venv
   (the later can be named anything, but venv is standard)

to activate venv:
   source venv/bin/activate
   (or equiv for other than unix)

to confirm that venv is activated:
   which python

to deactivate:
   deactivate

to import packages into this venv after it has been activated:
   pip3 install numpy
   pip3 install tensorflow
   pip3 install matplotlib

use python3 here (not python) due to tensorflow dependency and
many other packages

to run the code from command line and keep the python shell open
for interactive use:
   python3 -i <script>
