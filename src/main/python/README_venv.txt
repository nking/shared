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

Misc other:

to run the code from command line and keep the python shell open
for interactive use:
   python3 -i <script>

to make a jupytper notebook:
   pip3 install jupyter
   jupyter notebook

for easier use of open-images, can install fiftyone
   pip3 install fiftyone
   https://docs.voxel51.com/getting_started/install.html

interactive fiftyone plots for use in jupyter notebooks
   pip3 install ipywidgets>=7.5

pip3 install pycocotools
