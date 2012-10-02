Convert DICOM-RT Structures to a PLY Point Cloud
================================================

This is a very simple utility to convert any structures found in a DICOM-RT data set to individual PLY files for subsequent surface meshing. This will dump a bunch of PLY files into the current directory:

    python dicom2ply.py dicom_dir

Originally this code was written to quickly calculate descriptive statistics for each structure; some of this information is included in the PLY header.

The saved point clouds can be quite dense, so in some instances it is best to reduce the number of total points prior to meshing. Have a look at [MeshLab](http://meshlab.sourceforge.net/) for meshing and mesh refinement.

Dependencies
------------
* [pydicom](http://code.google.com/p/pydicom/)
* [Python Imaging Library](http://www.pythonware.com/products/pil/)
* numpy, scipy, pylab

TODO's
------
* convert into a proper Python module
* tidy up the statistics calculations
* remove hard coding of CT data file suffixes and prefixes 
