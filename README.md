# Bioformats Image Toolbox

This toolbox reimplements the OME Microscopy MATLAB toolbox as a MATLAB object. It is intended to be compatible with the MATLAB Add-On Manager to be as easy to use as possible.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development purposes. If you are looking for a released build for use in MATLAB, please have a look at the Wiki instead.

### Prerequisites

MATLAB R2016b or higher (although it is highly likely that an earlier version will work).

A copy of Git for your operating system

### Downloading a copy of the code


### File structure

The toolbox is structured following the Best Practices outlined in this [Mathworks Blogs post](http://blogs.mathworks.com/developer/2017/01/13/matlab-toolbox-best-practices/).

/BioformatsImage.prj 
/addsandbox.m
/rmsandbox.m
/README.md
/CHANGELOG.md
/.gitignore

/local
The local folder should be used to hold any data that will not be synced to the repository. It is ignored on .gitignore.

/local/.dummy
This file is used as a placeholder as Git does not like empty directories.

/tbx/bfimage
This is the main toolbox folder which contains all the code for the toolbox.

/tbx/bfimage/BioformatsImage.m
This is the main class file for the BioformatsImage object.

/tbx/bfimage/LookupTable.m
Holds data on the lookup table for the image.

/tbx/docs
Documentation for the toolbox is put in this directory.


### Running tests

To ensure compatibility, there are a number of tests that are included in the folder.

To run these tests

  runtests('testing')
  
To test specific parts of the code, 




### Building the toolbox

1. Make sure that the `/tbx` folder and all its subfolders are added to the MATLAB path. The easiest way to do this is to run the `addsandbox.m` script.

2. Load `BioformatsImage.prj`.

3. Check that the folders `bfimage`, `docs` and the file `Contents.m` is shown under Toolbox Files and Folders.

4. Under Install Actions, make sure that the following folders are present in the MATLAB Path frame:
  <Toolbox Folder>
  <Toolbox Folder>/bfimage
  <Toolbox Folder>/docs

If these folders are missing, check that the `/tbx` folder and its subfolders are added to the path, then click on 'Reset to the current MATLAB path'.

5. Click on the Package button on the toolbar.

6. To test the package, remove the /tbx folder and all its subfolders by running `rmsandbox.m`. Then run the tests described in the section above.

### Contribution guide

To report a bug, use the project's Gitlab issue tracker: https://biof-git.colorado.edu/core-code/bioformats-image-toolbox/issues

Note: you may need to have be connected to the internet on-campus internet access or via the VPN.

For information about connecting to the VPN, please consult [OIT](https://oit.colorado.edu/services/network-internet-services/vpn).