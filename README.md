# Bioformats Image Toolbox

This toolbox re-implements the OME Microscopy MATLAB toolbox as a MATLAB object. It is intended to be compatible with the MATLAB Add-On Manager to be as easy to use as possible.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for _development_. If you are looking for a released build for use in MATLAB, please have a look at the Wiki instead.

### Prerequisites

  - MATLAB R2016b or higher (although it is highly likely that an earlier version will work).
  - [A copy of Git](https://git-scm.com/downloads)

### Downloading a copy of the code

  git clone git@biof-git.colorado.edu:core-code/bioformats-image-toolbox.git


### Toolbox directory structure

The toolbox is structured following the Best Practices outlined in this [Mathworks Blogs post](http://blogs.mathworks.com/developer/2017/01/13/matlab-toolbox-best-practices/).

The main toolbox files are in `/tbx/bfimage`. Any supporting documentation (including MATLAB examples) are in `/tbx/docs/`.

The `/local` folder is used to hold data that will not be synced to the repository (i.e. ND2 files for testing).

### Using the toolbox

Add the `/tbx` folder along with all its subfolders to MATLAb's search path by running `addsandbox.m`. If for some reason this fails, you can do this manually by right clicking on `/tbx` → `Add to Path` → `Selected Folder and Subfolders` in the 'Current Folder' window.

Check that the toolbox is running by running

  bfReader = BioformatsImage;

If you do not have the OME Microscopy MATLAB toolbox installed, the object will download it on first run.
  
## Building the toolbox

1. Make sure that the `/tbx` folder and all its subfolders are added to the MATLAB path. The easiest way to do this is to run `addsandbox.m`.

2. Load `BioformatsImage.prj` in MATLAB's Toolbox Packager tool.

3. Check that the folders `bfimage` and `docs` are shown under 'Toolbox Files and Folders'.

4. Under 'Install Actions', make sure that the following folders are present in the 'MATLAB Path' frame:

        <Toolbox Folder>
        <Toolbox Folder>/bfimage
        <Toolbox Folder>/docs

    If these folders are missing, check that the `/tbx` folder and its subfolders are added to the path, then click on 'Reset to the current MATLAB path'.

5. Click on the 'Package' button on the toolbar.

6. To test the package, remove the `/tbx` folder and all its subfolders by running `rmsandbox.m`. Then run the tests described in the section above.

## Running tests

To ensure compatibility, there are a number of tests that are included in the folder.

To run these tests, use

    runtests('testing')
  
This project also uses [tags](https://www.mathworks.com/help/matlab/matlab_prog/tag-unit-tests.html) to group tests into different categories. For more details, consult the Developer's Guide of the Wiki.

> Please note: If you are intending to contribute code to the project, all tests have to be passed before the merge request can be approved.

## Contribution guide

To report a bug, use the project's [Gitlab issue tracker](https://biof-git.colorado.edu/core-code/bioformats-image-toolbox/issues). 

Note: you may need to have be connected to the internet on-campus internet access or via the VPN. For information about connecting to the VPN, please consult [OIT's website](https://oit.colorado.edu/services/network-internet-services/vpn).

## Attribution

This project uses the [OME Bio-Formats Library](https://www.openmicroscopy.org/site/products/bio-formats).
