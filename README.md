# Bioformats Image MATLAB

Welcome to the Bioformats Image MATLAB project!

## Description

This project was started as a way to implement the OME Bio-Formats library 
as an easy-to-use MATLAB toolbox. The [OME Bio-Formats library](https://www.openmicroscopy.org/bio-formats/)
is a Java library which reads over 160 microscope file formats and metadata.

This toolbox has several main differences with the official implementation:
- The toolbox is implemented as a MATLAB object
- Image planes can be read in one at a time, thereby reducing memory and 
  processing requirements
- (_Coming soon_) Image export as OME-TIFFs

### Recent changes

v1.2.4 - Updated to OME Bio-Formats v8.1.1

## Getting started

### Prerequisites

  - MATLAB R2016b or newer (it is highly likely that an earlier version
    will work).

### Installation

1. [Download the latest version](https://github.com/Biofrontiers-ALMC/bioformats-matlab/releases) of the toolbox.
   Note that you only need to download the .mltbx file.

2. Open the .mltbx file using MATLAB. The toolbox should automatically install.

### Usage

In MATLAB:

```matlab
%Create a new BioformatsImage object
reader = BioformatsImage('path/to/file.nd2');

%Get and display the image at (Z = 1, C = 'Cy5'1, T = 5)
I = getPlane(bfr, 1, 'Cy5', 5);
imshow(I, [])
```

For more details, see the [Documentation](https://github.com/Biofrontiers-ALMC/bioformats-matlab/wiki).

## Contributing

## Top contirbutors

<a href="https://github.com/Biofrontiers-ALMC/bioformats-matlab/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=Biofrontiers-ALMC/bioformats-matlab" alt="contrib.rocks image" />
</a>

### Issues or Bugs

If you encounter issues or have questions using this toolbox, please [send 
us an email](mailto:biof-imaging@colorado.edu) or 
[create an Issue](https://github.com/Biofrontiers-ALMC/bioformats-matlab/issues).

## License

This project is licensed under the [MIT License](https://github.com/Biofrontiers-ALMC/bioformats-matlab/blob/master/LICENSE).

## Acknowledgements

This project uses the [OME Bio-Formats Library](https://www.openmicroscopy.org/bio-formats/).