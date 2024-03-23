# VonMisesDesign

### Overview
The following code is an interactive airfoil generation tool. It uses conformal mapping of a circle to an airfoil via a Von Mises Transform. The underlying equations are based on a potential flow model, therefore all results assume inviscid, irrotational, incompressible flow. By shifting the center of the circle and the poles, unique airfoil shapes can be achieved. The airfoil shape, lift coefficient, and pressure/velocity distributions are updated in realtime.

![VMD-GUI](https://github.com/mbcapece3/VonMisesDesign/assets/104041016/d8dc2474-70d7-4059-9ebb-e12c34ede2ce)

### Usage Instructions
To generate an airfoil using VonMisesDesign perform the following steps:
- Run the VonMisesDesign.py file
- Drag the red dot in the center of the circle to the desired position.
  - With only two poles (represented by black dots), this represents a Joukowsky airfoil
  - For the airfoil to be physically valid, the X coordinate of the circle center must be negative
  - For the airfoil to be physically valid, all poles must be within the circle
- For additional control over the geometry, click the Add Pole button
  - Additional poles will appear at (0,0). They can be dragged to a different value as desired
  - Manipulating a single pole will cause other poles to shift to ensure the sum of all poles is (0,0)
- To remove a pole, click the Remove Pole button
  - Removing a pole will cause all remaining poles to reset their positions
- Adjust angle of attack using the slider at the bottom
- Toggle between pressure and velocity distribution display using the buttons on the right
- To save a file, click Create Save File. The terminal will prompt the user for a filename
  - Note: When saving a file, do not include an extension in the filename prompt
- To load a file, click Load Save File. The terminal will prompt the user for a filename
  - Note: When loading a file, the file extension must be included in the filename
- To export the airfoil as a set of coordinates, click Export DAT File. The terminal will prompt the user for a filename and an airfoil name
  - Note: When exporting a file, do not include an extension in the filename prompt

### How To Run
To run this program locally, python must be first be installed on your computer. It can be downloaded from here: [Python](https://www.python.org/downloads/)

To run this program locally, the following dependencies are required:
- [matplotlib](https://matplotlib.org/stable/)
  - To install, run "pip install matplotlib" in command prompt
- [numpy](https://numpy.org/install/)
  - To install, run "pip install numpy" in command prompt
- [pandas](https://pandas.pydata.org/docs/getting_started/index.html)
  - To install, run "pip install pandas" in command prompt
