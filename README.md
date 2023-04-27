# Quantum Tomography
This program allows to perform quantum tomography on quadrature marginal distributions to obtain the Wigner function. It also allows to calculate the density matrix in quadrature and the number representation from the generated Wigner function.  
These functionalities are presented in a user-friendly graphical interface and take advantage of multiprocessing techniques to provide fast results.

<img width="80%" alt="Screenshot 2022-10-30 at 7 52 04 PM" src="https://user-images.githubusercontent.com/44838230/198908531-4659e160-90d2-48f7-aea1-8a06d42f3ec7.png">


## Installation
### Install using binary files
You can use the installers for your system to easily install the software. 
- **[MacOS (Ventura 13.0)](https://github.com/amartinez1224/quantum-tomography/releases/download/v2.0/QuantumTomography.app.zip)**
- **Windows** (Cooming soon)

### Clone repo
You can also directly clone the repo and run the `gui.py` file (works in Windows and MacOS).
```
$ git clone https://github.com/amartinez1224/quantum-tomography.git
$ cd quantum-tomography
$ pip3 install -r requirements.txt
$ python3 gui.py
```

### Use a virtual environment (optional)
It's common practice to use a virtual environment to run Python applications. Virtual environments allow you to use an specific Python version (or the same one as the OS) and install specific versions of libraries that match those required by  the application without creating version conflicts with the global libraries in your OS. You can find more information on virtual environments [here](https://docs.python.org/3/library/venv.html).

To create a virtual environment using *venv* run
```
$ python3 -m venv DIR_OF_NEW_ENV
```
To run the `gui.py` file on the *venv* run
```
$ git clone https://github.com/amartinez1224/quantum-tomography.git
$ cd quantum-tomography
$ source DIR_OF_NEW_ENV/bin/activate
$ pip3 install -r requirements.txt
$ python3 gui.py
```

To exit the *venv* run
```
deactivate
```


## Usage
To perform a tomography, select an input file using the button in the lower right corner. Then click on the *Tomography* button. 

To tune the tomography, use the values on the right and then regenerate the tomography by clicking on the *Tomography* button:
- X min, X max: set the range of the X axis.
- Y min, Y max: set the range of the Y axis.
- Density: changes the amount of points generated for the Wigner function. The size of the generated Wigner function function is *density*\**density*.
- Kc: the value for the cutoff frequency.

You can also change the color of the graph and the angle from which it is viewed. For this you don't have to regenerate the tomography. In addition, the graph is fully interactive, so you can use the mouse to scroll.

At the bottom of the image you will find three options:
- Save data: saves the generated Wigner function in a JSON file.
- Density matrix quadrature: shows the density matrix in the quadrature space.
- Density matrix Fock: shows the density matrix in the Fock space. You can select the number of photons *n* and *m*.

### Input file
The input file has to be in the following JSON format:
```
{
  "phi": [],
  "x": [],
  "pr": [
    []
  ]
}
```
where *phi* is a one dimensional array that contains the angle of each marginal distribution, *x* another one dimensional array that contains the voltage range of the distributions and *pr* is a two dimensional array of size *phi*\**x* that contains the distribution.

## Samples
The folder *SampleStates* contains samples of simulated marginal distributions for the vacuum, coherent, squeezed and thermal states.

**sim_vacuum**
- Voltage range: -8V, 8V

**sim_thermal.json**:
- Voltage range: -8V, 8V

**sim_coherent.json**
- Voltage range: -8V, 8V
- Alpha: 3
- Angle: 180°

**sim_squeezed.json**
- Voltage range: -8V, 8V
- Alpha: 0
- Angle: 0°
- Eta: 0°

**sim_squeezed_vacuum.json**
- Voltage range: -8V, 8V
- Alpha: 3
- Angle: 90°
- Eta: 45°


## Simulate data
You can simulate your own data using the `simulateStates.py` script. The script has an inline interface for guidance.

## Read scope data
Data taken from a scope can be processed with the help of the `readScopeData.py` script.

## Wigner functions online gallery
Coming soon

## Documents
You can find the monograph that lead to this development [here](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/monograph.pdf).
### Lab guides
- [Lab guide 1](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/Lab_guide_1.pdf)
- [Lab guide 2](https://github.com/amartinez1224/quantum-tomography/blob/567bfe7e2c01903786542dab36ff7c925361e157/docs/Lab_guide_2.pdf)
