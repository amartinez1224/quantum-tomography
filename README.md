# Quantum Tomography
This program allows to perform quantum tomography on quadrature marginal distributions to obtain the Wigner function. It also allows to calculate the density matrix in quadrature and the number representation from the generated Wigner function.  
These functionalities are presented in a user-friendly graphical interface and take advantage of multiprocessing techniques to provide fast results.

<img width="80%" alt="Screenshot 2022-10-30 at 7 52 04 PM" src="https://user-images.githubusercontent.com/44838230/198908531-4659e160-90d2-48f7-aea1-8a06d42f3ec7.png">


## Installation
### Install using binary files
You can use the installers for your system to easily install the software. 
- **MacOS (Ventura 13.0):** 
- **Windows:** Cooming soon

### Clone repo
You can also directly clone the repo and run the *gui.py* file.
```
$ git clone https://github.com/amartinez1224/quantum-tomography.git
$ cd quantum-tomography
$ pip3 install -r requirements.txt
$ python3 gui.py
```

## Usage
To perform a tomography, select an input file using the button in the lower right corner. Then click on the *Tomography* button. 

To tune the tomography, use the values on the right and then regenerate the tomography by clicking on the *Tomography* button:
- Q min, Q max: set the range of the Q axis.
- P min, P max: set the range of the P axis.
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
