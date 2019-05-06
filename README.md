# jointforces

![MCF7-raw](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/gifs/mcf7-raw.gif)

A Python package for conducting 3D traction force microscopy on multicellular aggregates (so-called *spheroids*). `jointforces` provides high-level interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO), facilitating material simulations of contracting multicellular aggregates in highly non-linear biopolymer gels such as collagen. Additionally, `jointforces` provides an easy-to-use API for analyzing time-lapse images of contracting multicellular aggregates using the particle image velocimetry framework `OpenPIV`.

## Installation
TODO

## Documentation
This section details a complete walk-through of the analysis of a 3D traction force microscopy experiment. The data we use in this example can be downloaded here (TODO). The data consists of XXX consecutive images of a MCF7 tumor spheroid contracting for 12h within a 1.2mg/ml collagen gel. The contractility of a multicellular aggregate is estimated by comparing the measured deformation field induced by the contracting cells to a set of simulated deformation fields of varying contractility. While we provide precompiled simulation for various materials, this documentation not only covers the force reconstruction step, but also the material simulations.

The module is imported in Python as:

```python
import jointforces as jf
```

### 1. Setting up interfaces
The first step is to tell `jointforces` where `Gmsh` and `SAENO` are stored. This only has to be done once, or again whenever one of the programs os moved/re-installed.

```python
jf.set_gmsh_path(r'C:\Software\gmsh-4.3.0-Windows64-sdk')
jf.set_saeno_path(r'C:\Software\SAENO')
```

### 2. Mesh generation
Here, we create a spherical bulk of material (with a radius `r_outer=2cm`, emulating the biopolymer network) with a small, centered spherical inclusion (with a radius `r_inner=100µm`, emulating the multicellular aggregate). The keyword-argument `length_factor` determines the mesh granularity:

```python
jf.mesh.spherical_inclusion('spherical-inclusion.msh', 
                            r_inner=100, 
                            r_outer=20000, 
                            length_factor=0.2)
```

The resulting mesh is saved as `spherical-inclusion.msh` and can be displayed in `Gmsh` using the following command:

```python
jf.mesh.show_mesh('spherical-inclusion.msh')
```

![Gmsh](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/images/gmsh.png)

### 3. Material simulations
Having created the mesh geometry, we define appropriate boundary conditions and material parameters to emulate the contracting multicellular aggregate *in silico*. Here, we assume a constant in-bound pressure on the surface of the spherical inclusion (emulating the cells pulling on the surrounding matrix), and no material displacements on the outer boundary of the bulk material (emulating a hard boundary such a the walls of a Petri dish). The goal of the simulation is to obtain the displacements in the surrounding material that are the effect of the pressure exerted by the multicellular aggregate. The following command executes a simulation assuming a pressure of 100Pa and 1.2mg/ml collagen matrix as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685). After successful optimization, the results are stored in the output folder `simu`:

```python
jf.simulation.spherical_contraction('test.msh', 'simu', 100, jf.materials.collagen12)
```

More detailed information about the output files of a simulation can be found the [Wiki of the SAENO project](https://github.com/Tschaul/SAENO/wiki). The file `parameters.txt` contains all parameters used int he simulation. Note that `jointforces` provides functions that read in the resulting files of material simulations again to facilitate a comparison of simulated and measured deformation fields, as detailed below.

#### Material parameters
`jointforces` provides "pre-configured" material types for collagen gels of three different concentrations (0.6, 1.2, and 2.4mg/ml). Detailed protocols for reproducing these gels can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685) and [Condor et al. (2017)](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpcb.24). Furthermore, one can define a linear elastic material with a specified `stiffness` (in Pa) with:

```
jf.materials.linear(stiffness)
```

To define a non-linear elastic material, use the `custom` material type:

```
jf.materials.custom(K_0, D_0, L_S, D_S)
```

Non-linear materials are characterized by four parameters:
- `K_0`: the linear stiffness (in Pa)
- `D_0`: the rate of stiffness variation during fiber buckling
- `L_S`: the onset strain for strain stiffening
- `D_S`: the rate of stiffness variation during strain stiffening
A full description of the non-linear material model and the parameters can be found in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685)

#### Running simulations in parallel
To be able to estimate the contractility of a multicellular aggregate by relating the measured deformation field to simulated ones, we need to execute a set of simulations that cover a wide range of contractile pressures. To speed up this process, `jointforces` provides a parallelization method that distributes individual instances of the `SAENO` optimizer across different CPU cores. The following code snippet runs 100 simulations ranging from pressure values of 0.1Pa to 1000Pa (logarithmically spaced):

```
jf.simulation.distribute('jf.simulation.spherical_contraction',
                         const_args={'meshfile': 'spherical-contraction.msh', 'outfolder': 'simu',
                         'material': jf.materials.collagen12},
                         var_arg='pressure', start=0.1, end=1000, n=100, log_scaling=True)
```

The method automatically creates subfolders within the output-folder `simu`, called `simulation000000`, `simulation000001`, and so on, plus a file `pressure-values.txt` that contains the list of pressure values used in the simulations.

### 4. Pressure lookup tables

### 5. Particle image velocimetry

### 6. Force reconstruction

## Dependencies
*jointforces* is tested on Python 3.7 and a Windows 10 64bit system. It depends on ... All except ... are already included in the [Anaconda distribution](https://www.continuum.io/downloads) of Python. Windows users may also take advantage of pre-compiled binaries for all dependencies, which can be found at [Christoph Gohlke's page](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

## License
[The MIT License (MIT)](https://github.com/christophmark/jointforces/blob/master/LICENSE)
