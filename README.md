# jointforces

A Python package for conducting 3D traction force microscopy on multicellular aggregates (so-called *spheroids*). `jointforces` provides high-level interfaces to the open-source finite element mesh generator [`Gmsh`](http://gmsh.info/) and to the [Python port](https://github.com/rgerum/saenopy) of the network optimizer [`SAENO`](https://github.com/Tschaul/SAENO), facilitating material simulations of contracting multicellular aggregates in highly non-linear biopolymer gels such as collagen. Additionally, `jointforces` provides an easy-to-use API for analyzing time-lapse images of contracting multicellular aggregates using the particle image velocimetry framework `OpenPIV`.

![Loading GIF...](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/gifs/mcf7-raw.gif)

## Installation
The current version of this package can be downloaded as a zip file [here](https://github.com/christophmark/jointforces/zipball/master), or by cloning this repository. After unzipping, run the following command within the unzipped folder: `pip install -e .`. This will automatically download and install all other required packages.

### Troubleshooting
In some cases the installation of the required package `OpenPIV` fails due to missing compilers. In that case, Windows users may use the pre-compiled binaries provided [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/#openpiv). Download the approriate `*.whl`-File for your Python version and install via `pip install *.whl`. This works not only for `OpenPIV`, but for all packages hosted on [Christoph Gohlke's page](https://www.lfd.uci.edu/~gohlke/pythonlibs/). 

## Minimal example
`jointforces` provides [example data](https://github.com/christophmark/jointforces/tree/master/docs/data) and [pre-computed material simulations](https://github.com/christophmark/jointforces/tree/master/docs/data) for 1.2mg/ml collagen gels as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685). The following code snippet...
- loads a series of images
- segments the spheroid in each image
- applies particle image velocimetry to each pair of subsequent images
- compares the resulting matrix deformations to a range of material simulations
- computes the contractile pressure and overall contractility of the spheroid over time

```python
import jointforces as jf

jf.piv.compute_displacement_series('MCF7-time-lapse',  # image folder
                                   '*.tif',            # file pattern
                                   'MCF7-piv',         # output folder
                                   window_size=40,     # PIV window
                                   cutoff=650)         # PIV cutoff

jf.force.reconstruct('MCF7-piv',        # PIV output folder
                     'collagen12.pkl',  # lookup table
                     1.29,              # pixel size (µm)
                     'results.xlsx')    # output file

```

## Documentation
This section details a complete walk-through of the analysis of a 3D traction force microscopy experiment. The data we use in this example can be downloaded [here](https://github.com/christophmark/jointforces/tree/master/docs/data). The data consists of 145 consecutive images of a MCF7 tumor spheroid contracting for 12h within a 1.2mg/ml collagen gel. The contractility of a multicellular aggregate is estimated by comparing the measured deformation field induced by the contracting cells to a set of simulated deformation fields of varying contractility. While we provide precompiled simulations for various materials, this documentation not only covers the force reconstruction step, but also the material simulations.

The module is imported in Python as:

```python
import jointforces as jf
```

### 1. Mesh generation
Here, we create a spherical bulk of material (with a radius `r_outer=1cm`, emulating the biopolymer network) with a small, centered spherical inclusion (with a radius `r_inner=100µm`, emulating the multicellular aggregate). The keyword-argument `length_factor` determines the mesh granularity:

```python
jf.mesh.spherical_inclusion('spherical-inclusion.msh', 
                            r_inner=100, 
                            r_outer=10000, 
                            length_factor=0.05)
```

The resulting mesh is saved as `spherical-inclusion.msh` and can be displayed in `Gmsh` using the following command:

```python
jf.mesh.show_mesh('spherical-inclusion.msh')
```

![Gmsh](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/images/gmsh.png)

### 2. Material simulations
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
To be able to estimate the contractility of a multicellular aggregate by relating the measured deformation field to simulated ones, we need to execute a set of simulations that cover a wide range of contractile pressures. To speed up this process, `jointforces` provides a parallelization method that distributes individual instances of the `SAENO` optimizer across different CPU cores. The following code snippet runs 150 simulations ranging from pressure values of 0.1Pa to 10000Pa (logarithmically spaced):

```python
jf.simulation.distribute('jf.simulation.spherical_contraction',
                         const_args={'meshfile': 'spherical-contraction.msh', 'outfolder': 'simu',
                         'material': jf.materials.collagen12},
                         var_arg='pressure', start=0.1, end=10000, n=150, log_scaling=True)
```

The method automatically creates subfolders within the output-folder `simu`, called `simulation000000`, `simulation000001`, and so on, plus a file `pressure-values.txt` that contains the list of pressure values used in the simulations.

### 3. Pressure lookup tables
To compare a measured deformation field to a set of simulated ones, we need to create a lookup table that output the expected pressure for a given strain at a given distance from the spheroid (or the expected strain for a given pressure at a given distance). First, we convert the 3D displacement fields into a set of radial displacement curves:

```python
lookup_table = jf.simulation.create_lookup_table('simu')
```

Now we can interpolate between individual simulations to create lookup functions for pressure values and strain values:

```python
get_displacement, get_pressure = jf.simulation.create_lookup_functions(lookup_table)
```

Finally, we may save these functions to file, e.g. to load them again in another script:

```python
jf.simulation.save_lookup_functions(get_displacement, get_pressure, 'collagen12.pkl')
get_displacement, get_pressure = jf.simulation.load_lookup_functions('collagen12.pkl')
```

For example, we may now ask for the estimated pressure of a spheroid with radius `r` that induces a deformation of `0.2*r` at a distance of `2*r`:

```python
print(get_pressure(2, 0.2))
>>> 99.81708381387853
```

We provide a pre-computed lookup table for the standard 1.2mg/ml collagen gel [here](https://github.com/christophmark/jointforces/tree/master/docs/data). This lookup table has been created using the exact commands described above.

### 4. Particle image velocimetry

Up to this point, we have only covered material simulations, but not the analysis of measured time-lapse image series. To detect deformations in the material surrounding the spheroid, `jointforces` uses the [Particle Image Velocimetry](https://en.wikipedia.org/wiki/Particle_image_velocimetry) algorithm of the [`OpenPIV`](http://www.openpiv.net/openpiv-python/) package. The following command automatically reads in all image files that match a filterstring within a given folder, and computes the deformation fields between subsequent images, and saves overlay plots of the deformation fields. The exemplary data used in this example can be downloaded [here](https://github.com/christophmark/jointforces/tree/master/docs/data).

```python
jf.piv.compute_displacement_series('MCF7-time-lapse', '*.tif', 'MCF7-piv', 
                                   window_size=40, cutoff=650)
```

The command performs PIV on all `*.tif` files in the the folder `MCF7-time-lapse`, the results are saved int he folder `MCF7-piv`. The window size of the PIV algorithm should be chosen as small as possible to increase spatial resolution, but also large enough to contain multiple fiducial markers for an accurate detection of local material deformations. The `cutoff` parameter can be used to disregard all displacements that are detected further away from the center than the set value (e.g. because an optical coupler is visible in the corners of the image).

![Loading GIF...](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/gifs/mcf7-piv.gif)

### 5. Force reconstruction

Finally, we may use the lookup functions we have created above and use them to assign the best-fit pressure to each time step of the image series. Additionally, the user supplies the size of one pixel in the image in micrometers. With this information, the surface are of the spheroid is calculated to obtain the total contractility. The output is a [`Pandas`](https://pandas.pydata.org/) Dataframe containing mean values, median values and standard deviation of both pressure and contractility. If a filename is provided, the results are also saved as an Excel file:

```python
res = jf.force.reconstruct('MCF7-piv', 'collagen12.pkl', 6.45/5, 'MCF7-recon.xlsx')

t = np.arange(len(res))*5/60
mu = res['Mean Contractility (µN)']
std = res['St.dev. Contractility (µN)']

plt.figure(figsize=(6, 3))
plt.plot(t, mu, lw=2, c='C0')
plt.fill_between(t, mu-std, mu+std, facecolor='C0', lw=0, alpha=0.5)
plt.grid()
plt.xlabel('Time (h)')
plt.ylabel('Contractility (µN)')
plt.tight_layout()
plt.show()
```

![Plot](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/images/mcf7-recon.png)

![Table](https://raw.githubusercontent.com/christophmark/jointforces/master/docs/images/mcf7-excel.png)

## Dependencies
*jointforces* is tested on Python 3.7 and a Windows 10 64bit system. It depends on `numpy`, `pandas`, `matplotlib`, `scipy`, `scikit-image`, `openpiv`, `tqdm`, `natsort`, and `dill`. We recommend using the [Anaconda distribution](https://www.continuum.io/downloads) of Python. Windows users may also take advantage of pre-compiled binaries for all dependencies, which can be found at [Christoph Gohlke's page](http://www.lfd.uci.edu/~gohlke/pythonlibs/).

## License
The `jointforces` package itself is licensed under the [MIT License](https://github.com/christophmark/jointforces/blob/master/LICENSE). The [exemplary data](https://github.com/christophmark/jointforces/tree/master/docs/data) and the [precompiled binaries](https://github.com/christophmark/jointforces/tree/master/docs/bin) of SAENO are provided "as is", without any warranty of any kind, express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose and noninfringement.
