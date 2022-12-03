# Exemplary data

This folder contains relevant files (or the respective links) to test/apply the analysis workflow of `jointforces`.

## [MCF7-time-lapse.zip](https://www.dropbox.com/s/b6uztm3tgdo491p/MCF7-time-lapse.zip?dl=1)

Exemplary time-lapse images of a MCF7 spheroid containing 4000 cells (at the time of seeding) that is embedded in a 1.2mg/ml collagen gel together with silica beads as fiducial markers.

The image series consists of 145 images that are recorded with a 5 minute interval and cover a total time interval of 12 hours. The pixel size is 1.29 µm/pixel. Each image is a minimum projection of a brightfield image stack that covers 100µm across the equatorial plane of the spheroid.




## [k12_BatchA.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/k12_BatchA.pkl)

Lookup table for 1.2 mg/ml collagen measured for a new batch of collagen as described in (unpublished yet). These lookup-tables can be used between different python versions. The file was created by the following code:


```python
    import jointforces as jf
    
    out_folder = 'k12_mbioscience_2020_v2'
    out_table = 'k12_mbioscience_2020_table_v2.pkl'
    
    K_0 = 1449  
    D_0 = 0.00215
    L_S = 0.032
    D_S = 0.055 
      
    jf.simulation.distribute_solver('jf.simulation.spherical_contraction_solver',
    
                              const_args={'meshfile': meshfile_loc,     # path to the provided or the new generated mesh
                                         'outfolder': out_folder,    # output folder to store individual simulations
                                          'max_iter': 600,
                                          'step': 0.0033,  # before 0.033
                                          'material': jf.materials.custom(K_0, D_0, L_S, D_S) },      # Enter your own material parameters here
                                          var_arg='pressure', start=0.1, end=10000, n=150, log_scaling=True, n_cores=2, get_initial=True)
      
    lookup_table = jf.simulation.create_lookup_table_solver(out_folder, x0=1, x1=50, n=100)    # output folder for combining the individual simulations
    jf.simulation.save_lookup_table(lookup_table,out_table)
```

## [k06_BatchA.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/k06_BatchA.pkl)

Lookup table for 0.6 mg/ml collagen measured for a new batch of collagen as described in (unpublished yet). These lookup-tables can be used between different python versions.

## [k24_BatchA.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/k24_BatchA.pkl)

Lookup table for 2.4 mg/ml collagen measured for a new batch of collagen as described in (unpublished yet). These lookup-tables can be used between different python versions.




--------------------------------------------------------------------------------------------------------------------------------------------------
**The following lookup-tables might not be compatible with newer python versions anymore (use *save_lookup_table* instead of *save_lookup_functions* in order to use lookup-tables between different python versions).**




## [collagen12.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/collagen12.pkl)

Lookup table for 1.2 mg/ml collagen as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685). The file was created by the following code:

```python
import jointforces as jf

jf.mesh.spherical_inclusion('spherical-inclusion.msh',
                            r_inner=100,
                            r_outer=10000,
                            length_factor=0.05)

jf.simulation.distribute('jf.simulation.spherical_contraction',
                         const_args={'meshfile': 'spherical-inclusion.msh',
                                     'outfolder': 'D:/material-sims-collagen12',
                                     'material': jf.materials.collagen12},
                         var_arg='pressure', start=0.1, end=10000, n=150, log_scaling=True, n_cores=3)

lookup_table = jf.simulation.create_lookup_table('D:\material-sims-collagen12', x0=1, x1=50, n=150)

get_displacement, get_pressure = jf.simulation.create_lookup_functions(lookup_table)

jf.simulation.save_lookup_functions(get_displacement, get_pressure, 'collagen12.pkl')
```


## [collagen06.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/collagen06.pkl)

Lookup table for 0.6 mg/ml collagen as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685).

## [collagen24.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/collagen24.pkl)

Lookup table for 2.4 mg/ml collagen as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685).

## [fibrin40.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/fibrin40.pkl)

Lookup table for 4.0 mg/ml fibrin as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685).

## [matri10.pkl](https://github.com/christophmark/jointforces/blob/master/docs/data/matri10.pkl)

Lookup table for 10.0 mg/ml Matrigel as described in [Steinwachs et al. (2016)](https://www.nature.com/articles/nmeth.3685).

## Further Lookup Tables

For individual nonlinear materials, the material properties can be determined by using [saenopy](https://saenopy.readthedocs.io/en/latest/fit_material_parameters.html) and can then be used to create a new lookup table. Lookup tables for arbitrary linear elastic material of different stiffness can be easily created using an interpolation function as follows:

```python
jf.simulation.linear_lookup_interpolator(emodulus=250, output_newtable="linear-lookup-emodul-250Pa.pkl", 
```


