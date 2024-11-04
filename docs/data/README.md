# Exemplary data

This folder contains relevant files (or the respective links) to test/apply the analysis workflow of `jointforces`.

## [MCF7-time-lapse.zip](https://www.dropbox.com/s/b6uztm3tgdo491p/MCF7-time-lapse.zip?dl=1)

Exemplary time-lapse images of a MCF7 spheroid containing 4000 cells (at the time of seeding) that is embedded in a 1.2mg/ml collagen gel together with silica beads as fiducial markers.

The image series consists of 145 images that are recorded with a 5 minute interval and cover a total time interval of 12 hours. The pixel size is 1.29 µm/pixel. Each image is a minimum projection of a brightfield image stack that covers 100µm across the equatorial plane of the spheroid.




## [k12_BatchA.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/k12_BatchA.npy)

Lookup table for 1.2 mg/ml collagen of a collagen batch (A) as described in [Dynamic traction force measurements of migrating immune cells in 3D matrices](https://doi.org/10.1101/2022.11.16.516758). The file was created by the following code:


```python
    import jointforces as jf
    
    out_folder = 'Collagen12_BatchA'
    out_table = 'Collagen12_BatchA.npy'
    
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

## [Collagen06_BatchA.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen06_BatchA.npy)

Lookup table for 0.6 mg/ml collagen of a collagen batch (A) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 


## [Collagen06_BatchA.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen12_BatchA.npy)

Lookup table for 1.2 mg/ml collagen of a collagen batch (A) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 

## [Collagen24_BatchA.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen24_BatchA.npy)

Lookup table for 2.4 mg/ml collagen of a collagen batch (A) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 



## [Collagen06_BatchC.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen06_BatchC.npy)

Lookup table for 0.6 mg/ml collagen of a collagen batch (C) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 

## [Collagen12_BatchC.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen12_BatchC.npy)

Lookup table for 1.2 mg/ml collagen of a collagen batch (C) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 
## [Collagen24_BatchC.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen24_BatchC.npy)

Lookup table for 2.4 mg/ml collagen of a collagen batch (C) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 


## [Collagen_1mgml_BatchD.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen_1mgml_BatchD.npy)
Lookup table for 1.0 mg/ml collagen of a collagen batch (D) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 

## [Collagen_2mgml_BatchD.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen_2mgml_BatchD.npy)
Lookup table for 2.0 mg/ml collagen of a collagen batch (D) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 

## [Collagen_3mgml_BatchD.npy](https://github.com/christophmark/jointforces/blob/master/docs/data/Collagen_3mgml_BatchD.npy)
Lookup table for 3.0 mg/ml collagen of a collagen batch (D) as described in [Dynamic traction force measurements of migrating immune cells in 3D biopolymer matrices](https://doi.org/10.1101/2022.11.16.516758](https://doi.org/10.1038/s41567-024-02632-8 ). 





## Lookup Tables for linear elastic fiber materials

Lookup tables for arbitrary linear elastic material of different stiffness can be easily created using an interpolation function based on a pre-calculated lookuptable as follows:

```python
jf.simulation.linear_lookup_interpolator(emodulus=250, output_newtable="linear-lookup-emodul-250Pa.pkl", 
```


## Further Lookup Tables

For individual nonlinear materials, the material properties can be determined by using [saenopy](https://saenopy.readthedocs.io/en/latest/fit_material_parameters.html) and can then be used to create a new lookup table. 












