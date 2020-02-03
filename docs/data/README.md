# Exemplary data

This folder contains relevant files to test/apply the analysis workflow of `jointforces`.

## [MCF7-time-lapse.zip](https://github.com/christophmark/jointforces/blob/master/docs/data/MCF7-time-lapse.zip)

Exemplary time-lapse images of a MCF7 spheroid containing 4000 cells (at the time of seeding) that is embedded in a 1.2mg/ml collagen gel together with silica beads as fiducial markers.

The image series consists of 145 images that are recorded with a 5 minute interval and cover a total time interval of 12 hours. The pixel size is 1.29 µm/pixel. Each image is a minimum projection of a brightfield image stack that covers 100µm across the equatorial plane of the spheroid.

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

