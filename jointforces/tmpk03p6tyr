import os
import dill
import shutil
import numpy as np
from glob import glob
from subprocess import call
from subprocess import Popen
from subprocess import PIPE
from time import sleep
from tqdm import tqdm
from scipy.interpolate import LinearNDInterpolator
from saenopy import Solver
import saenopy
from saenopy.materials import SemiAffineFiberMaterial
import pandas as pd


def read_meshfile(meshfile, r_inner=None, r_outer=None):
    # open mesh file
    with open(meshfile, 'r') as f:
        lines = f.readlines()

    # if r_inner or r_outer are defined in the mesh-file, ignore user input
    meshinfo = {}
    try:
        paragraph = lines[lines.index('$Jointforces\n') + 1:lines.index('$EndJointforces\n')]
        for line in paragraph:
            key, value = line.strip().split('=')
            try:
                meshinfo[key] = float(value)
            except ValueError:
                meshinfo[key] = value

        r_inner = meshinfo['r_inner']
        r_outer = meshinfo['r_outer']

        print('Geometry for spherical contraction is defined in mesh file (r_inner={:.2f}, r_outer={:.2f}).'.format(r_inner, r_outer))
        if (r_inner is not None) or (r_outer is not None):
            print('Will ignore user-defined values of r_inner and r_outer.')
    except ValueError:
        if r_inner is None:
            raise ValueError('r_inner not defined')
        if r_outer is None:
            raise ValueError('r_outer not defined')
            
    # scale radii to meter
    r_inner *= 10**-6
    r_outer *= 10**-6
    
    # transform nodes and connection in SAENO format
    # nodes
    index_nodes = lines.index('$Nodes\n')
    n_nodes = int(lines[index_nodes + 1])

    coords = np.zeros((n_nodes, 3))
    for i in range(n_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i + index_nodes + 2].split()[1:]])
        
      
    # connections
    index_elements = lines.index('$Elements\n')
    n_elements = int(lines[index_elements + 1])

    tets = np.zeros((n_elements, 4))
    for i in range(n_elements):
        tets[i] = lines[i + index_elements + 2].split()[-4:]

    # to start with 0 and not 1
    tets -= 1
    
    return coords, tets, r_inner, r_outer
    

def spherical_contraction_solver(meshfile, outfolder, pressure, material, r_inner=None, r_outer=None, logfile=False, initial_displacenemts=None, max_iter = 600, step = 0.0033, conv_crit = 0.01):
    
    
    coords, tets, r_inner, r_outer = read_meshfile(meshfile, r_inner, r_outer)

    # read in material parameters
    K_0 = material['K_0']
    D_0 = material['D_0']
    L_S = material['L_S']
    D_S = material['D_S']

    # create output folder if it does not exist, print warning otherwise
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    else:
        print('WARNING: Output folder already exists! ({})'.format(outfolder))


    # Initialize saenopy solver opbject
    M = Solver()
    material_saenopy = SemiAffineFiberMaterial(K_0, D_0, L_S, D_S)
    M.setMaterialModel(material_saenopy)
    M.setNodes(coords)
    M.setTetrahedra(tets)
     
    # define boundary conditions
    distance = np.sqrt(np.sum(coords ** 2., axis=1))
    mask_inner = distance < r_inner * 1.001
    mask_outer = distance > r_outer * 0.999

    # Save Node Density at inner and outer sphere
    # Area per inner node
    A_node_inner = (np.pi*4*(r_inner)**2)/np.sum(mask_inner) 
    # simple sqrt as spacing
    inner_spacing = np.sqrt(A_node_inner)    
    
    # Area per outer node
    A_node_outer = (np.pi*4*(r_outer)**2)/np.sum(mask_outer)   
    # simple sqrt as spacing
    outer_spacing = np.sqrt(A_node_outer)
    
    print ('Inner node spacing: '+str(inner_spacing*1e6)+'µm')
    print ('Outer node spacing: '+str(outer_spacing*1e6)+'µm')

    
    # displacements are everywhere NaN
    bcond_displacement = np.zeros((len(coords), 3))*np.nan
    # except of a the outer border
    bcond_displacement[mask_outer] = 0
    
    # forces are everywhere 0
    bcond_forces = np.zeros((len(coords), 3))
    # except at the outer boder there they are NaN
    bcond_forces[mask_outer] = np.nan
    # and at the innter border they depend on the pressure
    bcond_forces[mask_inner] = coords[mask_inner]
    bcond_forces[mask_inner] /= distance[mask_inner, None]
    A_inner = 4 * np.pi * r_inner ** 2.
    force_per_node = pressure * A_inner / np.sum(mask_inner)
    bcond_forces[mask_inner, :3] *= force_per_node
    
    # give the boundary conditions to the solver
    M.setBoundaryCondition(bcond_displacement, bcond_forces)

    if initial_displacenemts is not None:
        M.setInitialDisplacements(initial_displacenemts)

    # create info file with all relevant parameters of the simulation
    parameters = r"""K_0 = {}
D_0 = {}
L_S = {}
D_S = {}
PRESSURE = {}
FORCE_PER_SURFACE_NODE = {}
INNER_RADIUS = {} µm
OUTER_RADIUS = {} µm
INNER_NODE_SPACING = {} µm
OUTER_NODE_SPACING = {} µm
SURFACE_NODES = {}
TOTAL_NODES = {}""".format(K_0, D_0, L_S, D_S, pressure, force_per_node, r_inner*1e6, r_outer*1e6, inner_spacing*1e6,
                           outer_spacing*1e6, np.sum(mask_inner), len(coords))

    with open(outfolder + "/parameters.txt", "w") as f:
        f.write(parameters)
        
    # solve the boundary problem
    M.solve_boundarycondition(stepper=step, i_max=max_iter, rel_conv_crit=conv_crit, relrecname=outfolder + "/relrec.txt") #, verbose=True
    M.save(outfolder + "/solver.npz")
    
  
def distribute_solver(func, const_args, var_arg='pressure', start=0.1, end=1000, n=120, log_scaling=True, n_cores=None, get_initial=True, max_iter = 600, step = 0.0033, conv_crit = 0.01, callback=None):
    # get_intial = True takes the deformationfield from previous simulation as start values for the next simulations, which reduces computation time
    
    # by default use spherical contraction as function
    func = spherical_contraction_solver
    
    if n_cores is None:
        n_cores = os.cpu_count()

    if log_scaling:
        values = np.logspace(np.log10(start), np.log10(end), n, endpoint=True)
    else:
        values = np.linspace(np.log10(start), np.log10(end), n, endpoint=True)

    outfolder = const_args['outfolder']
    del const_args['outfolder']

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    np.savetxt(outfolder+'/'+var_arg+'-values.txt', values)

    values = list(values)

    if n_cores == 1:
        U = None
        for index, v in enumerate(values):
            func(const_args["meshfile"], outfolder + '/simulation' + str(index).zfill(6),
                    v, const_args["material"], None, None, False, U, max_iter, step, conv_crit)
            # get the last displacement
            name = outfolder + '/simulation' + str(index).zfill(6) + "/solver.npz"
            if os.path.exists(name):
                M = saenopy.load(name)
                U = M.U
            # optionally call the callback
            if callback is not None:
                callback(index, len(values))
    else:
        index = 0
        processes = []
        import multiprocessing
        while True:
            processes = [p for p in processes if p.is_alive()]

            if len(processes) < n_cores and index < len(values):
                v = values[index]
                U = None
                if get_initial==True:
                    for i in range(index-3, -1, -1):
                        name = outfolder+'/simulation'+str(i).zfill(6) + "/solver.npz"
                        if os.path.exists(name):
                            M = saenopy.load(name)
                            U = M.U
                            break
                p = multiprocessing.Process(target=func,
                                            args=(const_args["meshfile"],  outfolder+'/simulation'+str(index).zfill(6),
                                                  v, const_args["material"], None, None, False, U, max_iter, step, conv_crit))

                p.start()
                processes.append(p)
                if callback is not None:
                    callback(index, len(values))
                index += 1
            sleep(1.)

            if len(processes) == 0:
                break
    if callback is not None:
        callback(index, len(values))
    return
    
  


def extract_deformation_curve_solver(folder, x):
    # get simulation parameters
    with open(folder+'/parameters.txt', 'r') as f:
        lines = f.readlines()

    parameters = {}
    for line in lines:
        try:
            key, value = line.split('= ')
            value = value.split(' ')[0]
            parameters[key.strip()] = float(value.strip())
        except:
            pass

    # load coordinates
    M = saenopy.load(folder + "/solver.npz")
    coords = M.R 
    # load displacements
    displ = M.U
    
    # compute binned normed displacements and normed coordinates
    u = np.sqrt(np.sum(coords ** 2., axis=1)) / (parameters['INNER_RADIUS']*10**-6)
    v = np.sqrt(np.sum(displ ** 2., axis=1)) / (parameters['INNER_RADIUS']*10**-6)

    y = np.array([np.nanmedian(v[(u >= x[i]) & (u < x[i + 1])]) for i in range(len(x) - 1)])

    # save results
    results = {'pressure': parameters['PRESSURE'], 'bins': x, 'displacements': y}
    return results
  
    
def create_lookup_table_solver(folder, x0=1, x1=50, n=100):
    subfolders = glob(folder+'/*/')

    x = np.logspace(np.log10(x0), np.log10(x1), n+1, endpoint=True)
    x_center = 10**(0.5*(np.log10(x[1:]) + np.log10(x[:-1])))

    pressure_values = []
    displacement_curves = []

    for subfolder in tqdm(subfolders):
        res = extract_deformation_curve_solver(subfolder, x)
        pressure_values.append(res['pressure'])
        displacement_curves.append(res['displacements'])

    pressure_values = np.array(pressure_values)
    displacement_curves = np.array(displacement_curves)

    return {'pressure': pressure_values, 'x': x_center, 'y': displacement_curves}
  

  
    
"""
Command_line_version 
"""   
    
    
    
    
def spherical_contraction(meshfile, outfolder, pressure, material, r_inner=None, r_outer=None, logfile=False,  max_iter = 600, step = 0.0033, conv_crit = 0.01):
    # open mesh file
    with open(meshfile, 'r') as f:
        lines = f.readlines()

    # if r_inner or r_outer are defined in the mesh-file, ignore user input
    meshinfo = {}
    try:
        paragraph = lines[lines.index('$Jointforces\n') + 1:lines.index('$EndJointforces\n')]
        for line in paragraph:
            key, value = line.strip().split('=')
            try:
                meshinfo[key] = float(value)
            except ValueError:
                meshinfo[key] = value

        r_inner = meshinfo['r_inner']
        r_outer = meshinfo['r_outer']

        print('Geometry for spherical contraction is defined in mesh file (r_inner={:.2f}, r_outer={:.2f}).'.format(r_inner, r_outer))
        if (r_inner is not None) or (r_outer is not None):
            print('Will ignore user-defined values of r_inner and r_outer.')
    except ValueError:
        if r_inner is None:
            raise ValueError('r_inner not defined')
        if r_outer is None:
            raise ValueError('r_outer not defined')
    
    # scale radii to meter
    r_inner *= 10**-6
    r_outer *= 10**-6

    # read in material parameters
    K_0 = material['K_0']
    D_0 = material['D_0']
    L_S = material['L_S']
    D_S = material['D_S']

    # create output folder if it does not exist, print warning otherwise
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
    else:
        print('WARNING: Output folder already exists! ({})'.format(outfolder))

    # transform nodes and connection in SAENO format
    # nodes
    index_nodes = lines.index('$Nodes\n')
    n_nodes = int(lines[index_nodes + 1])

    coords = np.zeros((n_nodes, 3))
    for i in range(n_nodes):
        coords[i] = np.array([np.float(x) for x in lines[i + index_nodes + 2].split()[1:]])
    np.savetxt(outfolder + '/coords.dat', coords)

    # connections
    index_elements = lines.index('$Elements\n')
    n_elements = int(lines[index_elements + 1])

    tets = np.zeros((n_elements, 4))
    for i in range(n_elements):
        tets[i] = lines[i + index_elements + 2].split()[-4:]
    np.savetxt(outfolder + '/tets.dat', tets, fmt='%i')

    # define boundary conditions
    distance = np.sqrt(np.sum(coords ** 2., axis=1))
    mask_inner = distance < r_inner * 1.001
    mask_outer = distance > r_outer * 0.999

    # Save Node Density at inner and outer sphere
    # Area per inner node
    A_node_inner = (np.pi*4*(r_inner)**2)/np.sum(mask_inner) 
    # simple sqrt as spacing
    inner_spacing = np.sqrt(A_node_inner)    
    
    # Area per outer node
    A_node_outer = (np.pi*4*(r_outer)**2)/np.sum(mask_outer)   
    # simple sqrt as spacing
    outer_spacing = np.sqrt(A_node_outer)
    
    print ('Inner node spacing: '+str(inner_spacing*1e6)+'µm')
    print ('Outer node spacing: '+str(outer_spacing*1e6)+'µm')

    bcond = np.zeros((len(coords), 4))
    bcond[:, 3] = 1.

    # fixed displacements for outer boundary
    bcond[mask_outer, 3] = 0

    # fixed non-zero force at spheroid surface
    bcond[mask_inner, :3] = coords[mask_inner, :3]
    bcond[mask_inner, :3] /= distance[mask_inner, None]

    A_inner = 4 * np.pi * r_inner ** 2.
    force_per_node = pressure * A_inner / np.sum(mask_inner)
    bcond[mask_inner, :3] *= force_per_node

    np.savetxt(outfolder + '/bcond.dat', bcond)

    # define initial configuration
    iconf = np.zeros((len(coords), 3))
    np.savetxt(outfolder + '/iconf.dat', iconf)

    # create config file for SAENO    
    config = r"""MODE = relaxation
BOXMESH = 0
FIBERPATTERNMATCHING = 0
REL_CONV_CRIT = {}
REL_ITERATIONS = {}
REL_SOLVER_STEP = {}
K_0 = {}
D_0 = {}
L_S = {}
D_S = {}
CONFIG = {}\config.txt
DATAOUT = {}""".format(conv_crit, max_iter, step, K_0, D_0, L_S, D_S, os.path.abspath(outfolder), os.path.abspath(outfolder))

    with open(outfolder + "/config.txt", "w") as f:
        f.write(config)

    # create info file with all relevant parameters of the simulation
    parameters = r"""K_0 = {}
D_0 = {}
L_S = {}
D_S = {}
PRESSURE = {}
FORCE_PER_SURFACE_NODE = {}
INNER_RADIUS = {} µm
OUTER_RADIUS = {} µm
INNER_NODE_SPACING = {} µm
OUTER_NODE_SPACING = {} µm
SURFACE_NODES = {}
TOTAL_NODES = {}""".format(K_0, D_0, L_S, D_S, pressure, force_per_node, r_inner*1e6, r_outer*1e6, inner_spacing*1e6,
                           outer_spacing*1e6, np.sum(mask_inner), len(coords))

    with open(outfolder + "/parameters.txt", "w") as f:
        f.write(parameters)

    # Create log file if activated
    if logfile:
        
        # create log file with system output
        logfile = open(outfolder + "/saeno_log.txt", 'w')
        cmd = Popen(["saenopy","CONFIG","{}//config.txt".format(os.path.abspath(outfolder))], stdout=PIPE, 
                    universal_newlines=True, shell=False)
        # print and save a reduced version of saeno log
        for line in cmd.stdout:
            if not '%' in line:
                print (line, end='')
                logfile.write(str(line))
        # close again to avoid loops            
        cmd.stdout.close()            
        
    # if false just show the non reduced system output    
    else:
        cmd = call(["saenopy","CONFIG", "{}//config.txt".format(os.path.abspath(outfolder))])

    # copy result files from "*_py2" folder
    for filename in glob(outfolder+'_py2/*.*'):
        shutil.copy(filename, outfolder)

    # remove "*_py2" folder
    shutil.rmtree(outfolder+'_py2')


def distribute(func, const_args, var_arg='pressure', start=0.1, end=1000, n=120, log_scaling=True, n_cores=None):
    if n_cores is None:
        n_cores = os.cpu_count()

    if log_scaling:
        values = np.logspace(np.log10(start), np.log10(end), n, endpoint=True)
    else:
        values = np.linspace(np.log10(start), np.log10(end), n, endpoint=True)

    outfolder = const_args['outfolder']
    del const_args['outfolder']

    if not os.path.exists(outfolder):
        os.makedirs(outfolder)

    np.savetxt(outfolder+'/'+var_arg+'-values.txt', values)

    values = list(values)
    index = list(np.arange(len(values)))
    processes = []

    
    
    while True:
        if len(processes) < n_cores and len(values) > 0:
            command = '''python -c "import jointforces; import jointforces as jf; {}('''.format(func)
            for key in const_args:
                if isinstance(const_args[key], str):
                    command += '''{}='{}','''.format(key, const_args[key])
                else:
                    command += '''{}={},'''.format(key, const_args[key])
            command += '''outfolder='{}', {}={})"'''.format(outfolder+'/simulation'+str(index[0]).zfill(6), var_arg, values[0])
            print ("Simulations in the queue: "+str(len(values)))
            processes.append(Popen(command, shell='True'))
            del values[0]
            del index[0]

        sleep(1.)

        processes = [p for p in processes if p.poll() is None]

        if len(processes) == 0:
            break


def extract_deformation_curve(folder, x):
    # get simulation parameters
    with open(folder+'/parameters.txt', 'r') as f:
        lines = f.readlines()

    parameters = {}
    for line in lines:
        try:
            key, value = line.split('= ')
            value = value.split(' ')[0]
            parameters[key.strip()] = float(value.strip())
        except:
            pass

    # load coordinates
    coords = np.genfromtxt(folder + '/R.dat')

    # load displacements
    displ = np.genfromtxt(folder + '/U.dat')

    # compute binned normed displacements and normed coordinates
    u = np.sqrt(np.sum(coords ** 2., axis=1)) / (parameters['INNER_RADIUS']*10**-6)
    v = np.sqrt(np.sum(displ ** 2., axis=1)) / (parameters['INNER_RADIUS']*10**-6)

    y = np.array([np.nanmedian(v[(u >= x[i]) & (u < x[i + 1])]) for i in range(len(x) - 1)])

    # save results
    results = {'pressure': parameters['PRESSURE'], 'bins': x, 'displacements': y}
    return results


def create_lookup_table(folder, x0=1, x1=50, n=100):
    subfolders = glob(folder+'/*/')

    x = np.logspace(np.log10(x0), np.log10(x1), n+1, endpoint=True)
    x_center = 10**(0.5*(np.log10(x[1:]) + np.log10(x[:-1])))

    pressure_values = []
    displacement_curves = []

    for subfolder in tqdm(subfolders):
        res = extract_deformation_curve(subfolder, x)
        pressure_values.append(res['pressure'])
        displacement_curves.append(res['displacements'])

    pressure_values = np.array(pressure_values)
    displacement_curves = np.array(displacement_curves)

    return {'pressure': pressure_values, 'x': x_center, 'y': displacement_curves}


def create_lookup_functions(lookup_table):
    pressure = lookup_table['pressure']
    distance = lookup_table['x']
    displacement = lookup_table['y']

    log_pressure = np.log(pressure)
    log_distance = np.log(distance)

    x, y = np.meshgrid(log_distance, log_pressure)

    mask = ~(np.isnan(x) | np.isnan(y) | np.isnan(displacement))

    x = x[mask]
    y = y[mask]
    displacement = displacement[mask]

    f = LinearNDInterpolator(np.array([x, y]).T, displacement)

    def get_displacement(distance, pressure):
        return f(np.log(distance), np.log(pressure))

    f_inv = LinearNDInterpolator(np.array([x, displacement]).T, y)

    def get_pressure(distance, displacement):
        return np.exp(f_inv(np.log(distance), displacement))

    return get_displacement, get_pressure


def save_lookup_functions(get_displacement, get_pressure, outfile):
    with open(outfile, 'wb') as f:
        dill.dump((get_displacement, get_pressure), f)


def save_lookup_table(lookup_table, outfile):
    with open(outfile, 'wb') as f:
        dill.dump(lookup_table, f)


def load_lookup_functions(file):
    try:
        with open(file, 'rb') as f:
            get_displacement, get_pressure = dill.load(f)
    except:
        with open(file, 'rb') as f:
            lookup_table = dill.load(f)
            get_displacement, get_pressure = create_lookup_functions(lookup_table)
            
    return get_displacement, get_pressure


def linear_lookup_interpolator(emodulus, output_newtable="new-lin-lookup.pkl", reference_folder=None):
    """
    Create individual lookup-tables for linear materials by shifting a reference lookuptable for a linear fiber material.
    
    For linear fiber materials the following relation is used:  k0 = 6 * E_Modulus  (for a possion ration of 0.25, see Steinwachs[2015]) 
    
    Original simulation reached up to 10 000 Pa for a simulated k0 = 2364 (emodul ~394 Pa) - interpolation should be usefull for a wide range of
     emoduli - however keep in mind that there might be limits for extreme Emoduli-pressures combination due to the range of the original simulations 
     (in such a case a constant maximal or minimal-pressure value will be returned since there are no better fitting simulations)
    
    emodulus:ArithmeticError Desired Youngs modulus for which the linear lookup table is created
    reference_folder: Folder containing the reference lookup functions and reference interpolators to create
    the new look up table; By default (None) will search for the reference files automatically  
    output_newtable: name for the new reference table (needs .pkl ending)
    """
    
    # if not specified differently we find the correct reference files automatically
    if not reference_folder:
        import jointforces
        reference_folder = os.path.join(jointforces.__file__,"..","..","docs","data","linear_reference_table")
        
    # load in  reference lookuptable to for a simulated k2364 (emodul ~394 Pa) up to 10 000 Pa 
    get_displacement_ref, get_pressure_ref = load_lookup_functions(os.path.join(reference_folder,'linear-ref-functions-k2364.pkl'))

    # load in in reference interpolators
    f_ref,f_inv_ref = load_lookup_functions(os.path.join(reference_folder,'linear-interp-f-finv-k2364.pkl')) 
    
    # shift the table accordingly 
    def get_displacement_new(distance, pressure):
        return f_ref(np.log(distance), np.log(pressure)) * 2364 /(emodulus*6)

    def get_pressure_new(distance, displacement):
        return np.exp(f_inv_ref(np.log(distance), displacement)) * (emodulus*6)/2364
    
    # save the new lookup functions
    save_lookup_functions(get_displacement_new, get_pressure_new, output_newtable)

    return get_displacement_new, get_pressure_new



def plot_lookup_table(lookup_table, pressure=[0,10000], log_scale = True, distance=[2,50], linewidth=2, n_lines = 1000, save_plot = None,
                      fig_size=(5,4), figure=None, show=True):
    """
    Create a figure of your (.pkl) material lookuptable
    
    lookup_table: path to .pkl file
    pressure: pressure range which will be plotted as [pressure_min,pressure_max] in pascal (use the range that was used to calculate the lookuptable)
    log_scale: plot logarithmically or linear if set to False
    distance: The distance to the spheroid which is plotted on the x-axis; Unit is spheroid radii
    linewidth: linewidth for the individual plots
    n_lines: number of lines plotted between the minimal and maximal pressure
    save_plot: saves the plot as png file if a path is given 
    """
    
    import matplotlib.cm as cm 
    import matplotlib.colors as colors
    import matplotlib.pyplot as plt
    
    # load lookup table
    get_displacement, get_pressure = load_lookup_functions(lookup_table)
    
    ## to avoid scaling issue for log(0) 
    if pressure[0] == 0:
        pressure[0] = 0.01     
   
    # define pressure range
    if log_scale:
        pressure_list = np.logspace(np.log(pressure[0])/np.log(10), np.log(pressure[1])/np.log(10),  num=n_lines, base=10)     
    else:
        pressure_list = np.arange(pressure[0],pressure[1], step = (pressure[1]-pressure[0])/n_lines )
      
       
      
    # make cmaps
    mycmap = cm.get_cmap('viridis')
    mynorm = colors.LogNorm(vmin=np.min(pressure_list),vmax=np.max(pressure_list))   
    c = mycmap(mynorm(pressure_list))

    # make a colorbar
    sm = plt.cm.ScalarMappable(cmap=mycmap, norm= mynorm)
    sm.set_array(c)
        
    # create distance list
    distance_list = np.arange(distance[0],distance[1], step =(distance[1]-distance[0]) / 1000 )
    
    
    # get displacements for pressures list
    displacement_list = [get_displacement(distance_list,i) for i in pressure_list]

    # create a figure if no figure is specified
    if figure is None:
        figure = plt.figure(figsize=fig_size)

    for i in range(len(displacement_list)):
        plt.plot( distance_list , displacement_list[i], c= c[i],linewidth=linewidth,alpha=0.5)
   
    # set x,y limit - go a bit above minimum and maximum
    plt.ylim(np.nanmin(displacement_list)-(0.1*np.nanmin(displacement_list)),
             np.nanmax(displacement_list)+(0.1*np.nanmax(displacement_list)))
    
    # set log if activated and go even more above minimum/maximum in log scale
    if log_scale:
        plt.loglog()
        plt.ylim(np.nanmin(displacement_list)-(0.8*np.nanmin(displacement_list)),
             np.nanmax(displacement_list)+(0.8*np.nanmax(displacement_list)))
    
    plt.grid(False) 
    plt.ylabel('Normalized matrix deformations (u/r₀)', fontsize=10)           # plt.ylabel('Normed Deformation', fontsize=14)
    plt.xlabel('Normalized distance from organoid center (r/r₀)', fontsize=10) # plt.xlabel('Normed Distance', fontsize=14)
    
    # make a colorbar
    cbar= plt.colorbar(sm, )
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel('Pressure (Pa)', rotation=270 , fontsize=14)
    plt.tight_layout()
    # save plot if specified
    if save_plot is not None:
        plt.savefig(save_plot,dpi=500)  
    # show figure if specified
    if show==True:
        plt.show()
    return figure


def plot_lookup_data(lookup_table, data_folder, timesteps=[10,30,60], distance=[2,50], linewidth=2, 
              color_line="k", color_raw="r", scatter_raw_data = True, marker_size_scatter=0.1,
              marker_size_mean=10,angle_filter=20, color_list = None,  plot_means = True,
              label_list = None, timesteps_scatter=None , result_file_name ="result.xlsx"): 
    """
    plot the pressure for certain timesteps into the lookup table as a grey line;
    scatter_raw_data option to vuisualize all deformtation-distance data at that timestep
    
    Use the function after calling "plot_lookup_table" to visualize a certain pressure within the created lookuptable
    
    lookup_table: path to .pkl filer
    data_folder: path to the evaluated folder containing result.xlsx and the dis*.npy & seg*.npy data
    timesteps: list of the timesteps to plot into the lookup function
    scatter_raw_data: option to scatter the individual deformations
    timesteps_scatter: might be used to scatter only several timesteps - if none identical to timesteps
    angle_filter: use same value as in evaluation
    if color_list & label_list is provided the given colors and labels are used when plotting the raw data
    e.g color_list=["C0","C1","C2"], label_list=["1h","3h","12h"],
    
    timesteps_scatter: might be used to scatter only several timesteps - if timesteps_scatter=None identical to timesteps;
    keep the same length as timesteps to asure same coloring :
    e.g timesteps=[1,2,3] --> timesteps_scatter= [None,2,None]  to show only second element
    """   
    import matplotlib.pyplot as plt
    from glob import glob
    from natsort import natsorted
    from .simulation import load_lookup_functions
    from .utils import load
    from .force import infer_pressure
    from tqdm import tqdm
    
    # scatter same timesteps if not spcified
    if not timesteps_scatter:
        timesteps_scatter = timesteps
    
    # read in the data
    pressure = pd.read_excel(os.path.join(data_folder,result_file_name))["Mean Pressure (Pa)"]
    
    # load lookup table
    get_displacement, get_pressure = load_lookup_functions(lookup_table)
    
    # define pressure range  # -1 since pressure list is 1 element shorter 
    # compared to image list (differences are taken)
    pressure_list = [pressure[i-1] for i in timesteps]
    print (pressure_list)
    
    # create distance list
    distance_list = np.arange(distance[0],distance[1], step =(distance[1]-distance[0]) / 1000 )

    # get displacements for pressures list
    displacement_list = [get_displacement(distance_list,i) for i in pressure_list]
   
    # draw simulations in uniform color
    for i in range(len(displacement_list)):
       plt.plot( distance_list , displacement_list[i], c= color_line, #linestyle="--",
                linewidth=linewidth,alpha=0.5,zorder=30)#, label="Simulation") 
       if i==0: # plot label once
           plt.plot( [] ,[] , c= color_line, #linestyle="--",
                    linewidth=linewidth,alpha=0.5,zorder=30, label="Simulations") 
        
    # read in segemntation
    seg_files = natsorted(glob(data_folder+'/seg*.npy'))[:np.nanmax(timesteps)]
    
    
    # load deformations
    # look for accumulated deformation files (new standard)
    d_accumulated_files = natsorted(glob(data_folder+'/def*.npy'))[:np.nanmax(timesteps)]  # do not calcualte more then the necessary time steps
    # look also for not-accummulated deformations (old standard)
    d_notaccumulated_files = natsorted(glob(data_folder+'/dis*.npy'))[:np.nanmax(timesteps)]  # do not calcualte more then the necessary time steps
    # if not-accumulated deformations are found chose different mode
    if len(d_notaccumulated_files) > len(d_accumulated_files):
        accumulated = False
        dis_files = d_notaccumulated_files
        print("Found not-accumulated deformation files (old standard) and will conduct calculations accordingly.")
    # else do the calcualtion with accumulated d eformations already
    else:
        accumulated = True
        dis_files = d_accumulated_files


    # initial spheroid radius and surface (used for force estimation)
    r0 = load(seg_files[0])['radius']
   
    u_sum = None
    v_sum = None
    
    distance_list_raw = []
    displacement_list_raw = []
    pressure_list_raw = []
    
    # loop over series of PIV results
    for (dis_file, seg_file) in tqdm(zip(dis_files, seg_files)):
        dis = load(dis_file)
        seg = load(seg_file)

        x_rav = np.ravel(dis['x'])
        y_rav = np.ravel(dis['y'])
        #print("data_points:"+str(len(x_rav)))
              
        # get deformations
        # sum up if we have not-accummulated deformations (old standard)
        if accumulated == False:
            try:
                u_sum += np.ravel(dis['u'])
                v_sum += np.ravel(dis['v'])
            except:
                u_sum = np.ravel(dis['u'])
                v_sum = np.ravel(dis['v'])
        # else read in accummulated deformations directly (new standard)
        else:
            u_sum = np.ravel(dis['u'])
            v_sum = np.ravel(dis['v'])    
            
 
        cx, cy = seg['centroid']
        distance_raw, displacement_raw, angle_raw, pressure_raw = infer_pressure(x_rav, y_rav, u_sum, v_sum, cx, cy, r0, get_pressure , angle_filter=angle_filter)
        #print (len(distance_raw)) # length of valid datapoints after angle filter      
          
        # create list with accumulated deformations 
        distance_list_raw.append(distance_raw)
        displacement_list_raw.append(displacement_raw)
        pressure_list_raw.append(pressure_raw)
        
        
    # SCATTERED RAW DATA
    #now plot the accumulated raw data at the corresponding timepoints if activated;
    # deformation and image data have a length difference of 1
    if scatter_raw_data:
        # in case colors and labels are defined; 
        if color_list and label_list:
          for ci,t in enumerate(timesteps_scatter):
            if t == None: #do not plot None elements 
                  continue
            plt.scatter(distance_list_raw[t-1],displacement_list_raw[t-1],s=marker_size_scatter,zorder=20,c = color_list[ci])      
        
            ## calculate CoV of pressure values for 72h and add to plot
            # if t == 3*70-2: #3*72:
            #     mask = (distance_list_raw[t-1]>=5)&(distance_list_raw[t-1]<=10)
                
            #     data_masked = pressure_list_raw[t-1][mask]
            #     CoV = np.std(data_masked)/np.nanmean(data_masked)
            #     print (len(data_masked))
            #     print (len(distance_list_raw[t-1]))
            #     print (CoV)
                
            #     plt.text(0.05,0.93,s=f"CoV (72h; 5-10r): {np.around(CoV,3)}",
            #              transform=plt.gca().transAxes, zorder=200,
            #              fontsize=11,c="darkred")
            
        
        # else same color for all here
        else: 
           for ci,t in enumerate(timesteps_scatter):
            if t == None: #do not plot None elements
                  continue
            plt.scatter(distance_list_raw[t-1],displacement_list_raw[t-1],s=marker_size_scatter,zorder=20,c = color_raw) # deformation and image data have a length difference of 1
        
        
    # PLOT THE MEAN RAW DATA
    if color_list and label_list:
        for ci,t in enumerate(timesteps):
            if plot_means == False:
                continue
            # calculate the mean in distance windows for timesteps
            mean_distance = []
            mean_displacement = []
            for i in range(distance[0],int(np.max(distance_list_raw[t-1]))):
                mean_distance.append(i+0.5)
                mean_disp = np.nanmean(displacement_list_raw[t-1][(distance_list_raw[t-1]>=i) & (distance_list_raw[t-1]<i+1)])
                mean_displacement.append(mean_disp)
            plt.plot(mean_distance,mean_displacement,"o-",ms=marker_size_mean, zorder=2000, 
                     linewidth=linewidth,  markerfacecolor="w",  markeredgecolor=color_list[ci], 
                     markeredgewidth=2 , label=label_list[ci], color="w" )     
        ax = plt.legend( markerscale=0.4, fontsize=6.5,loc="upper right")
        ax.set_zorder(2000)
    else: 
         for ci,t in enumerate(timesteps):
            if plot_means == False:
              continue
            # calculate the mean in distance windows for timesteps
            mean_distance = []
            mean_displacement = []
            for i in range(distance[0],int(np.max(distance_list_raw[t-1]))):
                mean_distance.append(i+0.5)
                mean_disp = np.nanmean(displacement_list_raw[t-1][(distance_list_raw[t-1]>=i) & (distance_list_raw[t-1]<i+1)])
                mean_displacement.append(mean_disp)
            plt.plot(mean_distance,mean_displacement,"o-",ms=marker_size_mean,zorder=2000,
                     linewidth=linewidth, markerfacecolor="w",color ="w",
                     markeredgecolor=color_raw, markeredgewidth=2)   
        #plt.xlim(2,20);plt.ylim(1e-4,3) 

    return 
