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


def spherical_contraction(meshfile, outfolder, pressure, material, r_inner=None, r_outer=None, logfile=False,  max_iter = 300, step = 0.066, conv_crit = 0.01):
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


def load_lookup_functions(file):
    with open(file, 'rb') as f:
        get_displacement, get_pressure = dill.load(f)
    return get_displacement, get_pressure



def linear_lookup_interpolator(emodulus, output_newtable="new-lin-lookup.pkl", reference_folder=None):
    """
    Create individual lookup-tables for linear materials by shifting a reference lookuptable for a linear fiber material.
    
    For linear fiber materials the following relation is used:  k0 = 6 * E_Modulus (see Steinwachs[2015]) 
    
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




