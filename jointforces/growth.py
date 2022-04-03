import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
from natsort import natsorted
from scipy.ndimage.morphology import distance_transform_edt
from .utils import load
import os
import matplotlib

def eval_angular_growth(folder, outfile=None):
    """
    Evaluates the growth in 5° angular sections (parelell to the agular pressure analyis) and the total growth of the 
    spheroid using the continuously tracked segmentation (therfore the "continous_segmentation" must be set to True
    in the function "compute_displacement_series")
    
    
    - folder: contains PIV results with segmentation
    - outfile: can be used to change the standard name of the excel file - if none excel files are created 
    in the given folder
    
    Returns:
    saves an excel file with the angular growth (angular area section at each timepoint / dived by 
                                                  angular area at timepoint 0)
    and another excel file with the total growth(total segmented area  at each timepoint / dived by 
                                                  total segmented area at timepoint 0)
    
    """
       
    ## load segmentations
    seg_files = natsorted(glob(folder+'/seg*.npy'))
       
    # initial spheroid mask 
    mask_t0 = load(seg_files[0])['mask']
    # intitial spheoir center for all further analysis
    cx_t0, cy_t0 = load(seg_files[0])['centroid']
    
    # Calculate Angle + Distances 
    y,x = np.indices(mask_t0.shape)
    dx = x - cx_t0
    dy = y - cy_t0
    #distance = np.sqrt(dx ** 2 + dy ** 2)  # dist to center
    angle = np.arctan2(dy, dx) 
           
    # initialize result dictionary
    results = {'total_growth': []}
    
    # create dict with all angles
    angles_dict = {str(a): [] for a in range(-180, 180, 5)}
    
      
    # loop over series of PIV segmentations
    for seg_file in tqdm(seg_files):
        # load current segmentation
        mask = load(seg_file)['mask']      
        ## total growth
        results['total_growth'].append(np.sum(mask)/np.sum(mask_t0))            
        # growth in angular section
        angular_sections = []
        growth_angular_sections = []
        for alpha in range(-180, 180, 5):
            mask2 = (angle >= (alpha-5)*np.pi/180.) & (angle < (alpha+5)*np.pi/180.)
            
           
            angular_sections.append(alpha)
            growth_angular_sections.append(np.nansum(mask[mask2])/np.nansum(mask_t0[mask2]))  #### growth per angle section
                 
        # append pressures for all angle data
        for i,a in enumerate(angles_dict):
            angles_dict[a].append(growth_angular_sections[i])
    
    ##### just for testing angular sections
    # plt.figure()    # check the segmentation
    # plt.imshow(mask_t0)
    # plt.scatter(cx_t0, cy_t0)
    # plt.show()      
    # plt.figure()    # check the angular section  
    # alpha= -70
    # mask2 = (angle >= (alpha-5)*np.pi/180.) & (angle < (alpha+5)*np.pi/180.)
    # plt.imshow(mask2)
    # plt.scatter(cx_t0, cy_t0)
    # plt.show()
    
    df = pd.DataFrame.from_dict(results)
    df.columns = ['total_growth',
                  ]
    if outfile is not None:
        df.to_excel(outfile+'_total_growth.xlsx')
    else:
        df.to_excel(folder+'//total_growth.xlsx')
     
    # save pressures for all angles
    an = pd.DataFrame.from_dict(angles_dict)
    an.columns = [a for a in range(-180, 180, 5)]
    if outfile is not None:
        an.to_excel(outfile[:-5]+'_angular_growth.xlsx')
    else:
        an.to_excel(folder+'//angular_growth.xlsx')    
    
    return



# Evaluate Angle dependet Pressures
def plot_growth(folder,  output, n_max=None,file_angulargrowth=None, file_totalgrowth=None, fontsize=7, dt=None, angle_legend=True):
    """
    Evaluate the angular growth over time for an spheroid (needs folder with PIV analysis in which also results from 'eval_angular_growth' are saved) and stores 
    the results in the output folder.
    
    - N_max may define the last image to evaluate 
    - Evaluated excel sheet with angular data is automatically searched in folder, if name differs from 'total_growth.xlsx' and 'angular_growth.xlsx' the
    user defined names can be specified using 'file_angulargrowth' or 'file_totalgrowth'
    - dt is time between consecutive images in seconds; If given a timestamp will be displayed 
    - angle_legend will draw a legend with angles if active 
    """
    from datetime import timedelta
    
   
    # read in total growth 
    if file_totalgrowth is None:   # by default looks in the PIV folder for default names
        total_growth  = pd.read_excel(folder + '//total_growth.xlsx')
    else:
        total_growth  = pd.read_excel(file_totalgrowth)  ## if stored somewhere else user can define the path
    # read in angular growth
    if file_angulargrowth is None:
        angular_growth  = pd.read_excel(folder + '//angular_growth.xlsx')
    else:
        angular_growth  = pd.read_excel(file_angulargrowth)
                
    # read in plots in folder
    plots = glob(folder+'//plot*.png')
    # make folder if not existing
    if not os.path.exists(output):
        os.makedirs(output)
    # list of evaluated angles
    angle_list = [a for a in range(-180, 180, 5)]
    # make figure
    f, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,5))
            
    #second plot with only the left plot
    fb,axb = plt.subplots( figsize=(3,3))  
    # for color norm
    max_a = np.nanmax([np.nanmax(np.array(angular_growth[z])[:n_max]) for z in angle_list])  
    plt.style.use('dark_background')
    # make result dictionary


    # loop over all timesteps
    # -1 as the deformation list is one entry less (between images) compared to
    # segmentation list
    for t in tqdm(range(len(angular_growth[0][:n_max])-1)):    


              
        ang_g = []
        # go through different angles 
        for z in angle_list:        
             # append for each angle the pressure of the current timestep
             ang_g.append(np.array(angular_growth[z])[t])    
                
        # norm colors for each image individually    viridis/RdBu      
        colors = matplotlib.cm.get_cmap('viridis')(ang_g / max_a ) 
        colorcycle = 'darkred'
    
        # do some statistics
        mean = np.nanmean(ang_g)
        sd = np.nanstd(ang_g)
        CoV = np.round(sd/mean, 2)

        if angle_legend == True:
				# position of label text
            shift = 0.9
            posx_legend = [np.cos(i*np.pi/180) * shift for i in angle_list] 
            posy_legend = [np.sin(i*np.pi/180) * shift for i in angle_list] 
				# plot text
            [ax1.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=4, horizontalalignment='center') for i in range(len(angle_list))] 
            [axb.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=3.6, horizontalalignment='center') for i in range(len(angle_list))]   
						
 			
		# now use 5 (500% growth) as maximum
        x = [np.cos(i*np.pi/180) for i in angle_list] * (ang_g / np.array([5] *len(ang_g)) ) 
        y = [np.sin(i*np.pi/180) for i in angle_list] * (ang_g / np.array([5] *len(ang_g)) ) 
				# combined figure
				# plot circles
        circle_1 = plt.Circle((0, 0), 1/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_2 = plt.Circle((0, 0), 2/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_3 = plt.Circle((0, 0), 3/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_4 = plt.Circle((0, 0), 4/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        ax1.text(0.5, 0.57, '1', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.67, '2', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.77, '3', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.87, '4', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.add_artist(circle_1)
        ax1.add_artist(circle_2)
        ax1.add_artist(circle_3)
        ax1.add_artist(circle_4)
        ax1.axis('equal')  
        ax1.scatter(x,y, s=22, c=colors)
        ax1.set_xlim([-1, 1])
        ax1.set_ylim([-1, 1])
        ax1.axis('off')
		

        # annotate global growth    
        ax1.text(0.01, 0.93,'Total Growth: '+str(np.round(total_growth["total_growth"][t],2)), 
				 horizontalalignment='center',
				 verticalalignment='center',
				 transform = ax1.transAxes, fontsize=7)  
            
        # annotate CoV    
        ax1.text(0.01, 0.91,'Coefficient of Variation: '+str(CoV), 
				 horizontalalignment='center',
				 verticalalignment='center',
				 transform = ax1.transAxes, fontsize=7)  
        
           
		# annotate time if timestep is given
        if dt is not None:
            ax1.text(0.03, 0.98,str(timedelta(seconds=t*dt)), 
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform = ax1.transAxes, fontsize=11)  

					
		# show quiver plot on the right side
        plot = plt.imread(plots[t])
        ax2.imshow(plot)
        plt.axis('off')
        ax2.axis('off')
				# save figure
        f.savefig(output+'//plot_{}.png'.format(str(t).zfill(4)), dpi=200)
        ax1.cla()
        ax2.cla()
		# Single figure
		# plot circles
        circle_1 = plt.Circle((0, 0), 1/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_2 = plt.Circle((0, 0), 2/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_3 = plt.Circle((0, 0), 3/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_4 = plt.Circle((0, 0), 4/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        axb.text(0.5, 0.57, '1', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=axb.transAxes,fontsize=fontsize)
        axb.text(0.5, 0.67, '2', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=axb.transAxes,fontsize=fontsize)
        axb.text(0.5, 0.77, '3', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=axb.transAxes,fontsize=fontsize)
        axb.text(0.5, 0.87, '4', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=axb.transAxes,fontsize=fontsize)
        axb.add_artist(circle_1)
        axb.add_artist(circle_2)
        axb.add_artist(circle_3)
        axb.add_artist(circle_4)
        axb.axis('equal')  
        axb.scatter(x,y, s=22, c=colors)
        axb.set_xlim([-1, 1])
        axb.set_ylim([-1, 1])
        axb.axis('off')
        plt.savefig(output+'//single_{}.png'.format(str(t).zfill(4)), dpi=200)
        axb.cla()
       

    return 



# Evaluate Angle dependet Pressures
def plot_growth_and_pressure(folder,  output, n_max=None,file_angulargrowth=None, file_totalgrowth=None, small_pressure=False, file_angularpressure=None, fontsize=7, dt=None, angle_legend=False):
    """
    Evaluate the angular growth over time for an spheroid (needs folder with PIV analysis in which also results from 'eval_angular_growth' are saved) and stores 
    the results in the output folder.
    
    - N_max may define the last image to evaluate 
    - Evaluated excel sheet with angular data is automatically searched in folder, if name differs from 'total_growth.xlsx' and 'angular_growth.xlsx' the
    user defined names can be specified using 'file_angulargrowth' or 'file_totalgrowth'
    - dt is time between consecutive images in seconds; If given a timestamp will be displayed 
    - angle_legend will draw a legend with angles if active 
    """
    from datetime import timedelta
    
   
    # read in total growth 
    if file_totalgrowth is None:   # by default looks in the PIV folder for default names
        total_growth  = pd.read_excel(folder + '//total_growth.xlsx')
    else:
        total_growth  = pd.read_excel(file_totalgrowth)  ## if stored somewhere else user can define the path
    # read in angular growth
    if file_angulargrowth is None:
        angular_growth  = pd.read_excel(folder + '//angular_growth.xlsx')
    else:
        angular_growth  = pd.read_excel(file_angulargrowth)
        
    # read in angular pressure
    if file_angularpressure is None:
        angular_pressure  = pd.read_excel(folder + '//result_angles.xlsx')
    else:
        angular_pressure  = pd.read_excel(file_angularpressure)

    
    # read in plots in folder
    plots = glob(folder+'//plot*.png')
    # make folder if not existing
    if not os.path.exists(output):
        os.makedirs(output)
    # list of evaluated angles
    angle_list = [a for a in range(-180, 180, 5)]
    # make figure
    f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(9+4.5,5))
            
    #second plot with only the left plot
    fb,axb = plt.subplots( figsize=(3-1,3-1))  
    # for color norm
    max_a = np.nanmax([np.nanmax(np.array(angular_growth[z])[:n_max]) for z in angle_list])  
    plt.style.use('dark_background')
    
    
    # read in pressures
    # for pressure color norm
    max_pr = np.nanmax([np.nanmax(np.array(angular_pressure[z])[:n_max]) for z in angle_list])  


    # loop over all timesteps
    for t in tqdm(range(len(angular_pressure[0][:n_max]))):   ### angular pressure has one entry less then angular growth since its calculated between images and segementation not

         
        ang_g = []
        # go through different angles 
        for z in angle_list:        
             # append for each angle the pressure of the current timestep
             ang_g.append(np.array(angular_growth[z])[t])    
                
        # norm colors for each image individually    viridis/RdBu      
        colors = matplotlib.cm.get_cmap('viridis')(ang_g / max_a ) 
        colorcycle = 'darkred'
    
        # do some statistics
        mean = np.nanmean(ang_g)
        sd = np.nanstd(ang_g)
        CoV = np.round(sd/mean, 2)

        if angle_legend == True:
				# position of label text
            shift = 0.9
            posx_legend = [np.cos(i*np.pi/180) * shift for i in angle_list] 
            posy_legend = [np.sin(i*np.pi/180) * shift for i in angle_list] 
				# plot text
            [ax1.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=4, horizontalalignment='center') for i in range(len(angle_list))] 
            [axb.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=3.6, horizontalalignment='center') for i in range(len(angle_list))]   
						
 			
		# now use 5 (500% growth) as maximum
        x = [np.cos(i*np.pi/180) for i in angle_list] * (ang_g / np.array([5] *len(ang_g)) ) 
        y = [np.sin(i*np.pi/180) for i in angle_list] * (ang_g / np.array([5] *len(ang_g)) ) 
				# combined figure
				# plot circles
        circle_1 = plt.Circle((0, 0), 1/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_2 = plt.Circle((0, 0), 2/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_3 = plt.Circle((0, 0), 3/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        circle_4 = plt.Circle((0, 0), 4/ 5 , color=colorcycle, zorder=10, fill=False, linestyle='--')
        ax1.text(0.5, 0.57, '1', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.67, '2', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.77, '3', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.text(0.5, 0.87, '4', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax1.transAxes,fontsize=fontsize)
        ax1.add_artist(circle_1)
        ax1.add_artist(circle_2)
        ax1.add_artist(circle_3)
        ax1.add_artist(circle_4)
        ax1.axis('equal')  
        ax1.scatter(x,y, s=22, c=colors)
        ax1.set_xlim([-1, 1])
        ax1.set_ylim([-1, 1])
        ax1.axis('off')
		

        # annotate global growth    
        ax1.text(0.01, 0.93,'Total Growth: '+str(np.round(total_growth["total_growth"][t],2)), 
				 horizontalalignment='center',
				 verticalalignment='center',
				 transform = ax1.transAxes, fontsize=7)  
            
        # annotate CoV    
        ax1.text(0.01, 0.91,'Coefficient of Variation: '+str(CoV), 
				 horizontalalignment='center',
				 verticalalignment='center',
				 transform = ax1.transAxes, fontsize=7)  
        
           
		# annotate time if timestep is given
        if dt is not None:
            ax1.text(0.03, 0.98,str(timedelta(seconds=t*dt)), 
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform = ax1.transAxes, fontsize=11)  

					
		### now plot presures as second plot
        pressures = []
        # go through different angles 
        for z in angle_list:        
             # append for each angle the pressure of the current timestep
             pressures.append(np.array(angular_pressure[z])[t])    
                
        # norm colors for each image individually    viridis/RdBu      
        colors = matplotlib.cm.get_cmap('viridis')(pressures / max_pr ) 
        colorcycle = 'darkred'
    
        # do some statistics
        mean_pr =np.nanmean(pressures)
        sd_pr = np.nanstd(pressures)
        CoV_pr = np.round(sd_pr/mean_pr, 2)
        # plot pressure
        if angle_legend == True:
            # position of label text
            shift = 0.9
            posx_legend = [np.cos(i*np.pi/180) * shift for i in angle_list] 
            posy_legend = [np.sin(i*np.pi/180) * shift for i in angle_list] 
            # plot text
            [ax1.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=4, horizontalalignment='center') for i in range(len(angle_list))] 
            [axb.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=3.6, horizontalalignment='center') for i in range(len(angle_list))]   
            [ax2.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=4, horizontalalignment='center') for i in range(len(angle_list))] 
            [axb.text(posx_legend[i],posy_legend[i], str(angle_list[i]), color='darkred', fontsize=3.6, horizontalalignment='center') for i in range(len(angle_list))]   
		
        
        # Mode that shows and norms up to 1000 Pa
        if small_pressure == False:

            # now use 1000 as maximum
            x = [np.cos(i*np.pi/180) for i in angle_list] * (pressures / np.array([1000]*len(pressures)) ) 
            y = [np.sin(i*np.pi/180) for i in angle_list] * (pressures / np.array([1000] *len(pressures)) ) 
            # combined figure
            # plot circles
            circle_100pa = plt.Circle((0, 0), 100/ 1000 , color=colorcycle, zorder=10, fill=False, linestyle='--')
            circle_250pa = plt.Circle((0, 0), 250/ 1000 , color=colorcycle, zorder=10, fill=False, linestyle='--')
            circle_500pa = plt.Circle((0, 0), 500/ 1000 , color=colorcycle, zorder=10, fill=False, linestyle='--')
            circle_750pa = plt.Circle((0, 0), 750/ 1000 , color=colorcycle, zorder=10, fill=False, linestyle='--')
            ax2.text(0.5, 0.57, '100 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
            ax2.text(0.5, 0.65, '250 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
            ax2.text(0.5, 0.77, '500 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
            ax2.text(0.5, 0.9, '750 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
            ax2.add_artist(circle_100pa)
            ax2.add_artist(circle_250pa)
            ax2.add_artist(circle_500pa)
            ax2.add_artist(circle_750pa)
            ax2.axis('equal')  
            ax2.scatter(x,y, s=22, c=colors)
            ax2.set_xlim([-1, 1])
            ax2.set_ylim([-1, 1])
            ax2.axis('off')
            # annotate CoV    
            ax2.text(0.01, 0.91,'Coefficient of Variation: '+str(CoV_pr), 
                     horizontalalignment='center',
                     verticalalignment='center',
                     transform = ax2.transAxes, fontsize=7)  

        # mode that shows and norms up to 100 Pa
        if small_pressure == True:
		# now use 100 as maximum
              x = [np.cos(i*np.pi/180) for i in angle_list] * (pressures / np.array([100]*len(pressures)) ) 
              y = [np.sin(i*np.pi/180) for i in angle_list] * (pressures / np.array([100] *len(pressures)) ) 
        		# combined figure
        		# plot circles
              circle_a = plt.Circle((0, 0), 10/ 100 , color=colorcycle, zorder=10, fill=False, linestyle='--')
              circle_b = plt.Circle((0, 0), 25/ 100 , color=colorcycle, zorder=10, fill=False, linestyle='--')
              circle_c = plt.Circle((0, 0), 50/ 100 , color=colorcycle, zorder=10, fill=False, linestyle='--')
              ax2.text(0.5, 0.57, '10 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
              ax2.text(0.5, 0.65, '25 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
              ax2.text(0.5, 0.77, '50 Pa', color = colorcycle ,horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes,fontsize=fontsize)
              ax2.add_artist(circle_a)
              ax2.add_artist(circle_b)
              ax2.add_artist(circle_c)
              ax2.axis('equal')  
              ax2.scatter(x,y, s=22, c=colors)
              ax2.set_xlim([-1, 1])
              ax2.set_ylim([-1, 1])
              ax2.axis('off')
		      # annotate CoV    
              ax2.text(0.01, 0.91,'Coefficient of Variation: '+str(CoV), 
                       horizontalalignment='center',
                           verticalalignment='center',
                          transform = ax2.transAxes, fontsize=9.2)  

        
        
        
        # show quiver plot on the right side
        plot = plt.imread(plots[t])
        ax3.imshow(plot)
        plt.axis('off')
        ax3.axis('off')
				# save figure
        f.savefig(output+'//plot_{}.png'.format(str(t).zfill(4)), dpi=200)
        ax1.cla()
        ax2.cla()
        ax3.cla()

    return 