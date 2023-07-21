import os
import gc
import shutil
import random
import numpy as np
import pandas as pd
import xml.etree.ElementTree as et

from matplotlib import pyplot as plt

import ppls1.imp.chp as imp
import ppls1.exp.chp as exp


def rho_droplet_vrabec2006(T, rad):
    '''
    Get saturated densities of DROPLETS by Vrabec et al., Molecular Physics 104 (2006). Equation numbers refer this paper.
    :param float T: Temperature
    :param float rad: Radius of droplet
    :return: float rhol, float rhov: Saturated liquid and vapor density
    '''
    
    tc,rc=1.0779,0.3190
    dt=tc-T
    a,b,c=0.5649,0.1314,0.0413
    rhol=rc+a*dt**(1/3.)+b*dt+c*dt**(3/2.)       # equation 4
    a,b,c=0.5649,0.2128,0.0702
    rhov=rc-a*dt**(1/3.)+b*dt+c*dt**(3/2.)       # equation 5
    
    a,b,c,d,e=0.1485,0.0471,1.272,0.4944,5.0090
    rhol=rhol+(a-b*dt**2)/(rad+c-d*dt**-1+e*dt)  # equation 18
    
    a,b,c,d,e=0.0541,0.1682,1.0810,0.2490,10.87
    rhov=rhov+(-a+b*T**2)/(rad-c-d*dt**-1+e*dt)  # equation 19
    
    return rhol,rhov


def calc_radius(chpTemp, box):
    '''
    Calculate radius of droplet based on checkpoint AS DATAFRAME; Radius is at position of inflection of density profile

    :param dataframe chpTemp: Checkpoint with one droplet in center
    :param list box: List of box sizes (x,y,z)
    :return rad: Radius of droplet
    :return profile: Density profile of droplet and vapor until 2*rad
    '''
    
    minBox = min(box[0],box[1],box[2])
    
    # Calculate distance of particle from center
    rx_dist = chpTemp['rx'] - 0.5*box[0]
    ry_dist = chpTemp['ry'] - 0.5*box[1]
    rz_dist = chpTemp['rz'] - 0.5*box[2]
    chpTemp['distance'] = (rx_dist*rx_dist + ry_dist*ry_dist + rz_dist*rz_dist)**0.5
    
    # Count particles within bins
    binWidth = 0.5
    bins = pd.cut(chpTemp['distance'], np.arange(0, minBox, binWidth))
    profile = chpTemp.groupby(bins)['distance'].count().to_frame()
    profile.rename(columns={'distance': 'numParts'}, inplace=True)
    
    # Calculate mean radius of bins
    profile['radius'] = 0.0
    for idx in profile.index:
        profile.loc[idx,'radius'] = idx.mid
    
    # Volume of bin/Kugelschale
    profile['binVolume'] = (4/3)*np.pi*((profile['radius']+0.5*binWidth)**3 - (profile['radius']-0.5*binWidth)**3)
    
    profile['rho'] = profile['numParts']/profile['binVolume']
    
    profile.set_index('radius', drop=True, inplace=True)
    
    profile = profile.loc[0:0.5*minBox]
    
    rad = profile.iloc[10:-10]['rho'].diff().abs().idxmax() - 0.5*binWidth
    
    return rad, profile


def calc_radius_chp(in_file_path):
    '''
    Calculate radius of droplet based on file path

    :param str in_file_path: Path to input checkpoint (binary)
    :return rad: Radius of droplet
    :return profile: Density profile of droplet and vapor until 2*rad
    '''
    
    in_file_path_header = in_file_path[:-4]+'.header.xml'
    
    # Read in checkpoint header
    headerXMLTree = et.parse(in_file_path_header)
    headerXML = headerXMLTree.getroot()
    
    # Read in checkpoint data
    chp = imp.imp_chp_bin_DF(in_file_path)
    
    # Get box lengths from xml header file
    xBoxLength = float(headerXML.find('headerinfo/length/x').text)
    yBoxLength = float(headerXML.find('headerinfo/length/y').text)
    zBoxLength = float(headerXML.find('headerinfo/length/z').text)
    
    # Check if box is cubic
    if xBoxLength == yBoxLength == zBoxLength:
        pass
    else:
        print('Warning! Box is not cubic!')
    
    # Calculate real radius
    rad, profile = calc_radius(chp, [xBoxLength, yBoxLength, zBoxLength])
    
    return rad, profile


def chp_replicate_cut(in_file_path, out_file_path):
    '''
    Function to generate a binary checkpoint including header which is num_reps times bigger in x, y and z direction than the input
    Additionally, cut out droplet

    :param str in_file_path: Path to input checkpoint (binary)
    :param str out_file_path: Path and filename of output checkpoint (binary)
    :return: 1 if successful
    '''
    
    in_file_path_header = in_file_path[:-4]+'.header.xml'
    out_file_path_header = out_file_path[:-4]+'.header.xml'
    
    # Read in checkpoint header
    headerXMLTree = et.parse(in_file_path_header)
    headerXML = headerXMLTree.getroot()
    
    # Read in checkpoint data
    dfChp = imp.imp_chp_bin_DF(in_file_path)
    
    # Get box lengths from xml header file
    xBoxLength = float(headerXML.find('headerinfo/length/x').text)
    yBoxLength = float(headerXML.find('headerinfo/length/y').text)
    zBoxLength = float(headerXML.find('headerinfo/length/z').text)

    lengthChp = len(dfChp)
    chpBigList = list()
    # Copy df num_reps*num_reps*num_reps times to create the large chp
    for xx in range(num_reps):
        for yy in range(num_reps):
            for zz in range(num_reps):
                chpBigList.append(dfChp)
    dfChpBig = pd.concat(chpBigList)
    # Free memory
    del dfChp
    del chpBigList
    gc.collect()  # Garbage collection to free memory
    
    # Modify positions
    for xx in range(num_reps):
        for yy in range(num_reps):
            for zz in range(num_reps):
                idx = lengthChp*(xx + yy*num_reps + zz*num_reps**2)
                dfChpBig.iloc[idx:idx+lengthChp]['rx'] += xx*xBoxLength
                dfChpBig.iloc[idx:idx+lengthChp]['ry'] += yy*yBoxLength
                dfChpBig.iloc[idx:idx+lengthChp]['rz'] += zz*zBoxLength
    
    xBoxLength *= num_reps
    yBoxLength *= num_reps
    zBoxLength *= num_reps
    
    # Cutout droplet
    chpDropList = []
    dfChpBig['distance'] = ((dfChpBig['rx']-0.5*xBoxLength)**2 + (dfChpBig['ry']-0.5*yBoxLength)**2 + (dfChpBig['rz']-0.5*zBoxLength)**2)**0.5
    chpDropList.append(dfChpBig[dfChpBig['distance']<=radius_shrink*radius])  # Inside droplet --> Liquid
    chpDropList.append(dfChpBig[dfChpBig['distance']>radius_shrink*radius].sample(frac=(rho_vapor/rho_liquid)))  # Outside droplet --> Vapor
    
    # Free memory
    del dfChpBig
    gc.collect()  # Garbage collection to free memory
    
    chpDrop = pd.concat(chpDropList)
    num_new = len(chpDrop)
    
    # Refresh particle ids
    chpDrop.reset_index(drop=True, inplace=True)
    chpDrop['pid'] = chpDrop.index + 1
    
    
    headerXML.find('headerinfo/length/x').text = str(xBoxLength)
    headerXML.find('headerinfo/length/y').text = str(yBoxLength)
    headerXML.find('headerinfo/length/z').text = str(zBoxLength)
    headerXML.find('headerinfo/number').text = str(num_new)
    headerXML.find('headerinfo/time').text = str(0.0)
    
    headerXMLTree.write(out_file_path_header)
    exp.exp_chp_bin_LD(out_file_path, chpDrop)
    
    return 1


def duplicate_droplet(in_file_path, out_file_path):
    '''
    Duplicate droplet and assign velocity

    :param str in_file_path: Path to input checkpoint (binary)
    :param str out_file_path: Path and filename of output checkpoint (binary)
    :return: 1 if successful
    '''
    
    
    in_file_path_header = in_file_path[:-4]+'.header.xml'
    out_file_path_header = out_file_path[:-4]+'.header.xml'
    
    # Read in checkpoint header
    headerXMLTree = et.parse(in_file_path_header)
    headerXML = headerXMLTree.getroot()
    
    num_old = int(headerXML.find('headerinfo/number').text)
    
    # Read in checkpoint data
    chp = imp.imp_chp_bin_DF(in_file_path)
    
    chpLeft  = chp.copy()
    chpRight = chp.copy()
    
    # Slightly modify velocities to prevent very identical speeds
    chpLeft['vx'] *= 0.999 + 0.002*np.random.rand(len(chpLeft))
    chpLeft['vz'] *= 0.999 + 0.002*np.random.rand(len(chpLeft))
    
    # Get box lengths from xml header file
    xBoxLength = float(headerXML.find('headerinfo/length/x').text)
    yBoxLength = float(headerXML.find('headerinfo/length/y').text)
    zBoxLength = float(headerXML.find('headerinfo/length/z').text)
    
    # Check if box is cubic
    if xBoxLength == yBoxLength == zBoxLength:
        pass
    else:
        print('Warning! Box is not cubic!')
    
    # Calculate real radius
    radius_real, profile = calc_radius(chp, [xBoxLength, yBoxLength, zBoxLength])
    
    print(f'   Real radius: {radius_real}')
    
    # Calculate distance of particle from center
    rx_dist = chp['rx'] - 0.5*xBoxLength
    ry_dist = chp['ry'] - 0.5*yBoxLength
    rz_dist = chp['rz'] - 0.5*zBoxLength
    chp['distance'] = (rx_dist*rx_dist + ry_dist*ry_dist + rz_dist*rz_dist)**0.5
    
    # Particles within droplet (radius plus approx. interface thickness (more or less arbitrarily set to 3.0))
    mask = chp['distance'] < (radius_real + 3.0)
    
    # Assign velocity
    chpLeft['vy']  = chpLeft['vy'].mask(mask, chpLeft['vy'] + 0.5*v_rel)
    chpRight['vy'] = chpRight['vy'].mask(mask, chpRight['vy'] - 0.5*v_rel)
    
    # Change componentid to make postprocessing easier
    chpLeft['cid']  = chpLeft['cid'].mask(mask, 2)
    chpRight['cid'] = chpRight['cid'].mask(mask, 3)
    
    # Cut both droplet checkpoints and keep a distance of 1 sigma between both
    sizeCutted = 0.5*yBoxLength + radius + 0.5*distance - 0.5  # Size in y direction of single box when cutted
    chpLeft = chpLeft[chpLeft['ry'] <= sizeCutted]
    chpRight = chpRight[chpRight['ry'] >= (yBoxLength - sizeCutted)]
    
    # Shift right checkpoint
    chpRight.loc[:,'ry'] = chpRight['ry'] + 2*radius + distance
    
    chpCollision= pd.concat([chpLeft,chpRight], axis=0, ignore_index=True)
    
    # Refresh particle ids
    chpCollision['pid'] = np.arange(len(chpCollision)) + 1
    
    
    headerXML.find('headerinfo/length/x').text = str(boxSize[0])
    headerXML.find('headerinfo/length/y').text = str(boxSize[1])
    headerXML.find('headerinfo/length/z').text = str(boxSize[2])
    headerXML.find('headerinfo/number').text = str(len(chpCollision))
    headerXML.find('headerinfo/time').text = str(0.0)
    
    headerXMLTree.write(out_file_path_header)
    exp.exp_chp_bin_DF(out_file_path, chpCollision)
    
    return  1




#%% Main program which can be called from the command line
if __name__ == '__main__':
    
    
    #%% Initial conditions
    
    radius = 30        # Approximate radius of droplet
    v_rel = 2.0        # Relative velocity of droplets; each droplet gets 0.5*v_rel assigned
    temperature = 0.7  # Initial temperature of system
    distance = 10      # Distance between droplets measured between surfaces
    
    
    # Other
    num_reps = 10
    timestep = 0.003647347
    
    radius_shrink = 1.00  # Cut out with larger radius to compensate shrinking due to surface tension
    
    # For reproducibility
    random.seed(10)
    
    #%% Computer settings
    
    # Path to folder where simulation files will be created
    work_folder = 'PATH'
    work_folder = os.path.join(work_folder, f'dropColl_r{radius}_v{str(v_rel).replace(".", "-")}_T{str(temperature).replace(".", "-")}')
    
    # Path to ls1 executable
    ls1_exec = 'PATH/build/src/MarDyn'
    
    # Path to templates
    template_folder = 'PATH/templates'
    
    
    #%% Derived values
    xLength = 10*radius
    boxSize = [xLength, xLength+2*radius+distance, xLength]  # Size of final simulation domain
    volume_droplet = (4/3)*np.pi*radius**3
    volume_liquid = 2*volume_droplet  # Liquid volume in final simulation domain
    volume_vapor = boxSize[0]*boxSize[1]*boxSize[2] - volume_liquid  # Vapor volume in final simulation domain
    
    rho_liquid,rho_vapor = rho_droplet_vrabec2006(temperature, radius)  # Get saturation properties at specified temperature

    numParts = int(rho_liquid*volume_liquid + rho_vapor*volume_vapor)  # Approximate number of particles in final simulation domain
    
    timesteps_collision = int(distance/(v_rel*timestep))  # Timesteps until droplets collide
    
    print(f'Work folder: {work_folder}')
    print(f'Setting up simulation with:')
    print(f'   Radius:              {radius}')
    print(f'   Relative velocity:   {v_rel}')
    print(f'   Temperature:         {temperature}')
    print(f'   At distance:         {distance}')
    print(f'   Number of particles: {numParts}')
    print(f'   Steps to collision:  {timesteps_collision}')

    
    
    # Copy xml files to work_folder
    shutil.copytree(os.path.join(template_folder, 'xml'), work_folder)
    
    
    #%% Step 1: Create and equilibrate bulk liquid
    print(f'Step 1:')
    # Create folder for current step
    folder_step1 = os.path.join(work_folder, 's1_Bulk-Liquid')
    os.mkdir(folder_step1)
    
    # Values to replace in template config xml
    replace_dict = {'ReplaceRedTemp': str(temperature), 'ReplaceBox': str(boxSize[0]/num_reps), 'ReplaceDensity': str(round(rho_liquid,8))}
    
    # Copy template file and replace placeholders
    with open(os.path.join(template_folder, 'tpl_1_generateLiquid.xml'), 'rt') as file_in:
        with open(os.path.join(folder_step1, 'config_step_1_generateLiquid.xml'), 'wt') as file_out:
            for line in file_in:
                for k in replace_dict.keys():
                    line = line.replace(k, replace_dict[k])
                file_out.write(line)
    
    
    # Run ls1 simulation
    print(f'   Running simulation')
    os.system(f'cd {folder_step1}; mpirun -np 8 {ls1_exec} config_step_1_generateLiquid.xml --final-checkpoint=0 > out.log')
    print(f'   Simulation done')
    
    
    #%% Step 2: Replicate bulk liquid and cut out spherical droplet
    print(f'Step 2:')
    # Create folder for current step
    folder_step2 = os.path.join(work_folder, 's2_Droplet-Equi')
    os.mkdir(folder_step2)
    
    # Values to replace in template config xml
    replace_dict = {'ReplaceRedTemp': str(temperature), 'ReplaceBox': str(boxSize[0]), 'ReplaceDensity': str(round(rho_liquid,8))}
    
    # Copy template file and replace placeholders
    with open(os.path.join(template_folder, 'tpl_2_dropletSingle.xml'), 'rt') as file_in:
        with open(os.path.join(folder_step2, 'config_step_2_dropletSingle.xml'), 'wt') as file_out:
            for line in file_in:
                for k in replace_dict.keys():
                    line = line.replace(k, replace_dict[k])
                file_out.write(line)
    
    
    # Replicate bulk liquid and cut out spherical droplet
    in_file_path = os.path.join(folder_step1, 'cp_binary_bulk-1.restart.dat')
    out_file_path = os.path.join(folder_step2, 'cp_binary_single-droplet.restart.dat')
    print(f'   Cutout droplet with radius: {radius_shrink*radius}')
    if chp_replicate_cut(in_file_path,out_file_path): print(f'   Replication and cutout successful')
    
    
    # Run ls1 simulation
    print(f'   Running simulation')
    os.system(f'cd {folder_step2}; mpirun -np 8 {ls1_exec} config_step_2_dropletSingle.xml --final-checkpoint=0 > out.log')
    print(f'   Simulation done')
    
    
    #%% Step 3: Duplicate droplet and set up initial checkpoint
    # Create folder for current step
    folder_step3 = os.path.join(work_folder, 's3_Droplet_Collision')
    os.mkdir(folder_step3)
    
    # Values to replace in template config xml
    replace_dict = {'ReplaceRedTemp': str(temperature), 'ReplaceBox': str(boxSize[0]), 'ReplaceYBox': str(boxSize[1]), 'ReplaceDensity': str(round(rho_liquid,8))}
    
    # Copy template file and replace placeholders
    with open(os.path.join(template_folder, 'tpl_3_dropletCollision.xml'), 'rt') as file_in:
        with open(os.path.join(folder_step3, 'config_step_3_dropletCollision.xml'), 'wt') as file_out:
            for line in file_in:
                for k in replace_dict.keys():
                    line = line.replace(k, replace_dict[k])
                file_out.write(line)
    
    
    # Duplicate droplet and set up initial checkpoint
    in_file_path = os.path.join(folder_step2, 'cp_binary_droplet_equi-1.restart.dat')
    out_file_path = os.path.join(folder_step3, 'cp_binary_dropletColl_init.restart.dat')
    if duplicate_droplet(in_file_path,out_file_path): print(f'   Duplication successful')
    
    
    
    #%% Read checkpoints and plot for check
    
    # Read checkpoint before and after equilibration
    radius_init, profile_init = calc_radius_chp(os.path.join(folder_step2, 'cp_binary_single-droplet.restart.dat'))
    
    radius_equi, profile_equi = calc_radius_chp(os.path.join(folder_step2, 'cp_binary_droplet_equi-1.restart.dat'))
    
    plt.figure()
    plt.plot(profile_init.index, profile_init['rho'], label='Cutout')
    plt.plot(profile_equi.index, profile_equi['rho'], label='Equilibrated')
    plt.legend()
    
    # Plot position of particles in final initial checkpoint
    chpCollision = imp.imp_chp_bin_DF(os.path.join(folder_step3, 'cp_binary_dropletColl_init.restart.dat'))

    fig, ax = plt.subplots(1, 1)
    plt.scatter(chpCollision['ry'],chpCollision['rx'],c=chpCollision['vy'])
    ax.set_aspect('equal', 'box')
    
    