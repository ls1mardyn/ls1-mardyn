import xml.etree.ElementTree as et

import ppls1.imp.chp as imp
import ppls1.exp.chp as exp


#%% Function to generate a binary checkpoint which is num_reps times bigger in x and z direction than the input
def chp_replicate(num_reps, in_file_path, out_file_path, createHeader=True):
    '''
    Function to generate a binary checkpoint including header which is num_reps times bigger in x and z direction than the input

    :param int num_reps: Number of replications in each direction
    :param str in_file_path: Path to input checkpoint (binary)
    :param str out_file_path: Path and filename of output checkpoint (binary)
    :param bool createHeader: Specify if header file should be created automatically
    :return: 1 if successful
    '''
    in_file_path_header = in_file_path[:-4]+'.header.xml'
    out_file_path_header = out_file_path[:-4]+'.header.xml'
    
    headerXMLTree = et.parse(in_file_path_header)
    headerXML = headerXMLTree.getroot()
    
    num_old = int(headerXML.find('headerinfo/number').text)
    num_new = num_reps*num_reps*num_old
    
    # Read in checkpoint
    chp = imp.imp_chp_bin_LD(in_file_path)
    chpNew = []
    
    # Get box lengths from xml header file
    xBoxLength = float(headerXML.find('headerinfo/length/x').text)
    zBoxLength = float(headerXML.find('headerinfo/length/z').text)
    
    # replicate particles
    for par in chp:
        for xx in range(num_reps):
            for zz in range(num_reps):
                tmpPar = par.copy()
                tmpPar['rx'] = xx*xBoxLength + par['rx']
                tmpPar['rz'] = zz*zBoxLength + par['rz']
                chpNew.append(tmpPar)
                
    # refresh particle ids
    for pi in range(num_new):
        chpNew[pi]['pid']=pi+1
                
    
    headerXML.find('headerinfo/length/x').text = str(num_reps*xBoxLength)
    headerXML.find('headerinfo/length/z').text = str(num_reps*zBoxLength)
    headerXML.find('headerinfo/number').text = str(num_new)
    headerXML.find('headerinfo/time').text = str(0.0)
    
    if createHeader: headerXMLTree.write(out_file_path_header)
    exp.exp_chp_bin_LD(out_file_path,chpNew)
    
    return 1


#%% Main program which can be called from the command line
if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('num_reps', help='Number of replications in each direction')
    parser.add_argument('in_file_path', help='Path to input checkpoint')
    parser.add_argument('out_file_path', help='Path and filename of output checkpoint')
    args = parser.parse_args()
    num_reps = int(float(args.num_reps))
    in_file_path = args.in_file_path
    out_file_path = args.out_file_path
    
    if chp_replicate(num_reps,in_file_path,out_file_path, createHeader=True): print('Replication ('+str(num_reps)+'x) successful')
