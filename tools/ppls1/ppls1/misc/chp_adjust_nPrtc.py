import ppls1.imp.chp as imp
import ppls1.exp.chp as exp
from random import shuffle


#%% Function to reduce number of particles in checkpoint
def chp_adjust_nPrtc(num_target,in_file_path,out_file_path):
    '''
    Function to generate dataframe out of simulation profile data files and, if wanted, export to path

    :param int num_target: Target number of particles
    :param str in_file_path: Path and name of the binary checkpoint to be read
    :param str out_file_path: Path and name of the binary checkpoint to be written
    '''
    
    chp=imp.imp_chp_bin_LD(in_file_path)
    num_particles_old = len(chp)
    
    # shuffle particle list
    shuffle(chp)
    
    # delete particles to match desired particle count
    chp=chp[:num_target]
    
    # refresh particle ids
    for pi in range(num_target):
        chp[pi]['pid']=pi+1
    
    # export checkpoint 
    exp.exp_chp_bin_LD(out_file_path,chp)
    
    print("   Number of particles successfully reduced ( "+str(num_particles_old)+" -> "+str(num_target)+" )")


#%% Main program which can be called from the command line
if __name__ == '__main__':
    
    flgInp=1 # 1: Parse args; 2: Hardcoded
    
    if flgInp==1:
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument('num_target', help='Target number of particles')
        parser.add_argument('in_file_path', help='Path to file to be changed')
        parser.add_argument('out_file_path', help='Path to file to be changed')
        args = parser.parse_args()
        num_target=int(float(args.num_target))
        in_file_path=args.in_file_path
        out_file_path=args.out_file_path
    
    if flgInp==2:
        num_target = 700
        work_folder = 'PATHTOSIMULATION'
        in_file_path = work_folder+'cp_binary-0.restart.dat'
        out_file_path = work_folder+'cp_binary_adjusted.restart.dat'
    
    chp_adjust_nPrtc(num_target,in_file_path,out_file_path)
