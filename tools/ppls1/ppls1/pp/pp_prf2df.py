import glob
import json
import os

import numpy as np
import pandas as pd


#%% Function to generate dataframe out of simulation profile data files and, if wanted, export to path
def prf2df(folderSims,flgExport=False,outFile='filenameForExport',configFile='config*.xml',quietMode=False):
    '''
    Function to generate dataframe out of simulation profile data files and, if wanted, export to path

    :param str folderSims: Path to the simulation folder
    :param bool flgExport: Flag if dataframe is exported as json in COLUMNS format (default: 0)
    :param str outFile: Path and name of the file to be exported (default: folderSims)
    :param str configFile: Path and name of the config.xml file (default: config*.xml)
    :param str quietMode: Do not print imported files (default: 0)
    :return: dataframe including the simulation data
    '''
    
    flgAsOneDF = 1
    
    ## Overwrite default value, if not set explicitly; not applicable, if called by __main__
    if flgExport == True and outFile == 'filenameForExport':
        outFile = folderSims+'/dataframe_simulation.json'
    
    ## Look in path 'folderSims' for desired files (profile data)
    flist_scal_all = sorted(glob.glob(folderSims+'/scal*all*'))
    flist_vect_all = sorted(glob.glob(folderSims+'/vect*all*'))
    flist_scal_pos = sorted(glob.glob(folderSims+'/scal*pos*'))
    flist_vect_pos = sorted(glob.glob(folderSims+'/vect*pos*'))
    flist_scal_neg = sorted(glob.glob(folderSims+'/scal*neg*'))
    flist_vect_neg = sorted(glob.glob(folderSims+'/vect*neg*'))
    
    ## Import data; set method, how data is read and processed
    flgMtdReadPrf = 2 ## Method 2 is recommended
    
    if flgMtdReadPrf == 1: # Legacy method
        import ppls1.imp.prf as imp_prf
        import ppls1.pp.prf as pp_prf
        
        fplists = {'scal_all': flist_scal_all, 'scal_pos': flist_scal_pos, 'scal_neg': flist_scal_neg, 'vect_all': flist_vect_all, 'vect_pos': flist_vect_pos, 'vect_neg': flist_vect_neg}
        
        data_raw2 = {key:imp_prf.imp_prf(val) for (key,val) in fplists.items()}
        for (_,val2) in data_raw2.items():
            pp_prf.split_comp_prf(val2)
        
        ## Format as Pandas Dataframe
        cidList = range(val2[0]['numC']+1)
        
        dfProfTemp = pd.DataFrame(columns=[int(s['timestep']) for s in val2], index=cidList)
        dfTC = pd.DataFrame()
        
        for cid in cidList:
            for timestep in list(dfProfTemp.columns):
                for d in data_raw2.items():
                    for timeData in d[1]:
                        if timeData['timestep']==timestep:
                            for dataKey in timeData['header']:
                                dfTC[dataKey+d[0][-4:]] = timeData[dataKey][cid]
                
                ## Delete old positions and add right one
                dfTC = dfTC.drop(['pos_all'], axis=1)
                dfTC = dfTC.drop(['pos_pos'], axis=1)
                dfTC = dfTC.drop(['pos_neg'], axis=1)
                dfTC.insert(loc=0, column='pos', value=timeData['pos'])
                
                dfProfTemp[timestep][cid] = dfTC
                dfTC = pd.DataFrame()
    
    if flgMtdReadPrf == 2:
        fplists = flist_scal_all + flist_scal_pos + flist_scal_neg + flist_vect_all + flist_vect_pos + flist_vect_neg
        
        dfTemp = pd.read_fwf(fplists[0])
        
        timesteps = list()
        
        for fp in fplists:
            timestep = int(fp.rpartition('TS')[-1].split('.')[0])
            if timestep not in timesteps:
                timesteps.append(timestep)
            
        components = np.arange(int((len(dfTemp.columns)-1)/13))
        
        dfProfTemp = pd.DataFrame(columns=timesteps, index=components)
        dfTC = pd.DataFrame()
        
        for fp in fplists:
            if not quietMode: print(fp)
            direction = 'NaN'
            direction = fp.rpartition('quant_')[-1].split('_reg')[0]
            dfTC = pd.read_fwf(fp)
            timestep = int(fp.rpartition('TS')[-1].split('.')[0])
            
            ## Iterate over columns of current dataframe
            for col in dfTC.columns:
                dfTemp = pd.DataFrame()
                if col == 'pos':
                    continue
                else:
                    cid = int(col.rpartition('[')[-1].split(']')[0])
                    if (not isinstance(dfProfTemp[timestep][cid], pd.DataFrame)) or (not 'pos' in dfProfTemp[timestep][cid]):
                        dfTemp['pos'] = dfTC['pos']
                    dfTemp[col[:-3]+'_'+direction] = dfTC[col]
                if isinstance(dfProfTemp[timestep][cid], pd.DataFrame):
                    dfProfTemp[timestep][cid] = pd.concat([dfProfTemp[timestep][cid], dfTemp], axis=1)
                    #dfProfTemp[timestep][cid] = dfProfTemp[timestep][cid].loc[:,~dfProfTemp[timestep][cid].columns.duplicated()]
                else:
                    dfProfTemp[timestep][cid] = dfTemp
                    
    if flgMtdReadPrf == 3: # NOT WORKING YET
        fplists = flist_scal_all + flist_scal_pos + flist_scal_neg + flist_vect_all + flist_vect_pos + flist_vect_neg
        
        dfTemp = pd.read_fwf(fplists[0])
        
        timesteps = list()
        
        for fp in fplists:
            timestep = int(fp.rpartition('TS')[-1].split('.')[0])
            if timestep not in timesteps:
                timesteps.append(timestep)
            
        components = np.arange(int((len(dfTemp.columns)-1)/13))
        
        dfTC = pd.DataFrame()
        dfProfTemp = pd.DataFrame()
        
        for fp in fplists:
            if not quietMode: print(fp)
            direction = 'NaN'
            direction = fp.rpartition('quant_')[-1].split('_reg')[0]
            timestep = int(fp.rpartition('TS')[-1].split('.')[0])
            dfTC = pd.read_fwf(fp)
            
            for comp in components:
                dfComp = dfTC.loc[:,dfTC.columns.str.contains('['+str(comp)+']')].copy()
                dfComp.columns = dfComp.columns.str.rstrip('['+str(comp)+']')
                dfComp.columns = dfComp.columns + '_' + direction
                dfComp['cid'] = comp
                dfComp['timestep'] = timestep
                idxStart = int( ((timestep/timesteps[-1]*len(timesteps)-1)*len(components)*len(dfComp))
                                 +(comp/components[-1]*(len(components)-1)*len(dfComp)) )
                idxEnde = idxStart + len(dfComp)
                dfComp.index = range(idxStart, idxEnde)
                dfProfTemp = dfProfTemp.merge(dfComp, left_index=True, right_index=True, how="right", sort=False)
                dfProfTemp = pd.concat([dfProfTemp, dfComp], axis=1)
                
        flgAsOneDF = False # Since already as one data frame using method 3
        
    
    if flgAsOneDF:
        for col in dfProfTemp.columns:
            for ind in dfProfTemp.index:
                dfProfTemp[col][ind]['timestep'] = col
                dfProfTemp[col][ind]['cid'] = ind
                # print(str(col)+" "+str(ind))
                # print(dfProfTemp[col][ind])

        # print('make frames')
        frames = list()

        for col in dfProfTemp.columns:
            # print(col)
            for ind in dfProfTemp.index:
                # print(ind)
                frames.append(dfProfTemp[col][ind])
                # print(frames)

        #print(frames)
        dfProfData = pd.concat(frames, ignore_index=True)
    else:
        dfProfData = dfProfTemp
    
    metaData = dict()
    metaData['dirPath'] = folderSims
    
    listConfigFiles = sorted(glob.glob(folderSims+'/'+configFile))
    if len(listConfigFiles) == 0:
        print('Warning: No config*.xml file found!')
    else:
        if len(listConfigFiles) > 1:
            print('Warning: More than one config*.xml file found! Choosing first one...')
        with open(listConfigFiles[0],'r') as datei:
            metaData['xmlConfig'] = datei.read()
    
    listOutputFiles = sorted(glob.glob(folderSims+'/*out*'))
    if len(listOutputFiles) == 0:
        print('Warning: No output file found!')
    else:
        if len(listOutputFiles) > 1:
            print('Warning: More than one output file found!')
        output = list()
        for entry in listOutputFiles:
            with open(entry,'r') as outputFile:
                outputAll = outputFile.readlines()
                output.append(''.join(outputAll[:500] + outputAll[-500:]))
        metaData['outputLog'] = output
    
    ## Export to JSON in 'columns' format; Only for all in one DF
    #flgExport = True ## Set by function input
    if flgExport == True:
        df2json(dfProfData, metaData, outFile)
    
    return dfProfData, metaData

def df2json(dfProfData, metaData, outFile):
    '''
    Function to export dataframe in json format

    :param df dfProfData: Dataframe to be exported
    :param str outFile: Path and name of the file to be exported into
    :return: -
    '''
    dfData_json = dfProfData.to_json(orient='columns',double_precision=15)
    dfMeta_json = json.dumps(metaData)
    
    json_export={
        'metadata':'aqdMETAkgsm' ,
        'profile_data':'aqdjnpjkgsm'  ## Random string to be replaced after json.dumps
    }

    top_level_json=json.dumps(json_export).replace('\"aqdMETAkgsm\"',dfMeta_json)
    top_level_json=top_level_json.replace('\"aqdjnpjkgsm\"',dfData_json)

    with open(outFile, 'w') as file:
        file.write(top_level_json+'\n')

    print('Exported successfully to '+str(outFile))


#%% Function to get dataframe back from json (import in Python)
def json2df(infile):
    '''
    Function to get dataframe back from json

    :param str infile: Path and name of the file to be read
    :return: dataframe including the simulation data
    '''
    print('Importing json file '+infile)
    with open(infile,'r') as file:
        data_json = json.load(file)
    
    ## OLD METHOD
    # with open(os.path.dirname(infile)+'/aqdMETAkgsm.tmp.json', 'w') as file:
    #     json.dump(data_json['metadata'], file)
    # metaDataImp = pd.read_json(os.path.dirname(infile)+'/aqdMETAkgsm.tmp.json',orient='columns')
    # os.remove(os.path.dirname(infile)+'/aqdMETAkgsm.tmp.json')
    #
    # with open(os.path.dirname(infile)+'/aqdjnpjkgsm.tmp.json', 'w') as file:
    #     json.dump(data_json['profiles_dataframe'], file)   
    # dfDataImp = pd.read_json(os.path.dirname(infile)+'/aqdjnpjkgsm.tmp.json',orient='columns')
    # dfDataImp = dfDataImp.sort_index() ## Sort indices
    # os.remove(os.path.dirname(infile)+'/aqdjnpjkgsm.tmp.json')
    
    metaDataImp = data_json['metadata']
    dfDataImp = pd.DataFrame(eval(str(data_json['profile_data'])))
    #dfDataImp['iNr'] = [float(i) for i in dfDataImp.index]
    #dfDataImp.sort_values(by=['iNr'],inplace=True)
    #dfDataImp.reset_index
    #dfDataImp.drop('iNr', 1, inplace = True)
    #dfDataImp.index.rename('pos')
    dfDataImp.index = dfDataImp.index.rename('pos')
    dfDataImp.index = dfDataImp.index.astype('float')
    dfDataImp = dfDataImp.sort_index(ascending=True)
    
    #dfDataImp = dfDataImp.sort_index() ## Sort indices
    
    
    ## Only if dataframe is nested
    if len(dfDataImp.index) < 100:
        print('json2df: Nested data frame')
        dfDataImp = dfDataImp.reindex(sorted(dfDataImp.columns), axis=1) ## Sort columns by name
        for cols in dfDataImp.columns:
            for inds in dfDataImp.index:
                dfDataImpCI = pd.DataFrame(dfDataImp[cols][inds]).T
                dfDataImp[cols][inds] = dfDataImpCI
    
    return dfDataImp, metaDataImp


#%% Main program which can be called from the command line; Will convert profile data to dataframe including export
if __name__ == '__main__':
    
    flgMtdMain = 2 ## Decide on input: 1 for command line argument; 2 for hard coded values
    
    if flgMtdMain == 1:
        import argparse
        argparser = argparse.ArgumentParser(description='Simulation profile data to dataframe, optional saving')
        argparser.add_argument('-s','--simPath', required=True, type=str, help='Path of simulation data')
        argparser.add_argument('-e','--exportFile', required=True, type=str, help='File and path to save dataframe into')
        
        args = vars(argparser.parse_args())
        simPath = args['simPath']
        expFile = args['exportFile']
        
        if expFile is None:
            flgExport = False
            print('Dataframe will not be saved!')
        else:
            flgExport = True
            print('Saving to: '+expFile)
    elif flgMtdMain == 2:
        print('Input set manually in script')
        simPath = '/home/pcfsuser/Simon/ls1/projects/heatflux/untersuchungDC_TT/parVar/supercritical/noDC_noTT_T1.98'
        flgExport = False
        expFile = simPath+'/dataframe_simulation.json'
    else:
        print('Requested input method not available')
    
    print('Path to simulation: '+simPath)
    
    dfDataProf, metaData = prf2df(simPath,flgExport,expFile)
    #dfDataImp, metaDataImp = json2df(expFile)
    
    ## Check equality of dataframes
    # import math
    # for col in dfDataProf.columns:
    #     for ind in dfDataProf.index:
    #         if math.isnan(dfDataImp[col][ind]): continue
    #         if round(dfDataImp[col][ind],10) != round(dfDataProf[col][ind],10):
    #             print(str(col)+" "+str(ind)+" "+str(dfDataImp[col][ind])+" "+str(dfDataProf[col][ind])+" ")
    
    print('Done')
    
