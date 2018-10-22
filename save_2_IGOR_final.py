import numpy as np
import scipy as sp
import math as ma
import inspect

def is_array(the_array) :
    return isinstance(the_array,((list,tuple)))
    
def retrieve_name(varName) :
    for fi in reversed(inspect.stack()) :
        names = [var_name for var_name, var_val in fi.frame.f_locals.items() if var_val is varName]
        if len(names) > 0 :
            return names[0]

def save_single_data(folderPath,DataWave,waveName,run,start_at,end_at) :
    waveName = waveName.replace('=','_')
    delta = (end_at - start_at)/DataWave.shape[0]
    with open(folderPath + 'data_output_2_igor_run%04d' % run + '.itx', 'w') as out_file : 
        for i in range(DataWave.shape[0]+1):
            out_string = ''
            if i == 0:
                out_string = 'IGOR\nWAVES/D/O/N=('+str(DataWave.shape[0])+') '+ waveName +'\nBEGIN\n'
                out_string += str(DataWave[i]) + '\n'
            elif i < len(DataWave):
                out_string += str(DataWave[i]) + '\n'
            else:
                out_string += 'END\n\nX SetScale/P y 0,1, \"Counts [a.u.]\", ' + waveName
                out_string += '\nX SetScale/P x {0:g},{1:g}, \"Kinetic Energy [eV]\",  '.format(start_at, delta)  + waveName
                out_string += '\nX Display ' + waveName
                out_string += '\nX Legend/RT/N=text0/F=0/B=1/A=LT\nX ModifyGraph nticks=10,minor=1\n'
                
            out_file.write(out_string)
            
def saveData(folderPath,DataWave,DataXWave,waveName,run,runType) :
    if is_array(DataWave) == True :
        DataLength = len(DataWave)
    else :
        DataLength = DataWave.shape[0]
    XfileName = retrieve_name(DataXWave)
    fileName = retrieve_name(DataWave)
    
    start_value = DataXWave[DataLength-1]
    end_value = DataXWave[0]
    
    if XfileName == 'Time' :
        xWave = 'timeWave_for_'
    else :
        xWave = 'EnergyWave_for_'
        
    SaveName = 'data2igor_' + fileName + '_t0_' + runType + '_%03d'
    
    with open(folderPath + SaveName % run + '.itx', 'w') as out_file :
        for i in range(DataLength+1):
            out_string = ''
            if i == 0:
                out_string = 'IGOR\nWAVES/D/O/N=(' + str(DataLength) + ')  '+ xWave + runType + '_%03d' % run 
                out_string += ', ' + waveName
                out_string += '\nBEGIN\n'
                out_string += str(DataXWave[i]) + '\t' + str(DataWave[i]) + '\n'
            elif i < DataLength:
                out_string += str(DataXWave[i]) + '\t' + str(DataWave[i]) + '\n'
            else:
                out_string += 'END\n'
                if XfileName == 'Time' :
                    out_string += '\nX SetScale d 0,0, \"Normalized Intensity\" , ' + waveName
                    out_string += '\nX SetScale d 0,0, \"Time [ps]\", ' + xWave + runType + '_%03d' % run
                else : 
                    out_string += '\nX SetScale d 0,0, \"Counts [a.u.]\" , ' + waveName 
                    out_string += '\nX SetScale d 0,0, \"Kinetic Energy [eV]\", ' + xWave + runType + '_%03d' % run
                out_string += '\nX Display ' + waveName + ' vs EnergyWave_for_' + runType + '_%03d' % run
                out_string += '\nX Legend/RT/N=text0/F=0/B=1/A=LT'
                out_string += '\nX ModifyGraph nticks=10,minor=1,lblMargin(left)=10' + '\n'
                if XfileName == 'Time' :
                    out_string += 'SetAxis left ' + str(start_value) + ', ' + str(end_value) + ';DelayUpdate'
                          
            out_file.write(out_string)

def saveDataForTimeDelay(folderPath,DataWave,DataXWave,waveName,run,runType) :
    if waveName.find('NP',0) == 0 :
        SaveName = 'data2igor_from_NP_' + runType + '_%03d'
        waveName = waveName + '_%03d' % run
        D_color = 1
        append = 1
    else :
        SaveName = 'data2igor_from_' + runType + '_%03d'
        waveName = waveName + '_%03d' % run
        D_color = 0
        append = 0
    #waveName = waveName.replace('=','_')
    if is_array(DataWave) == True :
        DataRow = len(DataWave)
        DataColum = len(DataWave[0])
    else :
        DataRow = DataWave.shape[0]
        DataColum = DataWave.shape[1]
        
    with open(folderPath + SaveName % run + '.itx', 'w') as out_file :
        for i in range(DataColum+1):
            out_string = ''
            if i == 0:
                out_string = 'IGOR\nWAVES/D/O/N=(' + str(DataColum) + ')    Energy_wave_for_' + runType + '_%03d' % run 
                for j in range(DataRow):
                    out_string += ', ' + waveName + '_wave%03d' % (j+1)
                out_string += '\nBEGIN\n'
                out_string += str(DataXWave[i])
                for j in range(DataRow):
                    out_string += '\t' + str(DataWave[j,i])
                out_string += '\n'
            elif i < DataColum:
                out_string += str(DataXWave[i])
                for j in range(DataRow):
                    out_string +=  '\t' + str(DataWave[j,i])
                out_string += '\n'
            else:
                out_string += 'END\n'
                for j in range(DataWave.shape[0]) :
                    out_string += '\nX SetScale d 0,0, \"Intensity in Arb. units\" , ' + waveName + '_wave%03d' % (j+1)
                    out_string += '\nX SetScale d 0,0, \"Kinetic Energy [eV]\", Energy_wave_for_' + runType + '_%03d' % run  # \"Intensity in arb. units\" or  \"Counts [a.u.]\"
                    if j == 0 and append == 0 :
                        out_string += '\nX Display ' + waveName + '_wave%03d' % (j+1) + ' vs Energy_wave_for_' + runType + '_%03d' % run
                    else :
                        out_string += '\nX AppendToGraph ' + waveName + '_wave%03d' % (j+1) + ' vs Energy_wave_for_' + runType + '_%03d' % run
                    if j == 0 :
                        if append == 0 :
                            out_string += '\nX Legend/RT/N=text%d/F=0/B=1/A=LT' %j
                    out_string += '\nX ModifyGraph grid(bottom)=2,nticks(left)=0,nticks(bottom)=10,minor=1,lblMargin(left)=10'
                    if D_color == 1 :
                        out_string += '\nX ModifyGraph rgb('+ waveName + '_wave%03d' % (j+1) +')=(0,0,0)'
                    out_string += '\nX ModifyGraph offset(' + waveName + '_wave%03d' % (j+1) + ')={0,%4f}' %(j*0.002) + '\n'
                                
            out_file.write(out_string)

def save_2D_data(folderPath,DataWave,waveName,start_at,end_at,run,runType) :
    if is_array(DataWave) == True :
        DataRow = len(DataWave)
        DataColum = len(DataWave[0])
    else :
        DataRow = DataWave.shape[0]
        DataColum = DataWave.shape[1]
        
    fileName = retrieve_name(DataWave)
    
    if fileName == 'Total_profile' :
        output_Name = 'output_2D_Igor_'
    else :
        output_Name = 'output_2D_Igor_Bin_'
    
    if end_at != 1 :
        delta = (start_at-end_at)/DataColum
    if end_at != 1 :
        print (start_at,end_at,delta)
    else :
        print (start_at,end_at)
    with open(folderPath + output_Name + runType + '_data_run%03d' % run + '.itx', 'w') as out_file :
        for i in range(DataRow+1):
            out_string = ''
            if i == 0 :
                if fileName == 'Total_profile' :
                    out_string = 'IGOR\nWAVES/D/O/N=('+str(DataRow)+','+str(DataColum)+') Intensity_2D_' + runType + '%03d' % run + '\nBEGIN\n'
                else :
                    out_string = 'IGOR\nWAVES/D/O/N=('+str(DataRow)+','+str(DataColum)+') Intensity_2D_Bin_' + runType + '%03d' % run + '\nBEGIN\n'
                for j in range(DataColum):
                    out_string += str(DataWave[i][j]) + '\t'
                
                out_string += '\n'
            elif i < DataRow:
                for j in range(DataColum):
                    out_string += str(DataWave[i][j]) + '\t'
                
                out_string += '\n'
            else:
                if fileName == 'Total_profile' :
                    out_string += 'END\n\nX SetScale/P y 0,1, \"Pixels\", Intensity_2D_' + runType + '%03d' % run   #Region Iteration [a.u.]
                    out_string += '\nX SetScale/P x {0:f},{1:f}, \"Shots\",  '.format(start_at, end_at) + 'Intensity_2D_' + runType + '%03d' % run #Kinetic Energy [eV]
                    out_string += '\nX NewImage/K=0 Intensity_2D_' + runType + '%03d' % run
                    out_string += '\nX ModifyImage Intensity_2D_' + runType + '%03d ctab= {*,*,PlanetEarth,0}' % run
                else :
                    out_string += 'END\n\nX SetScale/P x 0,1, \"Pixels\", Intensity_2D_Bin_' + runType + '%03d' % run
                    out_string += '\nX SetScale/P y {0:f},{1:f}, \"Time delay stage position [ps]\",  '.format(end_at, delta) + 'Intensity_2D_Bin_' + runType + '%03d' % run
                    out_string += '\nX NewImage/K=0 Intensity_2D_Bin_' + runType + '%03d' % run
                    out_string += '\nX ModifyImage Intensity_2D_Bin_' + runType + '%03d ctab= {*,*,PlanetEarth,0}' % run
                    out_string += '\nX SetAxis/R left {0:f},{1:f}'.format(start_at,end_at)
                out_string += '\nX ModifyGraph nticks(top)=20,nticks(left)=10,minor=1,sep(top)=5,sep(left)=5'
                out_string += '\nX ModifyGraph margin(left)=40,margin(right)=120,margin(top)=40'
                out_string += '\nX ModifyGraph lblMargin=5,lblLatPos=1,fSize=14'
                if fileName == 'Total_profile' :
                    out_string += '\nX ColorScale/C/N=text0/F=0/A=RT/X=1.50/Y=4.50/E=2 heightPct=100,frame=2.00,image=Intensity_2D_' + runType + '%03d,lblMargin=15,fsize=12,nticks=10,minor=1;DelayUpdate' % run
                else :
                    out_string += '\nX ColorScale/C/N=text0/F=0/A=RT/X=1.50/Y=4.50/E=2 heightPct=100,frame=2.00,image=Intensity_2D_Bin_' + runType + '%03d,lblMargin=15,fsize=12,nticks=10,minor=1;DelayUpdate' % run
                out_string += '\nX ColorScale/C/N=text0 \"\\\\Z20 Intensity\"\n'                    
            out_file.write(out_string)
