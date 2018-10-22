import numpy as np
import random
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import h5py
import sys
from os import walk
import os.path
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(0, '/Users/gladh/Documents/Projects/CO oxidation/src/')
import fit_data_tools
from fit_data_tools import *
import save_2_IGOR_final as sif

files = [5] # file number, we store the data into Time_zero folder
thr_camera = 20
Int_time = float (1) # interval delay gap, unit = ps
save2igor = 1

region_1s = 385
region_1e = 405
region_2s = 405
region_2e = 453
region_3s = 470
region_3e = 488
region_4s = 490
region_4e = 530

bin_regions = [region_1s,region_1e,region_2s,region_2e,region_3s,region_3e,region_4s,region_4e]

######################
# Load the data file #
######################
data_path = '/Volumes/Experiment_Data/DESY/PG2/BT2/'+'Time_zero%03d/' %files[0]
if not os.path.exists(data_path) :
     raise NameError('Check data_path! It does not exist!')
data_names = discover_datas(data_path)
print (data_names)

#############################################
# Creating empty matirx for collecting data #
#############################################
Total_profile = []
Total_V = []
Total_delay = []
Binned_profile = []
Binned_delay = []

################
# Data process #
################
for i in range(len(data_names)):
    h5file = h5py.File(data_path+data_names[i], 'r')
    Image = h5file['/FL1/Experiment/Camera/PG/Basler FW 1/image/'].value
    Voltage_num = h5file['/FL1/Experiment/PG/SIS8300 100MHz ADC/CH6/TD/'].value
    #Delay_stage0 = h5file['/FL1/Experiment/Pump probe laser/laser delay'].value
    Delay_stage = h5file['/FL1/Experiment/Pump probe laser/delay line IK220.0/ENC'].value
    
    
    Image[Image<thr_camera] = 0
    for m in np.arange((Image.shape[0])):
        Ind_profile = np.sum(Image[m,:,:],axis = 1)
        Total_profile.append(Ind_profile)
    
    Ind_voltage_num = np.average(Voltage_num, axis = 1)
    Ind_voltage_num = Ind_voltage_num*650/(2**16)
    Total_V.extend(Ind_voltage_num)
            
    Total_delay.extend(Delay_stage[:,1])

    
Total_profile = np.array(Total_profile) 
Total_V = np.array(Total_V)
Total_V = np.nan_to_num(Total_V)
Total_delay = np.array(Total_delay)
Total_delay = np.nan_to_num(Total_delay)

Total_profile = Total_profile[np.all([Total_V>10, Total_delay!=0], axis = 0)]
Total_delay = Total_delay[np.all([Total_V>10, Total_delay!=0], axis = 0)]

Delay_start = float(Total_delay.min())
Delay_end = float(Total_delay.max())
Binned_delay.append(Delay_start)
print (Delay_start)
print(Delay_end)
############################################
# Plot out the 2D map data wothout binning #
############################################

plt.figure(figsize=[14,6])   
ax = plt.subplot(1,1,1)
plt.imshow(Total_profile,vmin=0,vmax=200, origin='lower', interpolation='bicubic',cmap='hot',alpha = 0.9)
plt.xlabel('pixels',size=12)
plt.ylabel('shots',size=12)
plt.title('raw data')
plt.xlim(0,640)
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="4%", pad=0.1)
#plt.colorbar(cax=cax)
ax.set_aspect(aspect=0.01)
plt.show()

Spectrum = np.sum(Total_profile,axis = 0)
E_scale = np.arange(Spectrum.shape[0])
plt.figure(figsize=[9,6])
plt.plot(E_scale,Spectrum,ls = '-',linewidth ='3', c='black', alpha = 0.9)
plt.axvline(x=region_1s,ls='--',c='b')
plt.axvline(x=region_1e,ls='--',c='b')
plt.axvline(x=region_2s,ls='--',c='r')
plt.axvline(x=region_2e,ls='--',c='r')
plt.axvline(x=region_3s,ls='--',c='g')
plt.axvline(x=region_3e,ls='--',c='g')
plt.axvline(x=region_4s,ls='--',c='y')
plt.axvline(x=region_4e,ls='--',c='y')
plt.xlabel('pixels',size=12)
plt.ylabel('Summed Intensity',size=12)
plt.title('Summed raw data')
plt.xlim(0,640)
plt.show()

###########################################################
# Start to bin the data into desired delay time intervals #
###########################################################
Sum_0 = []
Sum_1 = []
Sum_2 = []
Sum_3 = []
d = 0
while (Delay_start + d*Int_time) < (Delay_end):
    mask = np.all([Total_delay >= (Delay_start + d*Int_time), Total_delay < (Delay_start + (d+1)*Int_time)], axis = 0)
    Int_Binned_profile = np.average(Total_profile[mask], axis = 0)
    Int_Binned_profile /= float(np.sum(Int_Binned_profile))
    Binned_profile.append(Int_Binned_profile)
    Binned_delay_mean = float(np.mean(Total_delay[mask]))
    Binned_delay.append(Binned_delay_mean)
    Sum_region_profile0 =float(np.sum(Int_Binned_profile[385:405]))
    Sum_region_profile1 =float(np.sum(Int_Binned_profile[405:453]))
    Sum_region_profile2 =float(np.sum(Int_Binned_profile[470:488]))
    Sum_region_profile3 =float(np.sum(Int_Binned_profile[490:530]))
    Sum_0.append(Sum_region_profile0)
    Sum_1.append(Sum_region_profile1)
    Sum_2.append(Sum_region_profile2)
    Sum_3.append(Sum_region_profile3)
    d += 1

Binned_profile = np.array(Binned_profile)
Binned_delay = np.array(Binned_delay)

##########################################
# Plot out the 2D map data after binning #
##########################################
a = np.arange(0,639.1,1)
b = Binned_delay
[x,y] = np.meshgrid(a,b)
z = Binned_profile
plt.figure(figsize=[12,6])   
ax = plt.subplot(1,1,1)
plt.pcolor(x,y,z, vmax = z.max(), cmap = 'hot')
plt.xlabel('Pixels',size=14)
plt.ylabel('Delay stage position (ps)',size=14)
plt.title('Binned Data', size = 16)
plt.colorbar()
ax.set_aspect(aspect=8)
plt.show()

Spectrum_1 = Binned_profile[44,:]
Spectrum_2 = Binned_profile[50,:]
Spectrum_3 = Binned_profile[47,:]
E_scale_bin = np.arange(Spectrum_1.shape[0])
plt.figure(figsize=[9,6])
plt.plot(E_scale_bin,Spectrum_1,ls = '-',linewidth ='1', c='b', alpha = 0.9, label = 'Before the pump')
plt.plot(E_scale_bin,Spectrum_2,ls = '-',linewidth ='1', c='r', alpha = 0.9, label = 'After the pump')
plt.plot(E_scale_bin,Spectrum_3,ls = '-',linewidth ='1', c='g', alpha = 0.9, label = 'On the pump')
plt.xlabel('Pixels',size=12)
plt.ylabel('Normalized Intensity',size=12)
plt.xlim(0,640)
plt.legend(frameon=False,loc='best', prop={'size': 12})
plt.tight_layout()
plt.show()

Time = Binned_delay[1:]
plt.figure(figsize=[8,6])
plt.plot(Time,Sum_0,'o-',linewidth ='2', c='b', label = '385-405')
plt.plot(Time,Sum_1,'o-',linewidth ='2', c='red', label = '405-453')
plt.plot(Time,Sum_2,'o-',linewidth ='2', c='g', label = '470-488')
plt.plot(Time,Sum_3,'o-',linewidth ='2', c='y', label = '490-530')
plt.xlim(Time.max(), Time.min())
plt.xlabel('Time (ps)',size=14)
plt.ylabel('Normalized Intensity',size=14)
plt.title('Selected region vs. Time', size = 16)
plt.grid(True)
plt.legend(frameon=False,loc='best', prop={'size': 12})
plt.tight_layout()
plt.show()

Min = int(np.argmin(Sum_1))
print(Min)
print(Time[Min])

save2folder = '/Volumes/Experiment_Data/DESY/DATA/export/' 
run = files[0]
if save2igor == 1 :
    waveName = 'unBinned_2D'
    start_at = 0
    runType = 'XPS'
    end_at = 1
    sif.save_2D_data(save2folder,Total_profile,waveName,start_at,end_at,run,runType)
    tp_Binned_profile = np.flip(Binned_profile.transpose(),axis=1)
    waveName = 'Binned_2D'
    sif.save_2D_data(save2folder,tp_Binned_profile,waveName,Delay_start,Delay_end,run,runType)
    waveName = 'unBinned_1D'
    sif.saveData(save2folder,Spectrum,E_scale,waveName,run,runType)
    waveName = 'Binned_01'
    sif.saveData(save2folder,Spectrum_1,E_scale_bin,waveName,run,runType)
    waveName = 'Binned_02'
    sif.saveData(save2folder,Spectrum_2,E_scale_bin,waveName,run,runType)
    waveName = 'Binned_03'
    sif.saveData(save2folder,Spectrum_3,E_scale_bin,waveName,run,runType)
    waveName = 'Sum_Binned_00'
    sif.saveData(save2folder,Sum_0,Time,waveName,run,runType)
    waveName = 'Sum_Binned_01'
    sif.saveData(save2folder,Sum_1,Time,waveName,run,runType)
    waveName = 'Sum_Binned_02'
    sif.saveData(save2folder,Sum_2,Time,waveName,run,runType)
    waveName = 'Sum_Binned_03'
    sif.saveData(save2folder,Sum_3,Time,waveName,run,runType)
    