#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import libraries
from pyopenms import * # main package used for handling MS data
import os # changing directories
import pandas as pd # creating and manipulating dataframe


# In[2]:


# change directory to find data files of interest
os.chdir(r'C:\Users\miar\Desktop\data')


# In[3]:


# experiment files 
mzML = 'HEK293T_De_Novo_061122_Glu-C_B_BP_anyLength_HCD10.mzML'
log = 'App-2022-06-12_14-16-26.log'
realtime = 'HEK293T_De_Novo_061122_Glu-C_B_BP_anyLength_HCD10_realtimesearch.tsv'


# In[4]:


# load the content of the mzML file into the exp variable of type MSExperiment
exp = MSExperiment()
MzMLFile().load(mzML, exp)


# In[5]:


# loop through the spectra to filter MS3 scans
specM3 = [] # list of MS3 spectra 
row_data = []
for s in exp.getSpectra():
    
    if s.getMSLevel() == 3:
        specM3.append(s)
        
        # get scan number
        s_number = s.getNativeID().split(' ')[-1]
        _, scan_number = s_number.split('=')
        
        # obtain mz and intensity values 
        mz, intensity = s.get_peaks()
        
        mz_mod = " ".join(str(m) for m in mz)
        intensity_mod = " ".join(str(i) for i in intensity)
        
        # create dict (rows of dataframe)
        data = {#'MS3_Scan':scan_number,
       'masses_raw':mz_mod,
       'intensities_raw':intensity_mod}
        
        row_data.append(data)


# In[6]:


# create series of all MS3 spectra
s_series = pd.Series(specM3)


# In[8]:


# set MS3 scans as index
df = pd.DataFrame(row_data)
#df.set_index('MS3_Scan', inplace=True)


# In[10]:


# parse functions
def parseScanLine(input):
    x = input.split(" For: ")
    [scan_number, mzs] = x[1].split(", ")
    [precursor_mz, fragment_mz] = mzs.split(";")
    trimmed_fragment_mz = fragment_mz.strip() # trim fragment strings to remove \n
    return [scan_number, precursor_mz, trimmed_fragment_mz]

def parseTargetIons(input):
    i = input.split('Target Fragment: ')
    ion = i[1].split(',')[0]
    return ion 


# In[11]:


# checking lines of log file and creating dictionary of scan numbers and fragment mzs
try:
  
    # words to search for
    search = ' Submitted Custom Scan For:'
    
    # dict for scan numbers and corresponding fragments 
    scan2frag = dict()
    with open(log) as f:
        for line in f:
            if search in line:
                scan_number, precursor_mz, trimmed_fragment_mz = parseScanLine(line)
                scan2frag[scan_number] = [float(precursor_mz), float(trimmed_fragment_mz)]
            
    # if the input string doesn't exist in the text file
    if len(scan2frag)==0:
        print("\n\"" + search + "\" is not found in \"" + log + "\"!")

except FileNotFoundError:
    print("The file does not exist!")


# In[12]:


# obtain MS3 scan numbers
# obtain precursor and fragment mzs directly from the MS3 spectrum

ms3scan2MZs = dict()
for s in specM3:
    s_number = s.getNativeID().split(' ')[-1]
    _, scan_number = s_number.split('=')
   
    fragment, precursor = s.getPrecursors()
    precursor_mz = precursor.getMZ()
    fragment_mz = fragment.getMZ()
    
    ms3scan2MZs[int(scan_number)] = [round(float(precursor_mz), 4), round(float(fragment_mz), 4)] # 4 decimal places, similar to log


# In[13]:


def matchingMS3s(ms2_mzs, ms3_mzs): # either fragment or precursor
    
    # making sure they are within 100 scans of each other
    too_far = []
    for ms2scan, ms3scan in zip(list(scan2frag), list(ms3scan2MZs)):
        scan_diff = int(ms3scan) - int(ms2scan)
        if scan_diff > 100:
            too_far.append('Scans are not within 100 scans of each other...' + 'MS2 = ' + str(ms2scan) + ' MS3 = ' + str(ms3scan))
    
    # do they not match off the bat?
    if ms2_mzs != ms3_mzs:
        # taking into consideration rounding discrepencies between the log and the spectrum
        mismatch = []
        for i in range(0, len(list(ms3scan2MZs))):
            
            precursor_diff = float(list(ms3scan2MZs.values())[i][0]) - float(list(scan2frag.values())[i][0]) 
            if precursor_diff < 0.000101 or (precursor_diff < 0 and precursor_diff > -0.000101): # because sometimes max number will be 0.0001000002 for example
                pass
            else:
                mismatch.append(i)
            
            fragment_diff = float(list(ms3scan2MZs.values())[i][1]) - float(list(scan2frag.values())[i][1]) 
            if fragment_diff < 0.000101 or (fragment_diff < 0 and fragment_diff > -0.000101): # because sometimes max number will be 0.0001000002 for example
                pass
            else:
                mismatch.append(i)
                
        # no mismatch after rounding and within 100 scans
        if len(mismatch) == len(too_far) == 0:
            print('Scans match up after taking rounding discrepencies into consideration')
            ms2_scans = list(scan2frag)
            return ms2_scans
        
        elif len(mismatch) != 0:
            print('There is mismatch at the following indicies:') # if this is the case, need to do more work...
            for i in mismatch:
                print(i) 
                
        elif len(too_far) != 0:
            print(too_far)
        
    # they match perfectly
    elif ms2_mzs == ms3_mzs:
        # within 100 scans
        if len(too_far) == 0:
            print('Scans match up perfectly!')
            ms2_scans = list(scan2frag)
            return ms2_scans
        else:
            print(too_far)


# In[14]:


# make sure that MS3 scans are in the same order as MS2 scans
ms2_scans = matchingMS3s(list(ms3scan2MZs.values()), list(scan2frag.values()))


# In[15]:


# use realtime file to obtain peptide sequence and charge
# read in peptide sequence from tsv
tsv = pd.read_csv(realtime, sep='\t')

# create dictionary with scan # as key and sequence/charge as values
scan2PeptideCharge = dict([(i, [x,y]) for i, x,y, in zip(tsv['Scan Number'], tsv['Peptide'], tsv['Charge State'])])

# removing all NaN sequences (not useful)
scan2PeptideCharge_modified = {k:v for k,v in scan2PeptideCharge.items() if str(v[0]) != 'nan'}


# In[16]:


# collect data for dataframe
seqs = []
charges = []
analyzer = []
collision = []

energy = int(realtime.split('_')[-2][-2:])

for scan in ms2_scans:
    if int(scan) in list(scan2PeptideCharge_modified):
        charge = scan2PeptideCharge_modified[int(scan)][1]
        charges.append(charge)
        
        sequence = scan2PeptideCharge_modified[int(scan)][0]     
        trimmed_sequence = sequence[2:-2] # remove first two and last two characters 
        seqs.append(trimmed_sequence)
        
        # all ms3 scans are orbit trap 
        # to be added to dataframe with other MS3 info
        analyzer.append('FTMS')
        
        # all scans have same collision energy
        collision.append(energy)


# In[17]:


# add all data columns to dataframe
df = df.assign(charge=charges, modified_sequence=seqs, mass_analyzer=analyzer, collision_energy=collision)


# In[20]:


# remove all modified sequences
for i in df.index:
    if ('[' or ']') in df['modified_sequence'][i]:
        df.drop(i, axis=0, inplace=True)
        s_series.drop(i, inplace=True)


# In[23]:


# remove all sequences more than 30 in length
for i in df.index:
    if len(df['modified_sequence'][i]) >= 30:
        df.drop(i, axis=0, inplace=True)
        s_series.drop(i, inplace=True)


# In[26]:


# importing variables and functions from other scripts
from constants import ION_TYPES, DEFAULT_MAX_CHARGE
from match import augment


# In[27]:


# running augment function
df_augmented = augment(df, ION_TYPES, DEFAULT_MAX_CHARGE)


# In[29]:


# remove all rows that are completely empty
for i in df_augmented.index:
    if df_augmented['matches_charge1'][i] == df_augmented['matches_charge2'][i] == df_augmented['matches_charge3'][i] == df_augmented['matches_charge4'][i] == df_augmented['matches_charge5'][i] == df_augmented['matches_charge6'][i]:
        df_augmented.drop(i, axis=0, inplace=True)
        s_series.drop(i, inplace=True)


# In[33]:


# to run csv function 'charge' needs to be replaced with 'precursor_charge' ?
df_augmented.rename(columns = {'charge':'precursor_charge'}, inplace = True)


# In[34]:


# import another important function from tensorize script
import tensorize as tens


# In[74]:


# for reloading with troubleshooting
import importlib
importlib.reload(tens)


# In[35]:


# create dictionary of data
data = tens.csv_training(df_augmented)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# the following code does not work properly 
# GPU configuration issues
# moved to computer cluster 


# In[26]:


# import another important function from io_local script
#from io_local import to_hdf5, from_hdf5

# convert dictionary data to hdf5 
#to_hdf5(data, 'hdf5_data2.hdf5')


# In[ ]:


# load in ML packages and created data matrix
#import tensorflow
#from tensorflow import keras
#from keras.utils import HDF5Matrix

# load in hdf5 snd create matrix
#tensor = from_hdf5('hdf5_data.hdf5', n_samples=None)


# In[ ]:


# import model 
#import training
#import model as model_lib
# load in model
#model, model_config = model_lib.load(r"C:\Users\miar\Downloads\model_fragmentation_prediction\prosit1", trained=True)


# In[ ]:


#import importlib
#importlib.reload(model_lib)


# In[ ]:


#checking for GPU


# In[ ]:


#from tensorflow.python.client import device_lib
#print(device_lib.list_local_devices())


# In[ ]:


#import tensorflow as tf 
#tf.test.gpu_device_name()


# In[ ]:


#from tensorflow.python.client import device_lib
#def get_available_devices():
#    local_device_protos = device_lib.list_local_devices()
#    return [x.name for x in local_device_protos]
#print(get_available_devices()) 

