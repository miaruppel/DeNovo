#!/usr/bin/env python
# coding: utf-8

# In[1]:


# import libraries
from pyopenms import *
import os


# In[2]:


# change directory to find file of interest
os.chdir(r'C:\Users\miar\Desktop\data')


# In[3]:


#------------------------------------------------------------------------------------------------------------------------------

# ONLY NEED TO RUN THE NEXT TWO BLOCKS ONCE

#-------------------------------------------------------------------------------------------------------------------------------


# In[4]:


# load the content of the mzML file into the exp variable of type MSExperiment
#exp = MSExperiment()
#MzMLFile().load("HEK293T_De_Novo_053122_Glu-C_B_correctRTSenzyme_BP.mzML", exp)

# access the raw data and spectra
#spectrum_data = exp.getSpectrum(0).get_peaks()
#spectrum_data


# In[5]:


# loop through the spectra to gather MS2 scans
#specM2 = []
#for s in exp.getSpectra():
#    if s.getMSLevel() == 2:
#        specM2.append(s)
        
#print("Number of MS2 scans: " + str(len(specM2)))

#exp.setSpectra(specM2) # keep only MS2

# store the modified data structure on disk
#MzMLFile().store("filtered_MS2.mzML", exp)


# In[3]:


# parse function
def parseScanLine(input):
    x = input.split(" For: ")
    [scan_number, mzs] = x[1].split(", ")
    [precursor_mz, fragment_mz] = mzs.split(";")
    trimmed_fragment_mz = fragment_mz.strip() # trim fragment strings to remove \n
    return [scan_number, precursor_mz, trimmed_fragment_mz]


# In[6]:


# checking lines of log file and creating dictionary of scan numbers and fragment mzs
try:
  
    # words to search for
    search = ' Submitted Custom Scan For:'
  
    # reading file content line by line
    search = ' Submitted Custom Scan For:'   # words to search for
    
    # dict for scan numbers and corresponding fragments 
    scan2frag = dict()
    with open('App-2022-05-31_20-49-35.log') as f:
        for line in f:
            if search in line:
                scan_number, precursor_mz, trimmed_fragment_mz = parseScanLine(line)
                scan2frag[scan_number] = trimmed_fragment_mz
            
    # if the input string doesn't exist in the text file
    if len(scan2frag)==0:
        print("\n\"" +search+ "\" is not found in \"" +'App-2022-05-31_20-49-35.log'+ "\"!")
    else:
        pass

except FileNotFoundError:
    print("The file does not exist!")


# In[9]:


# load in MS2 scans
exp1 = MSExperiment()
MzMLFile().load("filtered_MS2.mzML", exp1)


# In[10]:


# read in peptide sequence from tsv
import pandas as pd
tsv = pd.read_csv('HEK293T_De_Novo_053122_Glu-C_B_correctRTSenzyme_BP_realtimesearch1.tsv', sep='\t')

# create dictionary with scan # as key and sequence/charge as values
scan2PeptideCharge = dict([(i, [x,y]) for i, x,y, in zip(tsv['Scan Number'], tsv['Peptide'], tsv['Charge State'])])

# removing all NaN sequences (not useful)
scan2PeptideCharge_modified = {k:v for k,v in scan2PeptideCharge.items() if str(v[0]) != 'nan'}


# In[18]:


def findFragments(peptide_object, charge):
    # loop through each prefix and suffix (b and y ions, respectively)
    try:
        # y and b ions
        y_ions = []
        b_ions = []
        for i in range(1, (len(trimmed_sequence) - 1)): # start at index of 1, end at peptide length - 1
            y_ions.append(peptide_object.getSuffix(i))
            b_ions.append(peptide_object.getPrefix(i))

    except RuntimeError: # range above may be too large for indexing when considering modifications (ex. [15.9949])
        print('Modifications resulted in abnormal indexing for sequence: ' + str(peptide_object))

        # will have to remove "fragment" that is now the entire peptide length instead of peptide length - 1
        y_ions.pop()
        b_ions.pop()
        
    def loopChargeStates():
    
        # computing fragment mzs
        # compute all y_ion mzs
        y_ion_mzs = []
        for i in y_ions:
            for x in range(1, charge):
                mz_y = i.getMonoWeight(Residue.ResidueType.YIon, x) / x
                y_ion_mzs.append(mz_y)   
        # compute all b_ion mzs
        b_ion_mzs = []
        for i in b_ions:
            for x in range(1, charge):
                mz_b = i.getMonoWeight(Residue.ResidueType.BIon, x) / x
                b_ion_mzs.append(mz_b)
        return y_ion_mzs, b_ion_mzs
    
    y_ion_mzs, b_ion_mzs = loopChargeStates()


    y_indices = []
    for i in y_ion_mzs:
        y_indices.append(s.findNearest(i, 0.4))
    b_indices = []
    for i in b_ion_mzs:
        b_indices.append(s.findNearest(i, 0.4))

    return y_indices, b_indices


# In[24]:


for s in exp1:
    get_scan = s.getNativeID()
    c_type, c_number, s_number = get_scan.split(' ')
    words, scan_number = s_number.split('=')
    
    if scan_number in scan2frag and int(scan_number) in scan2PeptideCharge_modified:
        # isolate peptide sequence from dict
        sequence = scan2PeptideCharge_modified[int(scan_number)][0]     
        trimmed_sequence = sequence[2:-2] # remove first two and last two characters 
        
        # isolate charge from dict
        charge = scan2PeptideCharge_modified[int(scan_number)][1]
        
        # create peptide object 
        peptide_object = AASequence.fromString(trimmed_sequence)
        
        # call findFragments function
        y_indices, b_indices = findFragments(peptide_object, charge)
    
    else: 
        pass 
        


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[14]:


# read in peptide sequence from tsv
import pandas as pd
tsv = pd.read_csv('HEK293T_De_Novo_053122_Glu-C_B_correctRTSenzyme_BP_realtimesearch1.tsv', sep='\t')

# create dictionary with scan # as key and sequence/charge as values
dict2 = dict([(i, [x,y]) for i, x,y, in zip(tsv['Scan Number'], tsv['Peptide'], tsv['Charge State'])])

# removing all NaN sequences (not useful)
dict2_modified = {k:v for k,v in dict2.items() if str(v[0]) != 'nan'}


# In[15]:


# isolate sequences/scan numbers from dict2 for scan numbers found in both dict2 and filtered spectra
seqs = []
scan_nums = []
for i in spec_scans:
    if i in dict2_modified:
        seqs.append(dict2_modified[i][0])
        scan_nums.append(i)
        
# isolate portion of sequence between the periods       
trimmed_seqs = [] 
for i in seqs: 
    trimmed_seqs.append(i[2:-2]) # remove first two and last two characters in str 


# In[12]:


# create a peptide object for ID'd scans 
pept_objects = []
for i in trimmed_seqs:
    pept_objects.append(AASequence.fromString(i))


# In[16]:


def findFragments(input):
    # loop through each prefix and suffix (b and y ions, respectively)
    try:
        # y and b ions
        y_ions = []
        b_ions = []
        for i in range(1, (len(trimmed_seqs[input]) - 1)): # start at index of 1, end at peptide length - 1
            y_ions.append(pept_objects[input].getSuffix(i))
            b_ions.append(pept_objects[input].getPrefix(i))

    except RuntimeError: # range above may be too large for indexing when considering modifications (ex. [15.9949])
        print('Modifications resulted in abnormal indexing for sequence: ' + str(pept_objects[input])) 

        # will have to remove "fragment" that is now the entire peptide length instead of peptide length - 1
        del y_ions[-1]
        del b_ions[-1]
        
    def loopChargeStates():
        # isolate charges from dict2
        charges = []
        for i in scan_nums:
            charges.append(dict2_modified[i][1])
    
        # computing fragment mzs
        # compute all y_ion mzs
        y_ion_mzs = []
        for i in y_ions:
            for x in range(1, charges[input]):
                mz_y = i.getMonoWeight(Residue.ResidueType.YIon, x) / x
                y_ion_mzs.append(mz_y)   
        # compute all b_ion mzs
        b_ion_mzs = []
        for i in b_ions:
            for x in range(1, charges[input]):
                mz_b = i.getMonoWeight(Residue.ResidueType.BIon, x) / x
                b_ion_mzs.append(mz_b)
        return y_ion_mzs, b_ion_mzs
    
    y_ion_mzs, b_ion_mzs = loopChargeStates()
    
    # find fragments in associated spectrum
    # dict for ALL spec scan numbers and associated spectrum
    dict3 = dict(zip(spec_scans, all_spectra))
    
    # scans of interest are all in scan_nums
    dict3_modified = {k:v for k, v in dict3.items() if k in scan_nums}
    spectrum_list = list(dict3_modified.values()) # the spectra we care about

    y_indices = []
    for i in y_ion_mzs:
        y_indices.append(spectrum_list[input].findNearest(i, 0.4))
    b_indices = []
    for i in b_ion_mzs:
        b_indices.append(spectrum_list[input].findNearest(i, 0.4))

    return y_indices, b_indices


# In[14]:


# 1st sequence
y_indices, b_indices = findFragments(0)


# In[15]:


# 2nd sequence
y_indices, b_indices = findFragments(1)


# In[ ]:


# etc etc


# In[ ]:


# next....
# how many times we didn't see any peaks that correspond to the fragmentation at a specific peptide bond 


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




