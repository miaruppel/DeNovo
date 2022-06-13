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


# In[6]:


# parse function
def parseScanLine(input):
    x = input.split(" For: ")
    [scan_number, mzs] = x[1].split(", ")
    [precursor_mz, fragment_mz] = mzs.split(";")
    trimmed_fragment_mz = fragment_mz.strip() # trim fragment strings to remove \n
    return [scan_number, precursor_mz, trimmed_fragment_mz]


# In[7]:


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


# In[8]:


# load in MS2 scans
exp1 = MSExperiment()
MzMLFile().load("filtered_MS2.mzML", exp1)


# In[9]:


# read in peptide sequence from tsv
import pandas as pd
tsv = pd.read_csv('HEK293T_De_Novo_053122_Glu-C_B_correctRTSenzyme_BP_realtimesearch1.tsv', sep='\t')

# create dictionary with scan # as key and sequence/charge as values
scan2PeptideCharge = dict([(i, [x,y]) for i, x,y, in zip(tsv['Scan Number'], tsv['Peptide'], tsv['Charge State'])])

# removing all NaN sequences (not useful)
scan2PeptideCharge_modified = {k:v for k,v in scan2PeptideCharge.items() if str(v[0]) != 'nan'}


# In[10]:


def findFragments(peptide_object, charge):
    # loop through each prefix and suffix (b and y ions, respectively)
    # y and b ions
    y_ions = []
    b_ions = []
    for i in range(1, (peptide_object.size())): # start at index of 1, end at peptide length - 1
        y_ions.append(peptide_object.getSuffix(i))
        b_ions.append(peptide_object.getPrefix(i))
        
    def loopChargeStates():
        # computing fragment mzs
        # compute all y_ion mzs
        y_ion_mzs = []
        for i in y_ions[::-1]:
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
    
    # call nested function
    y_ion_mzs, b_ion_mzs = loopChargeStates()

    y_indices = []
    for i in y_ion_mzs:
        y_indices.append(s.findNearest(i, 0.4))
    b_indices = []
    for i in b_ion_mzs:
        b_indices.append(s.findNearest(i, 0.4))

    return y_indices, b_indices


# In[11]:


y_indices = []
b_indices = []
for s in exp1:
    s_number = s.getNativeID().split(' ')[-1]
    _, scan_number = s_number.split('=')
    
    if scan_number in scan2frag and int(scan_number) in scan2PeptideCharge_modified:
        # isolate peptide sequence from dict
        sequence = scan2PeptideCharge_modified[int(scan_number)][0]     
        trimmed_sequence = sequence[2:-2] # remove first two and last two characters 
        
        # isolate charge from dict
        charge = scan2PeptideCharge_modified[int(scan_number)][1]
        
        # create peptide object 
        peptide_object = AASequence.fromString(trimmed_sequence)
        
        # call findFragments function
        y, b = findFragments(peptide_object, charge)
        y_indices.append(y)
        b_indices.append(b)
    
    else: 
        pass 


# In[12]:


spectrum_list = []
scans_list = []
for s in exp1:
    s_number = s.getNativeID().split(' ')[-1]
    _, scan_number = s_number.split('=')
    
    if scan_number in scan2frag and int(scan_number) in scan2PeptideCharge_modified:
        spectrum_list.append(s)
        scans_list.append(scan_number)
        


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


# the following code is checking for how many fragments are missing from just the first peptide sequence 


# In[16]:


sequence = scan2PeptideCharge_modified[1707][0]     
trimmed_sequence = sequence[2:-2] # remove first two and last two characters       


# In[17]:


trimmed_sequence


# In[18]:


peptide_object = AASequence.fromString(trimmed_sequence)


# In[19]:


peptide_object.size()


# In[20]:


y_ions = []
b_ions = []
for i in range(1, (peptide_object.size())): # start at index of 1, end at peptide length - 1
    y_ions.append(peptide_object.getSuffix(i))
    b_ions.append(peptide_object.getPrefix(i))
        


# In[21]:


# b1
print(b_ions[0])


# In[22]:


# y16
print(y_ions[::-1][0]) # reverse the list indexing 


# In[23]:


charge = scan2PeptideCharge_modified[1707][1]
print(charge)


# In[ ]:





# In[ ]:





# In[ ]:





# In[24]:


# finding bs
def testingBs(fragment_index):
    global mz1, indb
    mz1 = []
    for charge_state in range(1, charge):
        mz1.append(b_ions[fragment_index].getMonoWeight(Residue.ResidueType.BIon, charge_state) / charge_state)
        
test1 = []
for i in range(len(b_ions)):
    testingBs(i)
    test1.append(mz1)
    
indb = []
for z in range(len(test1)):
    for i in test1[z]:
        indb.append(spectrum_list[0].findNearest(i, 0.4))

# call the function
testingBs(0)


# In[29]:


# finding ys
def testingYs(fragment_index):
    global mz2, indy
    mz2 = []
    for charge_state in range(1, charge):
        mz2.append(y_ions[fragment_index].getMonoWeight(Residue.ResidueType.YIon, charge_state) / charge_state)
        
test2 = []
for i in range(len(y_ions)):
    testingYs(i)
    test2.append(mz2)

indy = []
for z in range(len(test2)):
    for i in test2[z]:
        indy.append(spectrum_list[0].findNearest(i, 0.4))

# call the function
testingYs(0)


# In[39]:


# skim list for each corresponding y and b fragments based on charge states 
def skimList(fragment_index, charge):
    checkb = indb[fragment_index:fragment_index+charge-1]
    checky = indy[fragment_index:fragment_index+charge-1]
    return checkb, checky

# check after every two mz (two possible charge states), so index as such: 0,2,4,6,etc
checkb, checky = skimList(0, 3)

if all(item == checkb[0] for item in checkb and checky):
    print('All elements have a value of -1...we cannot uncover this fragment')
else: 
    print('At least one ion (b or y) has at least one charge state that was found in the MS2')


# In[38]:


# for every fragment 
for x in range(0, len(indy), 2):
    checkb, checky = skimList(x, charge)
    if all(item == checkb[0] for item in checkb and checky):
        print('All elements have a value of -1...we cannot uncover this fragment')
    else: 
        print('At least one ion (b or y) has at least one charge state that was found in the MS2')


# In[ ]:





# In[ ]:





# In[ ]:




