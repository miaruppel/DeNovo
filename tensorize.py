import collections
import numpy as np

import constants
import utils
import match
import annotate
import sanitize
from constants import (
    CHARGES,
    MAX_SEQUENCE,
    ALPHABET,
    MAX_ION,
    NLOSSES,
    CHARGES,
    ION_TYPES,
    ION_OFFSET,
)

def stack(queue):
    listed = collections.defaultdict(list)
    for t in queue.values():
        if t is not None:
            for k, d in t.items():
                listed[k].append(d)
    stacked = {}
    for k, d in listed.items():
        if isinstance(d[0], list):
            stacked[k] = [item for sublist in d for item in sublist]
        else:
            stacked[k] = np.vstack(d)
    return stacked


def get_numbers(vals, dtype=float):
    a = np.array(vals).astype(dtype)
    return a.reshape([len(vals), 1])


def get_precursor_charge_onehot(charges):
    array = np.zeros([len(charges), max(CHARGES)], dtype=int)
    for i, precursor_charge in enumerate(charges):
        array[i, precursor_charge - 1] = 1
    return array


def get_sequence_integer(sequences):
    array = np.zeros([len(sequences), MAX_SEQUENCE], dtype=int)
    for i, sequence in enumerate(sequences):
        for j, s in enumerate(utils.peptide_parser(sequence)):
            array[i, j] = ALPHABET[s]
    return array


def parse_ion(string):
    ion_type = ION_TYPES.index(string[0])
    if ("-") in string:
        ion_n, suffix = string[1:].split("-")
    else:
        ion_n = string[1:]
        suffix = ""
    return ion_type, int(ion_n) - 1, NLOSSES.index(suffix)


def get_mz_applied(df, ion_types="yb"):
    ito = {it: ION_OFFSET[it] for it in ion_types}

    def calc_row(row):
        array = np.zeros([MAX_ION, len(ION_TYPES), len(NLOSSES), len(CHARGES)])
        fw, bw = match.get_forward_backward(row.modified_sequence)
        for z in range(row.precursor_charge):
            zpp = z + 1
            annotation = annotate.get_annotation(fw, bw, zpp, ito)
            for ion, mz in annotation.items():
                it, _in, nloss = parse_ion(ion)
                array[_in, it, nloss, z] = mz
                #print(mz)
        return [array]

    mzs_series = df.apply(calc_row, 1)
    #print(mzs_series[0])
    out = np.squeeze(np.stack(mzs_series))
    if len(out.shape) == 4:
        out = out.reshape([1] + list(out.shape))
    return out


def csv(df):
    df.reset_index(drop=True, inplace=True)
    assert "modified_sequence" in df.columns
    assert "collision_energy" in df.columns
    assert "precursor_charge" in df.columns
    data = {
        "collision_energy_aligned_normed": get_numbers(df.collision_energy) / 100.0,
        "sequence_integer": get_sequence_integer(df.modified_sequence),
        "precursor_charge_onehot": get_precursor_charge_onehot(df.precursor_charge),
        "masses_pred": get_mz_applied(df),
    }
    nlosses = 1
    z = 3
    lengths = (data["sequence_integer"] > 0).sum(1)

    masses_pred = get_mz_applied(df)
    masses_pred = sanitize.cap(masses_pred, nlosses, z)
    masses_pred = sanitize.mask_outofrange(masses_pred, lengths)
    masses_pred = sanitize.mask_outofcharge(masses_pred, df.precursor_charge)
    masses_pred = sanitize.reshape_flat(masses_pred)
    data["masses_pred"] = masses_pred

    return data

def csv_training(df, series):
    df.reset_index(drop=True, inplace=True)
    assert "modified_sequence" in df.columns
    assert "collision_energy" in df.columns
    assert "precursor_charge" in df.columns
    data = {
        "collision_energy_aligned_normed": get_numbers(df.collision_energy) / 100.0,
        "sequence_integer": get_sequence_integer(df.modified_sequence),
        "precursor_charge_onehot": get_precursor_charge_onehot(df.precursor_charge),
        "masses_pred": get_mz_applied(df),
    }
    nlosses = 1
    z = 3
    lengths = (data["sequence_integer"] > 0).sum(1)

    masses_pred = get_mz_applied(df)
    masses_pred = sanitize.cap(masses_pred, nlosses, z)
    masses_pred = sanitize.mask_outofrange(masses_pred, lengths)
    masses_pred = sanitize.mask_outofcharge(masses_pred, df.precursor_charge)
    masses_pred = sanitize.reshape_flat(masses_pred)
    data["masses_pred"] = masses_pred
    print(masses_pred)

    # find matching mzs and replace with intensity values
    masses_pred_copy = masses_pred # to overwrite with intensities later
    # yeah python is weird about this actually so fix this!!!!
    
    
    
    df_masses = list(df['masses_raw']) # list of raw masses from dataframe
    df_intensities = list(df['intensities_raw']) # list of raw intensities from dataframe
    
    s_list = series.tolist() # list of MS3 spectra
    s_count = -1 # for indexing spectra
    #print(masses_pred_copy[0:3])
    
    # for finding closer mass to match spectra
    def closest(lst, K):
        return lst[min(range(len(lst)), key = lambda i: abs(lst[i]-K))]
    
    for mass_list in masses_pred_copy: # mass list for each spectrum
        s_count = s_count + 1
        mz_nums = df_masses[s_count].split(' ')
        df_mass_list = [float(num) for num in mz_nums]
        
        int_nums = df_intensities[s_count].split(' ')
        df_intens_list = [float(num) for num in int_nums]
        #print(df_intens_list)
        
        m_count = -1 # for indexing masses
        for mass in mass_list: # looking at individual masses
            m_count = m_count + 1
            #print(df_intens_list[s_count])
            if mass != -1: 
                test = s_list[s_count].findNearest(mass, 0.02) 
                if test != -1: # masses found 
                    
                    # bounds for matching with raw masses 
                    upper = round(mass + 0.02, 2)
                    lower = round(mass - 0.02, 2)
                    #print(mass, upper, lower)
                    index = [i for i, value in enumerate(df_mass_list) if (round(value, 2) <= upper and round(value, 2) >= lower)]
                    # if more than one index, we want the one that is closest to to the target
                    if len(index) > 1:
                        closer = closest(df_mass_list, mass)
                        for i in index:
                            if closer == df_mass_list[i]:
                                index = [i]
                    # use index to match intensity 
                    intensity = df_intens_list[int(str(index)[1:-1])] # removing brackets
                    # overwrite mass values with intensities 
                    masses_pred_copy[s_count][m_count] = intensity
                else:
                    masses_pred_copy[s_count][m_count] = -1
    print(masses_pred)
                 
    intensities_raw = masses_pred_copy
    #print(intensities_raw)
    #intensities_raw[intensities_raw < 0] = 0
    intensities_raw = sanitize.normalize_base_peak(intensities_raw)
    #intensities_raw = reshape_dims(intensities_raw)
    #intensities_raw = sanitize.mask_outofrange(intensities_raw, lengths)
    #intensities_raw = sanitize.mask_outofcharge(intensities_raw, df.precursor_charge)
    #intensities_raw = sanitize.reshape_flat(intensities_raw)
    data["intensities_raw"] = intensities_raw

    #for mass in masses_pred:
     #   print(x)
     #   for y in x:
     #       print(y)
        
    #intensities[intensities < 0] = 0
    #intensities = normalize_base_peak(intensities)
    #intensities = reshape_dims(intensities)
    #intensities = mask_outofrange(intensities, sequence_lengths)
    #intensities = mask_outofcharge(intensities, charges)
    #intensities = reshape_flat(intensities)
    #data["intensities_pred"] = intensities

    return data
