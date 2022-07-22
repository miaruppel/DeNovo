import numpy
import collections
from constants import AMINO_ACID, PROTON, ION_OFFSET, FORWARD, BACKWARD
import constants


def adjust_masses(method):
    if method == "SILAC":
        offsets = {"K": 8.01419881319, "R": 10.008268599}
    else:
        raise ValueError("Don't know method: " + method)

    for aa, offset in offsets.items():
        AMINO_ACID[aa] += offset


def get_mz(sum_, ion_offset, charge):
    return (sum_ + ion_offset + charge * PROTON) / charge


def get_mzs(cumsum, ion_type, z, fragment_type):
    
    if fragment_type == 'y':
        return [get_mz(s, ION_OFFSET[ion_type], z) for s in cumsum[:-1]]
    elif fragment_type == 'b':
        if ion_type == 'y':
            return [get_mz(s, ION_OFFSET[ion_type]-18.010564, z) for s in cumsum[:-1]]
        else:
            return [get_mz(s, ION_OFFSET[ion_type], z) for s in cumsum[:-1]]


def get_annotation(forward, backward, charge, ion_types, fragment_type):
    tmp = "{}{}"
    tmp_nl = "{}{}-{}"
    all_ = {}
    for ion_type in ion_types:
        if ion_type in constants.FORWARD:
            cummass = forward
        elif ion_type in constants.BACKWARD:
            cummass = backward
        else:
            raise ValueError("unknown ion_type: {}".format(ion_type))
        masses = get_mzs(cummass, ion_type, charge, fragment_type)
        #print(masses)
        d = {tmp.format(ion_type, i + 1): m for i, m in enumerate(masses)}
        all_.update(d)
        #for nl, offset in constants.NEUTRAL_LOSS.items():
            #nl_masses = get_mzs(cummass - offset, ion_type, charge)
            #d = {tmp_nl.format(ion_type, i + 1, nl): m for i, m in enumerate(nl_masses)}
            #all_.update(d)
    return collections.OrderedDict(sorted(all_.items(), key=lambda t: t[0]))
