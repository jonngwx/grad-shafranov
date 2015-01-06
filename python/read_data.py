import numpy as np
import io
import re

FLOAT='-?[0-9]*\.[0-9]*'
NAN = '-?nan'
DATA= '([A-Za-z]+): *((' + FLOAT+ ' *|'+ NAN +' *)*)'

def read_hdf5(filename):
    """Reads in an hdf5 file.
    """
    import h5py
    f = h5py.File(filename,'r')
    fields = {}
    for data in f.keys():
        if data == 'R' or data == 'z':
            fields[data]=np.array(f[data])
        else:
            fields[data]=np.transpose(np.array(f[data]))
            if data == 'psi':
                p0 = f[data].attrs['psi_0']
                pl = f[data].attrs['psi_l']
                fields['psilo'] = np.array([pl, p0])
    f.close()
    if is_invalid(fields):
        return None
    return fields


def read_text(filename):
    """ format assumed to be R\n,z\n, array\n, array\n"""
    f = io.open(filename,'r')
    prog = re.compile(DATA)
    fields = {}
    for line in f:
        result = prog.match(line)
        if result is None:
            continue
        fields[result.group(1)] = np.fromstring(result.group(2), sep=' ')
    f.close()
    if is_invalid(fields):
        return None
    nx = fields['R'].shape[0]
    ny = fields['z'].shape[0]
    for data in fields.keys():
        if data == 'R' or data == 'z' or data == 'psilo':
            continue
        else:
            fields[data] = np.transpose(np.reshape(fields[data],[nx,ny]))
#    print R
#    print z
#    print psi
#    print p
    return fields

def is_invalid(fields):
    k = fields.keys()
    invalid = False
    data = ['R','z','psi','p','g']
    for i in data:
        invalid = (invalid or check_data(i,k))
    return invalid

def check_data(name,keys):
    if name not in keys:
        print name + " not in fields"
        return True
    else:
        return False
