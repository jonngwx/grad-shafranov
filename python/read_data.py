import numpy as np
import io
import field
import h5py

FLOAT='[0-9]*\.[0-9]*'

DATA= '([A-Za-z]): ((' + FLOAT+ ' )*)'

def read_hdf5(filename):
    f = h5py.File(filename,'r')
    x = f['R'].shape[0]
    y = f['z'].shape[0]
    fields = field.Field(x,y)
    fields.R = np.array(f['R'])
    fields.z = np.array(f['z'])
    fields.psi = np.transpose(np.array(f['psi']))
    fields.p = np.transpose(np.array(f['p']))
    f.close()
    return fields


def read_text(filename):
    """ format assumed to be R\n,z\n, array\n, array\n"""
    f = io.open(filename,'r')
    s = f.readline()
    R = np.fromstring(s, sep=' ')
    s = f.readline()
    z = np.fromstring(s, sep=' ')
    x = R.shape[0]
    y = z.shape[0]
    fields = field.Field(x,y)
    fields.R = R
    fields.z = z
 #   print x
 #   print y
    for i in xrange(x):
        s = f.readline()
        fields.psi[i,::] = np.fromstring(s, sep = ' ')
    fields.psi = np.transpose(fields.psi)
    for i in xrange(x):
        s = f.readline()
        fields.p[i,::] = np.fromstring(s, sep = ' ')
    fields.p = np.transpose(fields.p)
    f.close()
#    print R
#    print z
#    print psi
#    print p
    return fields
