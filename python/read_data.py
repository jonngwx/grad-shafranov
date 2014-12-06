import numpy as np
import io
import field
#import h5py

def read_hdf5(filename):
    return 0


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
