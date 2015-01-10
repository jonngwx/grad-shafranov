import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import sys
#import field
import read_data
import re

def plot(filename,format):
    """
    This is a function which plots things. It takes two strings, one of which is the filename, followed by the format.
    @return F a dict containing all the data in the file.
    """
    try:
        if format == "txt" or format == "tsv":
            F = read_data.read_text(filename)
        elif format == "h5" or format == "hdf" or format == "hdf5":
            F = read_data.read_hdf5(filename)
        else:
            print "Invalid format\n"
            return None
    except IOError:
        print "Invalid filename\n"
        return None
    if F is None:
        print "Error reading file " + filename
        return None
    for i in F.keys():
        if np.isnan(np.sum(F[i])):
            print "Warning: nan in data, exiting."
            return F
    R0 = F['R'][0]
    Rend = F['R'][-1]
    z0 = F['z'][0]
    zend = F['z'][-1]
    nR = F['R'].shape[0]
    nz = F['z'].shape[0]
    dr = (Rend-R0)/(nR-1)
    dz = (zend-z0)/(nz-1)
    fig = plt.figure(1,figsize=(16,9),dpi = 80)
    cfig = plt.subplot(1,3,1)
    clist = np.linspace(np.min(F['psilo']),np.max(F['psilo']),10)
    matplotlib.rcParams['contour.negative_linestyle']='solid'
    plt.contour(F['R'],F['z'],F['psi'],clist,colors='k')
    plt.contour(F['R'],F['z'],F['psi'],[F['psilo'][0],F['psilo'][0]],colors='r')
    plt.pcolormesh(F['R'],F['z'],F['psi'])
    def format_coord_psi(x,y):
        col = int((x-R0)/dr)
        row = int((y-z0)/dz)
        if (col >=0 and col < nR and row >= 0 and row < nz):
            z = F['psi'][row,col]
            return 'R = %1.4f, z = %1.4f, Psi = %1.4f'%(x,y,z)
        else:
            return 'R = %1.4f, z = %1.4f'%(x,y)

    cfig.format_coord = format_coord_psi
    cb = plt.colorbar()
    plt.title('$\Psi$')
    plt.xlabel('$R$ (m)')
    plt.ylabel('z (m)')
    #plt.xticks([])
    cfig.set_aspect('equal')
    a = np.amin(np.abs(F['psi'] - F['psilo'][1]))
    ind = np.where(abs(F['psi'] - F['psilo'][1]) - a == 0)
    if ind is not None:
        try:
            plt.plot(F['R'][ind[1]], F['z'][ind[0]],'+',markersize=10,mew=5.0)
        except:
            print "Cannot find magnetic axis"
    pfig = plt.subplot(1,3,2)
    plt.pcolormesh(F['R'],F['z'],F['p'])#,shading='gouraud')
    def format_coord_p(x,y):
        col = int((x-R0)/dr)
        row = int((y-z0)/dz)
        if (col >=0 and col < nR and row >= 0 and row < nz):
            z = F['p'][row,col]
            return 'R = %1.4f, z = %1.4f, p = %1.4f'%(x,y,z)
        else:
            return 'R = %1.4f, z = %1.4f'%(x,y)
    pfig.format_coord = format_coord_p
    plt.title('p')
    plt.xlabel('$R$ (m)')
    plt.xlim([R0,Rend])
    plt.ylim([z0,zend])
    plt.colorbar()
    pfig.set_aspect('equal')


    gfig = plt.subplot(1,3,3)
    plt.pcolormesh(F['R'],F['z'],F['g'])#,shading='gouraud')
    plt.title('g')
    def format_coord_g(x,y):
        col = int((x-R0)/dr)
        row = int((y-z0)/dz)
        if (col >=0 and col < nR and row >= 0 and row < nz):
            z = F['g'][row,col]
            return 'R = %1.4f, z = %1.4f, g = %1.4f'%(x,y,z)
        else:
            return 'R = %1.4f, z = %1.4f'%(x,y)
    gfig.format_coord = format_coord_g
    plt.xlabel('$R$ (m)')
    plt.xlim([R0,Rend])
    plt.ylim([z0,zend])
    plt.colorbar()
    gfig.set_aspect('equal')
    #positional hackery
    cpos = cfig.get_position().get_points()
    gpos = gfig.get_position().get_points()
   # print cpos
   # print gpos
    cpos[0,0] = gpos[0,0]
    cpos[1,0] = gpos[1,0]-gpos[0,0]
    cpos[1,1] = cpos[1,1]-cpos[0,1]
   # print cpos
   # cfig.set_position(np.ndarray.flatten(cpos))
   # print cfig.get_position().get_points()
   # print gfig.get_position().get_points()
    fig.delaxes(fig.axes[1])
    plt.show()
    return F


def midplane_plot(F,x):
    """
    A radial plot of a quantity along the midplane of the tokamak
    @param F container holding fields and grid
    @param x the quantity to plot
    @return R the grid in the radial direction
    @return a the field x along the midplane, or None if invalid
    """
    if x not in F.keys():
        print "The specified quantity is not available"
        return F['R'],None
    a = F[x]
    if np.array(a.shape).shape[0] != 2:
        print "The specified quantity is not a 2d array"
        return F['R'],None
    nz = a.shape[0]
    plt.close(1)
    plt.figure(1)
    plt.clf()
    plt.plot(F['R'],a[nz/2,::])
    plt.show()
    return F['R'], a[nz/2,::]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'Format is viz <filename>'
    else:
        FILE = '.*\.(hdf5|tsv)'
        prog = re.compile(FILE)
        result = prog.match(sys.argv[1])
        print "Reading file %s ..."%sys.argv[1]
        if result is None:
            print "invalid file name %s"%sys.argv[1]
            sys.exit()
        plot(sys.argv[1],result.group(1))
