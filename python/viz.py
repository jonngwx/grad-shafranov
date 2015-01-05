import numpy as np
import matplotlib.pyplot as plt
import sys
#import field
import read_data

def plot(filename,format):
    """Please describe this function

    in some detail.
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
    fig = plt.figure(1,figsize=(4,9),dpi = 80)
    cfig = plt.subplot(3,1,1)
    plt.contour(F['R'],F['z'],F['psi'],colors='k')
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
    plt.xticks([])
    plt.ylabel('z')
    cfig.set_aspect('equal')

    pfig = plt.subplot(3,1,2)
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
    plt.ylabel('z')
    plt.xlim([R0,Rend])
    plt.colorbar()
    pfig.set_aspect('equal')


    gfig = plt.subplot(3,1,3)
    plt.pcolormesh(F['R'],F['z'],F['g'])#,shading='gouraud')
    plt.title('g')
    def format_coord_g(x,y):
        col = int((x-R0)/dr)
        row = int((y-z0)/dz)
        if (col >=0 and col < nR and row >= 0 and row < nz):
            z = F['psi'][row,col]
            return 'R = %1.4f, z = %1.4f, g = %1.4f'%(x,y,z)
        else:
            return 'R = %1.4f, z = %1.4f'%(x,y)
    gfig.format_coord = format_coord_g
    plt.xlabel('$R$')
    plt.ylabel('z')
    plt.xlim([R0,Rend])
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
    cfig.set_position(np.ndarray.flatten(cpos))
   # print cfig.get_position().get_points()
   # print gfig.get_position().get_points()
    fig.delaxes(fig.axes[1])
    plt.show()
    return F


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print 'wrong number of arguments'
    else:
        plot(sys.argv[1],sys.argv[2])
