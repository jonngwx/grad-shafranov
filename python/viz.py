import numpy as np
import matplotlib.pyplot as plt
import sys
#import field
import read_data

def plot(filename,format):
    try:
        if format == "txt":
            F = read_data.read_text(filename)
        elif format == "h5":
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
    plt.figure(1)
    cfig = plt.subplot(3,1,1)
    plt.contour(F['R'],F['z'],F['psi'],colors='k')
    plt.title('$\Psi$')
    plt.xticks([])
    plt.ylabel('z')
#    plt.colorbar()
    plt.subplot(3,1,2)
    plt.pcolormesh(F['R'],F['z'],F['p'])
    plt.title('p')
#    plt.xlabel('$R$')
    plt.ylabel('z')
    plt.colorbar()
    gfig = plt.subplot(3,1,3)
    plt.pcolormesh(F['R'],F['z'],F['g'])
    plt.title('g')
    plt.xlabel('$R$')
    plt.ylabel('z')
    plt.colorbar()
    #positional hackery
    cpos = cfig.get_position().get_points()
    gpos = gfig.get_position().get_points()
    cpos[1,0] = gpos[1,0]-gpos[0,0]
    cpos[1,1] = cpos[1,1]-cpos[0,1]
    cfig.set_position(np.ndarray.flatten(cpos))
    print cfig.get_position()
    plt.show()
    return F

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print 'wrong number of arguments'
    else:
        plot(sys.argv[1],sys.argv[2])
