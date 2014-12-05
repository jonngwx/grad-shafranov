import numpy as np
import matplotlib.pyplot as plt
import sys
#import field
import read_data

def plot(filename,format):
    if format == "txt":
        F = read_data.read_text(filename)
    else:
        1+1
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.pcolormesh(F.R,F.z,np.transpose(F.psi))
    plt.title('$\Psi$')
    plt.xticks([])
    plt.ylabel('z')
    plt.colorbar()
    plt.subplot(2,1,2)
    plt.pcolormesh(F.R,F.z,np.transpose(F.p))
    plt.title('p')
    plt.xlabel('$R$')
    plt.ylabel('z')
    plt.colorbar()
    plt.show()
    return F

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print 'wrong number of arguments'
    else:
        plot(sys.argv[1],sys.argv[2])
