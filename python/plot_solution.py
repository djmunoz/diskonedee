import numpy as np
import matplotlib.pyplot as plt



path='/home/dmunoz/Documents/REPOSITORIES/diskonedee/diskonedee/'
filename='output.txt'


if __name__ == '__main__':

    data=np.loadtxt(path+filename)
    print(data.shape)
    x = data[0,:]
    plt.plot(x,data[1,:],color='r')
    for sol in data[2:,:]:
        plt.plot(x,sol,color='k',lw=0.1)

    #plt.ylim(0,70.0)
    plt.show()
        
