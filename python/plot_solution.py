import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')


path='/home/dmunoz/Documents/REPOSITORIES/diskonedee/diskonedee/'
filename='output.txt'


if __name__ == '__main__':

    data=np.loadtxt(path+filename)
    x = data[0,:]

    fig = plt.figure(figsize=(8,5))
    fig.subplots_adjust(bottom=0.12,right=0.97,top=0.97,left=0.1)
    ax = fig.add_subplot(111)
    
    ax.plot(x,data[1,:],color='r',lw=1.8)
    for sol in data[2:,:]:
        ax.plot(x,sol,color='k',lw=0.3)

    ax.set_xlabel(r'$R$',size=24,labelpad=0)
    ax.set_ylabel(r'$\Sigma$',size=24,labelpad=0)
    ax.set_xlim(1,100)
    ax.set_ylim(0,0.04)
    ax.tick_params(axis='both',labelsize=16)
    
    fig.savefig('diffusion.png')
    plt.show()
        
