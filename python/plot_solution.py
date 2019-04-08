import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')


path='../diskonedee/'


ZAMFB = True

if (ZAMFB):
    filename='output_zamfb.txt'
else:
    filename='output_zerotorque.txt'

    
if __name__ == '__main__':

    data=np.loadtxt(path+filename)
    x = data[0,:]

    fig = plt.figure(figsize=(8,5))
    fig.subplots_adjust(bottom=0.12,right=0.97,top=0.97,left=0.13)
    ax = fig.add_subplot(111)

    ax.plot(x,data[1,:],color='r',lw=1.8)
    for sol in data[2:,:]:
        dsigmadr=np.gradient(sol)/np.gradient(x)
        #ax.plot(x,sol/dsigmadr,color='k',lw=0.3)
        ax.plot(x,sol,color='k',lw=0.3)

    ax.set_xlabel(r'$R$',size=24,labelpad=0)
    ax.set_ylabel(r'$\Sigma$',size=24,labelpad=0)
    ax.set_xlim(1,100)
    ax.set_ylim(0.001,0.05)
    ax.tick_params(axis='both',labelsize=16)
    ax.set_xscale('log')
    ax.set_yscale('log')
    if (ZAMFB):
        ax.plot(x,x**(-0.5)*0.03,color='royalblue',zorder=0,lw=4.0,alpha=0.6)
        ax.text(40,40**(-0.5)*0.03,r'$\propto R^{-1/2}$',
                transform=ax.transData,size=22,color='royalblue')
        
    if (ZAMFB): 
        fig.savefig('diffusion_zamfb.png')
    else:
        fig.savefig('diffusion_zerotorque.png')
    plt.show()
        
