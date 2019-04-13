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
    print(data.shape)
    x = data[0,1:]

    fig = plt.figure(figsize=(18,4))
    fig.subplots_adjust(bottom=0.15,right=0.99,top=0.97,left=0.05,wspace=.25)

    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    omega = data[1,1:]
    coeff = -2 * np.pi * x**3 * data[4,1:] * data[2,1:]
    fnu = coeff * data[6,1:]
    mdot = np.gradient(fnu,x) / data[3,1:]
    flux = mdot * omega * x**2 - fnu
    ind = mdot != 1e30

    
    ax1.plot(x,data[6,1:],color='r',lw=1.8)
    ax2.plot(x[ind],mdot[ind],color='r',lw=1.8)
    ax3.plot(x[ind],flux[ind],color='r',lw=1.8)
    
    for sol in data[7:,1:]:
        ax1.plot(x,sol,color='k',lw=0.3)
        fnu = coeff * sol
        mdot = np.gradient(fnu,x) / data[3,1:]
        flux = mdot * omega * x**2 - fnu
        ax2.plot(x[ind],mdot[ind],color='k',lw=0.3)
        ax3.plot(x[ind],flux[ind],color='k',lw=0.3)

        
    ax1.set_xlabel(r'$R$',size=24,labelpad=0)
    ax1.set_ylabel(r'$\Sigma$',size=24,labelpad=0)
    ax1.set_xlim(1,100)
    ax1.set_ylim(0.001,0.05)
    ax1.tick_params(axis='both',labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax2.set_xlabel(r'$R$',size=24,labelpad=0)
    ax2.set_ylabel(r'$\dot{M}$',size=24,labelpad=0)
    ax2.set_xlim(1,x[ind].max())
    ax2.set_ylim(1.e-5,1e-3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    
    ax3.set_xlabel(r'$R$',size=24,labelpad=0)
    ax3.set_ylabel(r'$\dot{M}l - F_\nu$',size=24,labelpad=0)
    ax3.set_xlim(1,x[ind].max())
    ax3.set_ylim(-0.0001,0.0005)
    ax3.set_xscale('log')
    #ax3.set_yscale('log')
    

    if (ZAMFB):
        ax1.plot(x,x**(-0.5)*0.03,color='royalblue',zorder=0,lw=4.0,alpha=0.6)
        ax1.text(40,40**(-0.5)*0.03,r'$\propto R^{-1/2}$',
                transform=ax1.transData,size=18,color='royalblue')




        
    if (ZAMFB): 
        fig.savefig('diffusion_zamfb.png')
    else:
        fig.savefig('diffusion_zerotorque.png')
    plt.show()
        
