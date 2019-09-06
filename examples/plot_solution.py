
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')


path = './'


ZAMFB = True

filename = 'output.txt'


if __name__ == '__main__':

    data = np.loadtxt(path+filename)
    print(data.shape)
    x = data[0, 1:]

    fig = plt.figure(figsize=(18, 4))
    fig.subplots_adjust(bottom=0.15, right=0.99,
                        top=0.97, left=0.05, wspace=.25)

    ax0 = fig.add_subplot(141)
    ax1 = fig.add_subplot(142)
    ax2 = fig.add_subplot(143)
    ax3 = fig.add_subplot(144)

    omega = data[1, 1:]
    omegaprime = data[2, 1:]
    nu = data[4, 1:]
    coeff = -2 * np.pi * x**3 * nu * omegaprime
    fnu = coeff * data[8, 1:]
    mdot = (np.gradient(fnu, x) + 2 * np.pi * x**3 *
            (data[7, 1:] - data[6, 1:]) * data[8, 1:] * omega)/data[3, 1:]
    #mdot = (np.gradient(fnu,x))/data[3,1:]
    flux = mdot * omega * x**2 - fnu
    ind = mdot != 1e30

    print(x.shape, omega.shape)
    ax0.plot(x, omega, color='orange', lw=1.8)
    ax0.plot(x, -omegaprime, color='b', lw=1.8)
    ax0.plot(x, nu, color='g', lw=1.8)
    # ax0.plot(x[::20],-np.gradient(omega,x)[::20],'ks',ms=3.5)
    ax1.plot(x, data[8, 1:], color='r', lw=1.8)
    ax2.plot(x[ind], mdot[ind], color='r', lw=1.8)
    ax3.plot(x[ind], flux[ind], color='r', lw=1.8)

    for sol in data[9:, 1:]:
        ax1.plot(x, sol, color='k', lw=0.3)
        fnu = coeff * sol
        mdot = (np.gradient(fnu, x) + 2 * np.pi * x**3 *
                (data[7, 1:] - data[6, 1:]) * sol * omega)/data[3, 1:]
        #mdot = (np.gradient(fnu,x))/data[3,1:]
        flux = mdot * omega * x**2 - fnu
        ax2.plot(x[ind], mdot[ind], color='k', lw=0.3)
        ax3.plot(x[ind], flux[ind], color='k', lw=0.3)

    ax0.plot([2.8, 2.8], [1e-5, 1], 'gray', ls='--', zorder=0)
    ax0.set_xlabel(r'$R$', size=24, labelpad=0)
    # ax0.set_ylabel(r'$\Sigma$',size=24,labelpad=0)
    ax0.set_xlim(1, 100)
    # ax0.set_ylim(0.001,0.05)
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.text(0.15, 0.92, r'$\Omega$', transform=ax0.transAxes,
             size=18, color='orange')
    ax0.text(0.04, 0.72, r'$-\frac{d\Omega}{dR}$', transform=ax0.transAxes,
             size=24, color='b')
    ax0.text(0.3, 0.35, r'$\nu$', transform=ax0.transAxes,
             size=18, color='g')
    ax0.text(2.8, 2e-5, 'softening', transform=ax0.transData,
             size=12, color='gray', rotation=90, va='bottom', ha='right')

    ax1.set_xlabel(r'$R$', size=24, labelpad=0)
    ax1.set_ylabel(r'$\Sigma$', size=24, labelpad=0)
    ax1.set_xlim(1, 100)
    ax1.set_ylim(0.001, 0.05)
    ax1.tick_params(axis='both', labelsize=16)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax2.set_xlabel(r'$R$', size=24, labelpad=0)
    ax2.set_ylabel(r'$\dot{M}$', size=24, labelpad=0)
    ax2.set_xlim(1, x[ind].max())
    ax2.set_ylim(1.e-5, 1e-3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax3.set_xlabel(r'$R$', size=24, labelpad=0)
    ax3.set_ylabel(r'$\dot{M}l - F_\nu$', size=24, labelpad=0)
    ax3.set_xlim(1, x[ind].max())
    ax3.set_ylim(-0.0001, 0.0005)
    ax3.set_xscale('log')
    # ax3.set_yscale('log')

    if (ZAMFB):
        ax1.plot(x, x**(-0.5)*0.03, color='royalblue',
                 zorder=0, lw=4.0, alpha=0.6)
        ax1.text(30, 30**(-0.5)*0.03, r'$\propto R^{-1/2}$',
                 transform=ax1.transData, size=18, color='royalblue')

    if (ZAMFB):
        fig.savefig('diffusion_zamfb.png')
    else:
        fig.savefig('diffusion_zerotorque.png')
    plt.show()
