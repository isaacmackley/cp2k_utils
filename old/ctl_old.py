import numpy as np
import matplotlib.pyplot as plt


def plot_int_ctl(defects,bulk,vbm,cbm,ymin=0,ymax=8):
    """
    """

    Ha_to_eV = 27.211386

    E_cbm = cbm * Ha_to_eV
    E_vbm = vbm * Ha_to_eV

    xmax=round((E_cbm-E_vbm)+1)
    #E_fermi = np.linspace(vbm-1,cbm+1,)
    E_fermi = np.linspace(-1,xmax,(xmax+2))


    
    E_={}
    
    n=0
    while n < len(defects):
        name = str(defects[n][0])+'_'+defects[n][5]+str(defects[n][1])
        E_[name] = defects[n][2]*Ha_to_eV - bulk*Ha_to_eV - defects[n][3]*Ha_to_eV + defects[n][7] + defects[n][1] * (E_vbm + E_fermi) + defects[n][4]
        n=n+1


    plt.figure(figsize=(15,10))
    plt.title("Charge Levels", size=30, pad=30) 
    plt.xlabel('Fermi Energy (eV)', size=20, labelpad=30) 
    plt.ylabel('Formation Energy (eV)', size=20, labelpad=30)

    #VBM and CBM 
    plt.axvline(x=0, color='tab:orange', linestyle='-', alpha=0.5)
    plt.axvline(x=E_cbm-E_vbm, color='tab:blue', linestyle='-', alpha=0.5)
    plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:orange',[E_cbm-E_vbm,xmax,xmax,E_cbm-E_vbm], [ymin,ymin,ymax,ymax], 'tab:blue', alpha=0.2)

    #Defect Plots
    n=0
    while n < len(defects):

        name = str(defects[n][0])+'_'+defects[n][5]+str(defects[n][1])

        if defects[n][5] == 'i':
            if defects[n][1] > 0:
                if abs(defects[n][1]) == 1:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{defects[n][0]}$_{'i'}^{'+'}$')
                else:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{defects[n][0]}$_{'i'}^{{{str(abs(defects[n][1]))}}}{{^{'+'}}}$')
            elif defects[n][1] < 0:
                if abs(defects[n][1]) == 1:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{defects[n][0]}$_{'i'}^{'-'}$')
                else:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{defects[n][0]}$_{'i'}^{{{str(abs(defects[n][1]))}}}{{^{'-'}}}$')
            else:
                plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{defects[n][0]}$_{'i'}^{str(abs(defects[n][1]))}$')
        
        if defects[n][5] == 'v':
            if defects[n][1] > 0:
                if abs(defects[n][1]) == 1:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'${defects[n][6]}_{defects[n][0]}^{'+'}$')
                else:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'${defects[n][6]}_{defects[n][0]}^{{{str(abs(defects[n][1]))}}}{{^{'+'}}}$')
            elif defects[n][1] < 0:
                if abs(defects[n][1]) == 1:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'${defects[n][6]}_{defects[n][0]}^{'-'}$')
                else:
                    plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'${defects[n][6]}_{defects[n][0]}^{{{str(abs(defects[n][1]))}}}{{^{'-'}}}$')
            else:
                plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'${defects[n][6]}_{defects[n][0]}^{str(abs(defects[n][1]))}$')

        n=n+1


    #Lowest state line
    #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

    #Other details
    plt.xlim([-1,xmax])
    plt.ylim([ymin,ymax])

    plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

    plt.legend(markerscale=5.0, fontsize=15, loc='upper left')

    plt.savefig('CTLs')

    plt.show()

