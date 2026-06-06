import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dataclasses import dataclass

colors = plt.cm.tab10.colors
Ha_to_eV = 27.211386

@dataclass
class Defect:
    name: str
    charge: int
    energy: float       # Ha
    mu_added: float     # Ha
    correction: float   # eV
    defect_type: str    # 'i', 's', 'v'
    site: str
    mu_removed: float   # Ha


def plot_ctl(defects, bulk, vbm, cbm, xlim=None, ylim=None):
    """
    Plots CTL diagrams for input defects.

    defects: list of Defect dataclass instances
    bulk: energy of pristine cell (Ha)
    vbm: VBM energy (Ha)
    cbm: CBM energy (Ha)
    """
    
    Ha_to_eV = 27.211386

    E_cbm = cbm * Ha_to_eV
    E_vbm = vbm * Ha_to_eV

    xmax = round((E_cbm - E_vbm) + 1)
    E_fermi = np.linspace(-1, xmax, ((xmax + 1) * 100) + 1)

    type_labels = {
        'i': lambda d: f"{d.name}_i",
        's': lambda d: f"{d.name}_{d.site}",
        'v': lambda d: f"V_{d.site}",
    }

    types = {}
    for defect in defects:
        type_key = defect.name + defect.defect_type + defect.site
        if type_key not in types:
            types[type_key] = []
        types[type_key].append(defect)

    E_low = {}
    defect_types = []

    for type_key, typegroup in types.items():

        first = typegroup[0]
        defect_type = type_labels[first.defect_type](first)

        E_ = {}
        ymin1 = float('inf')
        ymin2 = float('inf')
        ymax_group = float('-inf')

        for i, defect in enumerate(typegroup):
            name = defect.name + '_' + defect.defect_type + str(defect.charge)
            base_E = (defect.energy * Ha_to_eV - bulk * Ha_to_eV - defect.mu_added * Ha_to_eV + defect.mu_removed * Ha_to_eV + defect.correction)

            E_[name] = base_E + defect.charge * (E_vbm + E_fermi)

            if i == 0:
                defect_types.append(defect_type)
                E_low[defect_type] = E_[name]
            else:
                E_low[defect_type] = np.minimum(E_low[defect_type], E_[name])

            e_at_vbm = base_E + defect.charge * E_vbm
            e_at_cbm = base_E + defect.charge * E_cbm

            ymin1 = min(ymin1, e_at_vbm)
            ymax_group = max(ymax_group, e_at_vbm)
            ymin2 = min(ymin2, e_at_cbm)

        ymin = round(min(ymin1, ymin2) - 1)
        ymax = round(ymax_group + 1)

        title = type_labels[first.defect_type](first)

        plt.figure(figsize=(15, 10))
        plt.title(title, size=40, pad=30)
        plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30)
        plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

        plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
        plt.axvline(x=E_cbm - E_vbm, color='tab:orange', linestyle='-', alpha=0.5)
        plt.fill([-1, 0, 0, -1], [ymin, ymin, ymax, ymax], 'tab:green', [E_cbm - E_vbm, xmax, xmax, E_cbm - E_vbm], [ymin, ymin, ymax, ymax], 'tab:orange', alpha=0.2)

        for defect in typegroup:
            name = defect.name + '_' + defect.defect_type + str(defect.charge)
            if defect.charge > 0:
                label = f'{abs(defect.charge)}+'
            elif defect.charge < 0:
                label = f'{abs(defect.charge)}-'
            else:
                label = 'neutral'
            plt.plot(E_fermi, E_[name], marker='', lw=3, label=label)

        if xlim is not None:
            plt.xlim(xlim)
        else:
            plt.xlim([-1, xmax])

        plt.xticks(fontsize=20)
        if ylim is not None:
            plt.ylim(ylim)
        else:
            plt.ylim([ymin, ymax])
            if ymin == 0 or ymax == 0:
                plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)
            else:
                if ymin / abs(ymin) == ymax / abs(ymax):
                    plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) - 1), fontsize=20)
                else:
                    plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)

        plt.grid(visible=0)
        plt.legend(markerscale=8.0, fontsize=20, loc='upper left')
        plt.savefig(f'{title} CTL Diagram', bbox_inches='tight')

    plt.figure(figsize=(15, 10))
    plt.title('Charge Levels', size=40, pad=30)
    plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30)
    plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

    ymin = float('inf')
    ymax = float('-inf')

    for defect_type in defect_types:
        ymin = min(ymin, np.min(E_low[defect_type]))
        ymax = max(ymax, np.max(E_low[defect_type]))

    ymin = round(ymin - 1)
    ymax = round(ymax + 1)

    plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
    plt.axvline(x=E_cbm - E_vbm, color='tab:orange', linestyle='-', alpha=0.5)
    plt.fill([-1, 0, 0, -1], [ymin, ymin, ymax, ymax], 'tab:green',
             [E_cbm - E_vbm, xmax, xmax, E_cbm - E_vbm], [ymin, ymin, ymax, ymax], 'tab:orange', alpha=0.2)

    for defect_type in defect_types:
        plt.plot(E_fermi, E_low[defect_type], marker='', lw=3, label=f'{defect_type}')

    plt.xlim([-1, xmax])
    plt.ylim([ymin, ymax])
    plt.xticks(np.linspace(-1, xmax, abs(xmax) + 2), fontsize=20)
    if ymin == 0 or ymax == 0:
        plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)
    else:
        if ymin / abs(ymin) == ymax / abs(ymax):
            plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) - 1), fontsize=20)
        else:
            plt.yticks(np.linspace(ymin, ymax, abs(ymin) + abs(ymax) + 1), fontsize=20)

    plt.legend(markerscale=5.0, fontsize=20, loc='upper left')
    plt.savefig('Charge Levels', bbox_inches='tight')



def plot_range_ctl(defects,bulks,legend=False):
    """
    Plots CTL range diagrams for input defects

    defects: list of lists with format ['Defect Specie', Charge, Defect Cell Energy, Chemical Potential of added species, Err correction, Defect type ('i', 's', 'v' for interstial, substional or vacancy), Site name (e.g. O for an O site. Can be left as '' for interstial), Chemical potential for removed species)
    bulk: Energy of pristine cell
    vmb: vbm energy
    cbm: cbm energy

    note - all energies should be given in Ha
    """

    Ha_to_eV = 27.211386

    E_vbm = 0
    E_cbm = 0
    n=0
    while n < len(bulks):
        E_vbm = E_vbm + (bulks[n][2] * Ha_to_eV)
        E_cbm = E_cbm + (bulks[n][3] * Ha_to_eV)
        n=n+1
    E_vbm_ave = E_vbm/len(bulks)
    E_cbm_ave = E_cbm/len(bulks)
    xmax=round((E_cbm_ave-E_vbm_ave)+1)
    E_fermi = np.linspace(-1,xmax,((xmax+1)*100)+1)


    types = {}
    n=0
    while n < len(defects):
        type_key = defects[n][2]+defects[n][7]+defects[n][8]

        if type_key not in types:
            types[type_key] = []
        types[type_key].append(defects[n])

        n=n+1
    
    E_low={}
    E_lo={}
    E_hi={}
    defect_types = []
    defect_charges = []

    for type_key, typegroup in types.items():

        
    
        E_={}
        xequal0_={}
        xequalcbm_={}
        ymin1 = float(1000)
        ymin2 = float(1000)
        ymax = float(-1000)

        m=0

        

        for defect in typegroup:   
            
            n=0
            while n<len(bulks):
                if bulks[n][0] == defect[0]:
                    bulk=bulks[n][1]
                    vbm=bulks[n][2]*Ha_to_eV
                    cbm=bulks[n][3]*Ha_to_eV
                    n=len(bulks)
                n=n+1
            
            name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
            E_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (vbm + E_fermi) + defect[6]
            
            if typegroup[0][7] == 'i':
                defect_type = str(defect[2])+'_'+str(defect[7])
            if typegroup[0][7] == 's':
                defect_type = str(defect[2])+'_'+str(defect[8])
            if typegroup[0][7] == 'v':
                defect_type = 'V_'+str(defect[8])
            
            
            defect_charge = defect_type+str(defect[3])
            
            if defect_charge not in defect_charges:
                defect_charges.append(defect_charge)
            if defect_charge not in E_lo:
                E_lo[defect_charge] = E_[name]
            else:
                E_lo[defect_charge] = np.minimum(E_lo[defect_charge], E_[name])
            if defect_charge not in E_hi:
                E_hi[defect_charge] = E_[name]
            else:
                E_hi[defect_charge] = np.maximum(E_hi[defect_charge], E_[name])

            if m == 0:
                defect_types.append(defect_type)
                E_low[defect_type] = E_[name]
            else:
                E_low[defect_type] = np.minimum(E_low[defect_type], E_[name])


            xequal0_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (vbm) + defect[6]
            ymin1 = np.minimum(ymin1, xequal0_[name])
            ymax = np.maximum(ymax, xequal0_[name])

            xequalcbm_[name] = defect[4]*Ha_to_eV - bulk*Ha_to_eV - defect[5]*Ha_to_eV + defect[9]*Ha_to_eV + defect[3] * (cbm) + defect[6]
            ymin2 = np.minimum(ymin2, xequalcbm_[name])

            m=m+1
           

        ymin = np.minimum(ymin1, ymin2)
        ymin = round(ymin - 1)
        ymax = round(ymax + 1)

        plt.figure(figsize=(15,10))
        if typegroup[0][7] == 'i':
            plt.title(f'{typegroup[0][2]}$_{'i'}$ Charge Levels', size=40, pad=30)
        if typegroup[0][7] == 's':
            plt.title(f'{typegroup[0][2]}$_{{{typegroup[0][8]}}}$ Charge Levels', size=40, pad=30)
        if typegroup[0][7] == 'v':
            plt.title(f'V$_{{{typegroup[0][8]}}}$ Charge Levels', size=40, pad=30)
        plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30) 
        plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

        #VBM and CBM 
        plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
        plt.axvline(x=E_cbm_ave-E_vbm_ave, color='tab:orange', linestyle='-', alpha=0.5)
        plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm_ave-E_vbm_ave,xmax,xmax,E_cbm_ave-E_vbm_ave], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

        #Defect Plots
        for defect in typegroup:
            
            name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
            
            q=-5
            while q < 6:
                color = colors[q % len(colors)]
                if defect[3] == q:
                    if defect[3] == 0:
                        plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label='neutral')
                    elif defect[3] > 0:
                        plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label=f'{str(abs(defect[3]))}{'+'}')
                    elif defect[3] < 0:
                        plt.plot(E_fermi, E_[name], marker = '', lw=3, color=color, label=f'{str(abs(defect[3]))}{'-'}')
                q=q+1


        #Lowest state line
        #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

        #Other details
        plt.xlim([-1,xmax])
        plt.ylim([ymin,ymax])
        plt.xticks(np.linspace(-1,xmax,abs(xmax)+2),fontsize=20)
        if ymin == 0 or ymax == 0:
            plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)
        else:
            if ymin/abs(ymin) == ymax/abs(ymax):
                plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)-1),fontsize=20)
            else:
                plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)

        plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

        #if legend==True:
            #plt.legend(markerscale=5.0, fontsize=15, loc='upper left')

        if typegroup[0][7] == 'i':
            plt.savefig(f'{typegroup[0][2]}_i Multi CTL Diagram', bbox_inches='tight')
        if typegroup[0][7] == 's':
            plt.savefig(f'{typegroup[0][2]}_{typegroup[0][8]} Multi CTL Diagram', bbox_inches='tight')
        if typegroup[0][7] == 'v':
            plt.savefig(f'V_{typegroup[0][8]} Multi CTL Diagram', bbox_inches='tight')
    



        plt.figure(figsize=(15,10))
        if typegroup[0][7] == 'i':
            plt.title(f'{typegroup[0][2]}$_{'i'}$ Charge Levels', size=30, pad=30)
        if typegroup[0][7] == 's':
            plt.title(f'{typegroup[0][2]}$_{{{typegroup[0][8]}}}$ Charge Levels', size=30, pad=30)
        if typegroup[0][7] == 'v':
            plt.title(f'V$_{{{typegroup[0][8]}}}$ Charge Levels', size=30, pad=30)
        plt.xlabel('Fermi Energy (eV)', size=20, labelpad=30) 
        plt.ylabel('Formation Energy (eV)', size=20, labelpad=30)

        #VBM and CBM 
        plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
        plt.axvline(x=E_cbm_ave-E_vbm_ave, color='tab:orange', linestyle='-', alpha=0.5)
        plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm_ave-E_vbm_ave,xmax,xmax,E_cbm_ave-E_vbm_ave], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

        #Defect Plots
        for defect in typegroup:

            name = 'cell'+str(defect[0])+'site'+str(defect[1])+defect[2]+'_'+defect[7]+str(defect[3])
            
            n=0
            while n < len(defect_charges):
                color = colors[n % len(colors)]
                plt.plot(E_fermi, E_lo[defect_charges[n]], marker = '', lw=1, color=color, label=f'{defect_charges[n]}')
                plt.plot(E_fermi, E_hi[defect_charges[n]], marker = '', lw=1, color=color)

                plt.fill([-1,xmax,xmax,-1], [np.interp(-1, E_fermi, E_lo[defect_charges[n]]),np.interp(xmax, E_fermi, E_lo[defect_charges[n]]),np.interp(xmax, E_fermi, E_hi[defect_charges[n]]),np.interp(-1, E_fermi, E_hi[defect_charges[n]])], color=color, alpha=0.2)

                n=n+1

        #Lowest state line
        #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

        #Other details
        plt.xlim([-1,xmax])
        plt.ylim([ymin,ymax])
        plt.xticks(np.linspace(-1,xmax,abs(xmax)+2))
        if ymin == 0 or ymax == 0:
            plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1))
        else:
            if ymin/abs(ymin) == ymax/abs(ymax):
                plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)-1))
            else:
                plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1))

        plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

        #if legend==True:
            #plt.legend(markerscale=5.0, fontsize=15, loc='upper left')

        if typegroup[0][7] == 'i':
            plt.savefig(f'{typegroup[0][2]}_i Range CTL Diagram', bbox_inches='tight')
        if typegroup[0][7] == 's':
            plt.savefig(f'{typegroup[0][2]}_{typegroup[0][8]} Range CTL Diagram', bbox_inches='tight')
        if typegroup[0][7] == 'v':
            plt.savefig(f'V_{typegroup[0][8]} Range CTL Diagram', bbox_inches='tight')
    




def plot_multi_mat_ctl(defects,materials,fermi_line=None):
    """
    Plots all given CTLs in given materials with band alignment.

    nb
        Ensure that the VBM and CBM energies are given relative to one another between different materials.
        Ensure that the 'Material Name' is the same in the defects and materials lists.
        CTL level energies to be given relative to the VBM of the material.

    args:

    defects:
            Each row: [Defect Name, Material Name, List of CTL tuples [Name, Energy] ]
            e.g.
                species = [
                    [f'N$_{'i'}$', f'HfO$_2$', [
                        ['+1/0', 1.3469378771661913],
                        ['0/-1', 2.8837208880316783],
                        ['-1/-2', 2.369073993248366],
                        ['-2/-3', 2.6696762876888496],
                    ]],
                    [f'N$_{'O'}$', f'HfO$_2$', [
                        ['3+/2+', 1.7950036300732979],
                        ['2+/1+', 3.0614261580351125],
                        ['+1/0', 2.0040361516655647],
                        ['0/-1', 2.164946988290245],
                    ]],
                    [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
                        ['0/-1', 3.380293837285865],
                        ['-1/-2', 4.198164515355801],
                        ['-2/-3', 4.441756959205842],
                        ['-3/-4', 3.7187798290658316]
                    ]]
                ]
    materials:
                Each row: [Material Name, VBM, CBM]
                e.g.
                    materials = [
                        ['Si', Sivbm, Sicbm],
                        [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
                    ]
    fermi_line:
                Draws dashed line at desired fermi level
    """

    defects_list = []
    materials_list = []

    for row in defects:
        defects_list.append(row[:2])
    for i in range(len(materials)):
        materials_list.append(materials[i][0])

    ymin = round(min(materials, key=lambda x: x[1])[1]-1)
    ymax = round(max(materials, key=lambda x: x[2])[2]+1)

    xmax = len(defects)+len(materials)

    plt.figure(figsize=(20,20))

    plt.ylim(ymin,ymax)
    plt.xlim(0,xmax)

    n=0
    for i in range(len(materials)):
        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        defects_list.insert(n+no_of_defects,['','',''])
        n=n+no_of_defects+1
    defects_list.insert(0,['','',''])

    plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=20)
    plt.yticks(fontsize=20)

    plt.ylabel('Fermi Energy (eV)', size=30, labelpad=30) 

    n=0
    for i in range(len(materials)):

        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        xmin_i = n
        xmax_i = n+no_of_defects+1
        vbm_i = materials[i][1]
        cbm_i = materials[i][2]

        # plt.hlines(
        #     y=vbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:green',
        #     linestyle='-',
        #     alpha=0.5
        # )
        
        # plt.hlines(
        #     y=cbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:orange',
        #     linestyle='-',
        #     alpha=0.5
        # )

        plt.fill(
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [ymin,vbm_i,vbm_i,ymin],
            'tab:green',
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [cbm_i,ymax,ymax,cbm_i],
            'tab:orange',
            alpha=0.2
        )

        x_i = n + (no_of_defects+1)/2

        plt.annotate(
            materials[i][0],
            [x_i,ymin+0.5],
            fontsize=20,
            ha='center'
        )
        
        n=n+no_of_defects+1
    
    if fermi_line is not None:
        plt.axhline(y=fermi_line, color='tab:cyan', linestyle='--', alpha=0.5)

    n=1
    for i in range(len(materials)):
        for row in defects:
            if row[1] == materials[i][0]:

                used_y = []
                for j in range(len(row[2])):
                    y_value = materials[i][1]+row[2][j][1]
                    offset = 0
                    for used in used_y:
                        if abs(y_value-used) < 0.1:
                            offset += 0.1
                    
                    plt.hlines(y_value, n-0.1, n+0.1, lw=3)
                    plt.annotate(
                        f'{row[2][j][0]}',
                        [n,y_value+offset],
                        va='center',
                        xytext=(30,0),
                        textcoords='offset points',
                        fontsize=20
                    )
                    used_y.append(y_value+offset)
                n=n+1
        n=n+1
            
    plt.savefig('Multi_Mat_CTLs.png')
    plt.show()









def plot_multi_mat_ctl_errors(defects,materials,fermi_line=None):
    """
    Plots all given CTLs in given materials with band alignment.

    nb
        Ensure that the VBM and CBM energies are given relative to one another between different materials.
        Ensure that the 'Material Name' is the same in the defects and materials lists.
        CTL level energies to be given relative to the VBM of the material.

    args:

    defects:
            Each row: [Defect Name, Material Name, List of CTL tuplesenergies [Name, Ave_E, Min_E, Max_E] ]
            e.g.
                species = [
                    [f'N$_{'i'}$', f'HfO$_2$', [
                        ['+1/0', 1.3469378771661913],
                        ['0/-1', 2.8837208880316783],
                        ['-1/-2', 2.369073993248366],
                        ['-2/-3', 2.6696762876888496],
                    ]],
                    [f'N$_{'O'}$', f'HfO$_2$', [
                        ['3+/2+', 1.7950036300732979],
                        ['2+/1+', 3.0614261580351125],
                        ['+1/0', 2.0040361516655647],
                        ['0/-1', 2.164946988290245],
                    ]],
                    [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
                        ['0/-1', 3.380293837285865],
                        ['-1/-2', 4.198164515355801],
                        ['-2/-3', 4.441756959205842],
                        ['-3/-4', 3.7187798290658316]
                    ]]
                ]
    materials:
                Each row: [Material Name, VBM, CBM]
                e.g.
                    materials = [
                        ['Si', Sivbm, Sicbm],
                        [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
                    ]
    fermi_line:
                Draws dashed line at desired fermi level
    """

    #color_list = ["tab:orange", "tab:pink", "tab:blue", "tab:brown", "tab:red", "tab:green"]

    color_list = {
        '3+/2+': 'tab:red',
        '2+/1+': 'tab:olive',
        '1+/0': 'tab:cyan',
        '0/1-': 'tab:blue',
        '1-/2-': 'tab:pink',
        '2-/3-': 'tab:purple',
        '3-/4-': 'tab:brown'
    }

    defects_list = []
    materials_list = []

    for row in defects:
        defects_list.append(row[:2])
    for i in range(len(materials)):
        materials_list.append(materials[i][0])

    ymin = round(min(materials, key=lambda x: x[1])[1]-1)
    ymax = round(max(materials, key=lambda x: x[2])[2]+1)

    xmax = len(defects)+len(materials)

    plt.figure(figsize=(30,20))
    plt.title('Charge Transition Levels', size=60, pad=60)

    plt.ylim(ymin,ymax)
    plt.xlim(0,xmax)

    n=0
    for i in range(len(materials)):
        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        defects_list.insert(n+no_of_defects,['','',''])
        n=n+no_of_defects+1
    defects_list.insert(0,['','',''])

    plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=30)
    plt.yticks(fontsize=30)

    plt.ylabel('Fermi Energy (eV)', size=40, labelpad=20)
    plt.xlabel('Defects', size=40, labelpad=30) 

    n=0
    for i in range(len(materials)):

        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        xmin_i = n
        xmax_i = n+no_of_defects+1
        vbm_i = materials[i][1]
        cbm_i = materials[i][2]

        # plt.hlines(
        #     y=vbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:green',
        #     linestyle='-',
        #     alpha=0.5
        # )
        
        # plt.hlines(
        #     y=cbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:orange',
        #     linestyle='-',
        #     alpha=0.5
        # )

        plt.fill(
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [ymin,vbm_i,vbm_i,ymin],
            'tab:green',
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [cbm_i,ymax,ymax,cbm_i],
            'tab:orange',
            alpha=0.2
        )

        x_i = n + (no_of_defects+1)/2

        plt.annotate(
            materials[i][0],
            [x_i,ymin+0.5],
            fontsize=30,
            ha='center'
        )
        
        n=n+no_of_defects+1
    
    if fermi_line is not None:
        plt.axhline(y=fermi_line, color='tab:cyan', linestyle='--', alpha=0.5)

    n=1
    for i in range(len(materials)):
        for row in defects:
            if row[1] == materials[i][0]:

                used_y = []

                for j in range(len(row[2])):
                    w = 0.4/(len(row[2])-1)
                    m = n + w * (j - (len(row[2]) - 1) / 2)

                    min_y = materials[i][1]+row[2][j][2]
                    max_y = materials[i][1]+row[2][j][3]
                    plt.fill(
                        [m-0.05,m+0.05,m+0.05,m-0.05],
                        [min_y,min_y, max_y,max_y],
                        color_list.get(f'{row[2][j][0]}'),
                        alpha=0.2
                    )

                for j in range(len(row[2])):
                    w = 0.4/(len(row[2])-1)
                    m = n + w * (j - (len(row[2]) - 1) / 2)

                    min_y = materials[i][1]+row[2][j][2]
                    max_y = materials[i][1]+row[2][j][3]
                    y_value = materials[i][1]+row[2][j][1]
                    offset = 0
                    for used in used_y:
                        if abs(y_value-used) < 0.1:
                            offset += 0.1
                    
                    o = n-0.4 if m<n else n+0.4 if m>n else n+0.4
                    plt.hlines(min_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    plt.hlines(max_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    plt.hlines(y_value, m-0.05, m+0.05, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    # plt.annotate(
                    #     f'{row[2][j][0]}',
                    #     [o,y_value+offset],
                    #     va='center',
                    #     ha='center',
                    #     xytext=(0,0),
                    #     textcoords='offset points',
                    #     fontsize=20,
                    #     color=color_list.get(f'{row[2][j][0]}')
                    # )
                    used_y.append(y_value+offset)
                n=n+1
        n=n+1
    
    custom_leg = [[],[]]

    for key in color_list:
        custom_leg[0].append(key)
        custom_leg[1].append(Line2D([0], [0], color=color_list[key], lw=6))
    
    
    plt.legend(custom_leg[1], custom_leg[0], loc='center right', fontsize=30)

    plt.savefig('Multi_Mat_CTLs_error.png', bbox_inches='tight')
    plt.show()






def temp_plot_multi_mat_ctl_errors(defects,materials,fermis):
    """
    Plots all given CTLs in given materials with band alignment.

    nb
        Ensure that the VBM and CBM energies are given relative to one another between different materials.
        Ensure that the 'Material Name' is the same in the defects and materials lists.
        CTL level energies to be given relative to the VBM of the material.

    args:

    defects:
            Each row: [Defect Name, Material Name, List of CTL tuplesenergies [Name, Ave_E, Min_E, Max_E] ]
            e.g.
                species = [
                    [f'N$_{'i'}$', f'HfO$_2$', [
                        ['+1/0', 1.3469378771661913],
                        ['0/-1', 2.8837208880316783],
                        ['-1/-2', 2.369073993248366],
                        ['-2/-3', 2.6696762876888496],
                    ]],
                    [f'N$_{'O'}$', f'HfO$_2$', [
                        ['3+/2+', 1.7950036300732979],
                        ['2+/1+', 3.0614261580351125],
                        ['+1/0', 2.0040361516655647],
                        ['0/-1', 2.164946988290245],
                    ]],
                    [f'N$_{'2'}$$_{'i'}$', f'HfO$_2$', [
                        ['0/-1', 3.380293837285865],
                        ['-1/-2', 4.198164515355801],
                        ['-2/-3', 4.441756959205842],
                        ['-3/-4', 3.7187798290658316]
                    ]]
                ]
    materials:
                Each row: [Material Name, VBM, CBM]
                e.g.
                    materials = [
                        ['Si', Sivbm, Sicbm],
                        [f'HfO$_{'2'}$', HfO2vbm, HfO2cbm]
                    ]
    fermi_line:
                Draws dashed line at desired fermi level
    """

    #color_list = ["tab:orange", "tab:pink", "tab:blue", "tab:brown", "tab:red", "tab:green"]

    color_list = {
        #'3+/2+': 'tab:red',
        '2+/1+': 'tab:olive',
        '1+/0': 'tab:cyan',
        '0/1-': 'tab:blue',
        '1-/2-': 'tab:pink',
        #'2-/3-': 'tab:purple',
        #'3-/4-': 'tab:brown'
    }

    defects_list = []
    materials_list = []

    for row in defects:
        defects_list.append(row[:2])
    for i in range(len(materials)):
        materials_list.append(materials[i][0])

    ymin = round(min(materials, key=lambda x: x[1])[1]-1)
    ymax = round(max(materials, key=lambda x: x[2])[2]+1)

    xmax = len(defects)+len(materials)

    plt.figure(figsize=(20,20))
    #plt.title('Charge Transition Levels', size=60, pad=60)

    plt.ylim(ymin,ymax)
    plt.xlim(0,xmax)

    n=0
    for i in range(len(materials)):
        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        defects_list.insert(n+no_of_defects,['','',''])
        n=n+no_of_defects+1
    defects_list.insert(0,['','',''])

    plt.xticks(np.linspace(0,xmax,xmax+1), [row[0] for row in defects_list], fontsize=30)
    plt.yticks(fontsize=30)

    plt.ylabel('Fermi Energy (eV)', size=40, labelpad=20)
    #plt.xlabel('Defects', size=40, labelpad=30) 

    n=0
    for i in range(len(materials)):

        no_of_defects = sum(1 for row in defects if row[1] == materials[i][0])
        xmin_i = n
        xmax_i = n+no_of_defects+1
        vbm_i = materials[i][1]
        cbm_i = materials[i][2]

        # plt.hlines(
        #     y=vbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:green',
        #     linestyle='-',
        #     alpha=0.5
        # )
        
        # plt.hlines(
        #     y=cbm_i,
        #     xmin=xmin_i,
        #     xmax=xmax_i,
        #     color='tab:orange',
        #     linestyle='-',
        #     alpha=0.5
        # )

        plt.fill(
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [ymin,vbm_i,vbm_i,ymin],
            'tab:green',
            [xmin_i,xmin_i,xmax_i,xmax_i],
            [cbm_i,ymax,ymax,cbm_i],
            'tab:orange',
            alpha=0.2
        )

        x_i = n + (no_of_defects+1)/2

        plt.annotate(
            materials[i][0],
            [x_i,ymin+0.5],
            fontsize=30,
            ha='center'
        )
        
        n=n+no_of_defects+1
    
    for n in range(len(fermis)):
        plt.axhline(y=fermis[n][1], label=fermis[n][0], linestyle='--', lw=3, alpha=0.8)

    n=1
    for i in range(len(materials)):
        for row in defects:
            if row[1] == materials[i][0]:

                used_y = []

                for j in range(len(row[2])):
                    w = 0.4/(len(row[2])-1)
                    m = n + w * (j - (len(row[2]) - 1) / 2)

                    min_y = materials[i][1]+row[2][j][2]
                    max_y = materials[i][1]+row[2][j][3]
                    plt.fill(
                        [m-0.05,m+0.05,m+0.05,m-0.05],
                        [min_y,min_y, max_y,max_y],
                        color_list.get(f'{row[2][j][0]}'),
                        alpha=0.2
                    )

                for j in range(len(row[2])):
                    w = 0.4/(len(row[2])-1)
                    m = n + w * (j - (len(row[2]) - 1) / 2)

                    min_y = materials[i][1]+row[2][j][2]
                    max_y = materials[i][1]+row[2][j][3]
                    y_value = materials[i][1]+row[2][j][1]
                    offset = 0
                    for used in used_y:
                        if abs(y_value-used) < 0.1:
                            offset += 0.1
                    
                    o = n-0.4 if m<n else n+0.4 if m>n else n+0.4
                    plt.hlines(min_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    plt.hlines(max_y, m-0.075, m+0.075, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    plt.hlines(y_value, m-0.05, m+0.05, lw=3, colors=color_list.get(f'{row[2][j][0]}'))
                    # plt.annotate(
                    #     f'{row[2][j][0]}',
                    #     [o,y_value+offset],
                    #     va='center',
                    #     ha='center',
                    #     xytext=(0,0),
                    #     textcoords='offset points',
                    #     fontsize=20,
                    #     color=color_list.get(f'{row[2][j][0]}')
                    # )
                    used_y.append(y_value+offset)
                n=n+1
        n=n+1
    
    # custom_leg = [[],[]]

    # for key in color_list:
    #     custom_leg[0].append(key)
    #     custom_leg[1].append(Line2D([0], [0], color=color_list[key], lw=6))
    
    
    #plt.legend(custom_leg[1], custom_leg[0], loc='center right', fontsize=30)
    #plt.legend(fontsize=20)
    plt.savefig('Mats_fermis.png', bbox_inches='tight')
    plt.show()