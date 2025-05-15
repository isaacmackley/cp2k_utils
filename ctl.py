import numpy as np
import matplotlib.pyplot as plt

colors = plt.cm.tab10.colors

def plot_ctl(defects,bulk,vbm,cbm,customx=False,xxmin=1,xxmax=4,customy=False,yymin=1,yymax=4):
    """
    Plots CTL diagrams for input defects

    defects: list of lists with format ['Defect Specie', Charge, Defect Cell Energy, Chemical Potential of added species, Err correction, Defect type ('i', 's', 'v' for interstial, substional or vacancy), Site name (e.g. O for an O site. Can be left as '' for interstial), Chemical potential for removed species)
    bulk: Energy of pristine cell
    vmb: vbm energy
    cbm: cbm energy

    note - all energies should be given in Ha
    """

    Ha_to_eV = 27.211386

    E_cbm = cbm * Ha_to_eV
    E_vbm = vbm * Ha_to_eV

    xmax=round((E_cbm-E_vbm)+1)
    #E_fermi = np.linspace(vbm-1,cbm+1,)
    E_fermi = np.linspace(-1,xmax,((xmax+1)*100)+1)

    
    types = {}
    n=0
    while n < len(defects):
        type_key = defects[n][0]+defects[n][5]+defects[n][6]

        if type_key not in types:
            types[type_key] = []
        types[type_key].append(defects[n])

        n=n+1
    
    E_low={}
    defect_types = []

    for type_key, typegroup in types.items():

        
    
        E_={}
        xequal0_={}
        xequalcbm_={}
        ymin1 = float(1000)
        ymin2 = float(1000)
        ymax = float(-1000)

        m=0

        

        for defect in typegroup:
            name = str(defect[0])+'_'+defect[5]+str(defect[1])
            E_[name] = defect[2]*Ha_to_eV - bulk*Ha_to_eV - defect[3]*Ha_to_eV + defect[7]*Ha_to_eV + defect[1] * (E_vbm + E_fermi) + defect[4]

            if typegroup[0][5] == 'i':
                defect_type = str(defect[0])+'_'+str(defect[5])
            if typegroup[0][5] == 's':
                defect_type = str(defect[0])+'_'+str(defect[6])
            if typegroup[0][5] == 'v':
                defect_type = 'V_'+str(defect[6])
            
            
            if m == 0:
                defect_types.append(defect_type)
                E_low[defect_type] = np.minimum(10000, E_[name])
            else:
                E_low[defect_type] = np.minimum(E_low[defect_type], E_[name])

            xequal0_[name] = defect[2]*Ha_to_eV - bulk*Ha_to_eV - defect[3]*Ha_to_eV + defect[7]*Ha_to_eV + defect[1] * (E_vbm) + defect[4]
            ymin1 = np.minimum(ymin1, xequal0_[name])
            ymax = np.maximum(ymax, xequal0_[name])

            xequalcbm_[name] = defect[2]*Ha_to_eV - bulk*Ha_to_eV - defect[3]*Ha_to_eV + defect[7]*Ha_to_eV + defect[1] * (E_cbm) + defect[4]
            ymin2 = np.minimum(ymin2, xequalcbm_[name])

            m=m+1

        ymin = np.minimum(ymin1, ymin2)
        ymin = round(ymin - 1)
        ymax = round(ymax + 1)

        plt.figure(figsize=(15,10))
        if typegroup[0][5] == 'i':
            plt.title(f'{typegroup[0][0]}$_{'i'}$ Charge Levels', size=40, pad=30)
        if typegroup[0][5] == 's':
            plt.title(f'{typegroup[0][0]}$_{{{typegroup[0][6]}}}$ Charge Levels', size=40, pad=30)
        if typegroup[0][5] == 'v':
            plt.title(f'V$_{{{typegroup[0][6]}}}$ Charge Levels', size=40, pad=30)
        plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30) 
        plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)

        #VBM and CBM 
        plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
        plt.axvline(x=E_cbm-E_vbm, color='tab:orange', linestyle='-', alpha=0.5)
        plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm-E_vbm,xmax,xmax,E_cbm-E_vbm], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

        #Defect Plots
        for defect in typegroup:

            name = str(defect[0])+'_'+defect[5]+str(defect[1])

            if defect[1] > 0:
                plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{str(abs(defect[1]))}{'+'}')
            elif defect[1] < 0:
                plt.plot(E_fermi, E_[name], marker = '', lw=3, label=f'{str(abs(defect[1]))}{'-'}')
            else:
                plt.plot(E_fermi, E_[name], marker = '', lw=3, label='neutral')



        #Lowest state line
        #plt.plot(E_fermi, x_int, linestyle='--', color='k', lw=2, label='')

        #Other details
        if customx == True:
            plt.xlim([xxmin,xxmax])
            plt.xticks(np.linspace(xxmin,xxmax,abs(xxmax-xxmin)+1),fontsize=20)
        if customx == False:
            plt.xlim([-1,xmax])
            plt.xticks(np.linspace(-1,xmax,abs(xmax)+2),fontsize=20)
        if customy == True:
            plt.ylim([yymin,yymax])
            plt.yticks(np.linspace(yymin,yymax,abs(yymax-yymin)+1),fontsize=20)
        if customy == False:
            plt.ylim([ymin,ymax])
            if ymin == 0 or ymax == 0:
                plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)
            else:
                if ymin/abs(ymin) == ymax/abs(ymax):
                    plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)-1),fontsize=20)
                else:
                    plt.yticks(np.linspace(ymin,ymax,abs(ymin)+abs(ymax)+1),fontsize=20)

        plt.grid(visible=0)#, which='major', axis='both',linestyle='--')

        plt.legend(markerscale=8.0, fontsize=20, loc='upper left')

        if typegroup[0][5] == 'i':
            plt.savefig(f'{typegroup[0][0]}_i CTL Diagram', bbox_inches='tight')
        if typegroup[0][5] == 's':
            plt.savefig(f'{typegroup[0][0]}_{typegroup[0][6]} CTL Diagram', bbox_inches='tight')
        if typegroup[0][5] == 'v':
            plt.savefig(f'V_{typegroup[0][6]} CTL Diagram', bbox_inches='tight')


    

    plt.figure(figsize=(15,10))
    plt.title('Charge Levels', size=40, pad=30)
    plt.xlabel('Fermi Energy (eV)', size=30, labelpad=30) 
    plt.ylabel('Formation Energy (eV)', size=30, labelpad=30)


    ymin=float(1000)
    ymax=float(-1000)
    n=0
    while n < len(defect_types):
        ymin = np.minimum(ymin, np.min(E_low[defect_types[n]]))
        ymax = np.maximum(ymax, np.max(E_low[defect_types[n]]))
        n=n+1
    
    ymin = round(ymin - 1)
    ymax = round(ymax + 1)
    #VBM and CBM 
    plt.axvline(x=0, color='tab:green', linestyle='-', alpha=0.5)
    plt.axvline(x=E_cbm-E_vbm, color='tab:orange', linestyle='-', alpha=0.5)
    plt.fill([-1,0,0,-1], [ymin,ymin,ymax,ymax], 'tab:green',[E_cbm-E_vbm,xmax,xmax,E_cbm-E_vbm], [ymin,ymin,ymax,ymax], 'tab:orange', alpha=0.2)

    n=0
    while n < len(defect_types):
        plt.plot(E_fermi, E_low[defect_types[n]], marker = '', lw=3, label=f'{defect_types[n]}')
        n=n+1

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

        plt.hlines(
            y=vbm_i,
            xmin=xmin_i,
            xmax=xmax_i,
            color='tab:green',
            linestyle='-',
            alpha=0.5
        )
        
        plt.hlines(
            y=cbm_i,
            xmin=xmin_i,
            xmax=xmax_i,
            color='tab:orange',
            linestyle='-',
            alpha=0.5
        )

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

        plt.hlines(
            y=vbm_i,
            xmin=xmin_i,
            xmax=xmax_i,
            color='tab:green',
            linestyle='-',
            alpha=0.5
        )
        
        plt.hlines(
            y=cbm_i,
            xmin=xmin_i,
            xmax=xmax_i,
            color='tab:orange',
            linestyle='-',
            alpha=0.5
        )

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