import re
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt
from cp2k_utils.Elements import ELEMENT_COLOURS
from cp2k_utils.tools import lighten_hex_colour

HEADER_MATCH = re.compile(
    r'\# Projected DOS for atomic kind (?P<element>\w+) at iteration step i = \d+, E\(Fermi\) = [ \t]* (?P<Efermi>[^\t ]+) a\.u\.')

EIGENVALUE_COLUMN = 1
DENSITY_COLUMN = 3


def smeared_pdos(pdosfilenames, sigma=0.02, de=0.001, scale=1, total_sum=False, no_header=False, output=None):
    alldata = []
    orb_headers = []

    for pdosfilename in pdosfilenames:
        with open(pdosfilename, 'r') as fhandle:
            match = HEADER_MATCH.match(fhandle.readline().rstrip())
            if not match:
                raise ValueError(f"The file '{pdosfilename}' does not look like a CP2K PDOS output.")
            efermi = float(match.group('Efermi'))
            header = fhandle.readline().rstrip().split()[1:]
            header[1:3] = [' '.join(header[1:3])]
            data = np.loadtxt(fhandle)
        alldata.append(data)
        orb_headers += header[DENSITY_COLUMN:]

    margin = 10 * sigma
    emin = min(np.min(data[:, EIGENVALUE_COLUMN]) for data in alldata) - margin
    emax = max(np.max(data[:, EIGENVALUE_COLUMN]) for data in alldata) + margin
    ncols = sum(data.shape[1] - DENSITY_COLUMN for data in alldata)
    nmesh = int((emax - emin) / de) + 1
    xmesh = np.linspace(emin, emax, nmesh)
    ymesh = np.zeros((nmesh, ncols))

    fact = de / (sigma * np.sqrt(2.0 * np.pi))

    coloffset = 0
    for data in alldata:
        ncol = data.shape[1] - DENSITY_COLUMN
        for idx in range(nmesh):
            func = np.exp(-(xmesh[idx] - data[:, EIGENVALUE_COLUMN]) ** 2 / (2.0 * sigma ** 2)) * fact
            ymesh[idx, coloffset:(coloffset + ncol)] = func.dot(data[:, DENSITY_COLUMN:])
        coloffset += ncol

    if total_sum:
        finalsum = np.sum(ymesh, 0) * de
        print("Sum over all meshpoints, per orbital:")
        print(("{:16.8f}" * ncols).format(*finalsum))

    xmesh -= efermi
    xmesh *= 27.211384
    ymesh *= scale

    result = []
    if not no_header:
        header_row = ["Energy_[eV]"] + orb_headers
        result.append(header_row)
    for idx in range(nmesh):
        data_row = [xmesh[idx]] + list(ymesh[idx, :])
        result.append(data_row)

    if output:
        with open(output, 'w') as f:
            if not no_header:
                f.write("{:>16}".format("Energy_[eV]") + "".join("{:>16}".format(header) for header in orb_headers) + "\n")
            for row in result[1:]:
                f.write("{:16.8f}".format(row[0]) + "".join("{:16.8f}".format(val) for val in row[1:]) + "\n")

    return result


def pdos_all(elements,spin=True,sigma=0.003,output=True):
    """
    
    """

        
    k=0

    for element in elements:
        k=k+1
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*ALPHA_k'+str(k)+'*'):
                alpha=file
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*BETA_k'+str(k)+'*'):
                beta=file
            

        element_name="element"+str(k)
        names={}
        names[element_name] = smeared_pdos([alpha,beta],sigma)

        if output:
            smeared_pdos([alpha,beta],sigma,output=f'{element}_pdos')
            







def pdos_plot(elements,spin=True,sigma=0.003):
    """
    
    """

        
    k=0
    names={}
    plots={}


    for element in elements:
        k=k+1
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*ALPHA_k'+str(k)+'*'):
                alpha=file
        for file in os.listdir('.'):
            if fnmatch.fnmatch(file, '*BETA_k'+str(k)+'*'):
                beta=file

        element_no=f"element{k}"
        
        names[element_no] = smeared_pdos([alpha,beta],sigma)


    n=1
    while n <= len(elements):
    
        x = []
        y = []

        s = []
        p = [] 
        d = [] 
        f = []

        y_2 = []

        s_2 = []
        p_2 = [] 
        d_2 = [] 
        f_2 = []

        i=1

        orbs = (len(names[f'element{n}'][0])-1)/2

        while i < len(names[f'element{n}']):

            x.append(names[f'element{n}'][i][0])
            
            s.append(names[f'element{n}'][i][1])

            if orbs == 1:
                s_2.append(names[f'element{n}'][i][2]*-1)
                y = [s]
                y_2 = [s_2]

            if orbs == 2:
                p.append(names[f'element{n}'][i][2])
                s_2.append(names[f'element{n}'][i][3]*-1)
                p_2.append(names[f'element{n}'][i][4]*-1)
                y = [s, p]
                y_2 = [s_2, p_2]

            if orbs == 3:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                s_2.append(names[f'element{n}'][i][4]*-1)
                p_2.append(names[f'element{n}'][i][5]*-1)
                d_2.append(names[f'element{n}'][i][6]*-1)
                y = [s, p, d]
                y_2 = [s_2, p_2, d_2]

            if orbs == 4:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                f.append(names[f'element{n}'][i][4])
                s_2.append(names[f'element{n}'][i][5]*-1)
                p_2.append(names[f'element{n}'][i][6]*-1)
                d_2.append(names[f'element{n}'][i][7]*-1)
                f_2.append(names[f'element{n}'][i][8]*-1)
                y = [s, p, d, f]
                y_2 = [s_2, p_2, d_2, f_2]
                        

            i=i+1

        element_name = f'{elements[n-1]}'
        plots[f'{element_name}_x'] = x
        plots[f'{element_name}_alpha_tot'] = sum(map(np.array, y))
        plots[f'{element_name}_beta_tot']  = sum(map(np.array, y_2))


        plt.figure(figsize=(25,10))
        plt.title(f"{elements[n-1]} Density of States", size=30) 
        plt.xlabel('Energy (eV)', size=20) 
        plt.ylabel('Denisty of States (arb.)', size=20)

        xmin=-3
        xmax=6

        plt.plot(x, s, marker = '',c=ELEMENT_COLOURS[elements[n-1]],linestyle='solid', label=f'{elements[n-1]}_s')
        plt.plot(x, s_2, marker = '',c=ELEMENT_COLOURS[elements[n-1]],linestyle='solid')
        
        ymax = max(max([y for x, y in zip(x, s) if xmin <= x <= xmax]), max([y for x, y in zip(x, s_2) if xmin <= x <= xmax], key=abs))
        if orbs >= 2:
            plt.plot(x, p, marker = '', label=f'{elements[n-1]}_p',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=0.5),linestyle='solid')
            plt.plot(x, p_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=0.5),linestyle='solid')
            ymax = max(ymax, max([y for x, y in zip(x, p) if xmin <= x <= xmax]), max([y for x, y in zip(x, p_2) if xmin <= x <= xmax], key=abs))
            if orbs >= 3:
                plt.plot(x, d, marker = '', label=f'{elements[n-1]}_d',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.2),linestyle='solid')
                plt.plot(x, d_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.2),linestyle='solid')
                ymax = max(ymax, max([y for x, y in zip(x, d) if xmin <= x <= xmax]), max([y for x, y in zip(x, d_2) if xmin <= x <= xmax], key=abs))
                if orbs >= 4:
                    plt.plot(x, f, marker = '', label=f'{elements[n-1]}_f',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.5),linestyle='solid')
                    plt.plot(x, f_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.5),linestyle='solid')
                    ymax = max(ymax, max([y for x, y in zip(x, f) if xmin <= x <= xmax]), max([y for x, y in zip(x, f_2) if xmin <= x <= xmax], key=abs))

        plt.axvline(x=0, color='k', linestyle='--')

        plt.grid(visible=1)

        ymax = np.ceil(ymax*10)/10

        if ymax > 1.0:
            plt.yticks(np.linspace(-ymax,ymax,int(ymax*10)+1))
        if 0.5 < ymax <= 1.0:
            plt.yticks(np.linspace(-ymax,ymax,int(ymax*20)+1))
        if 0.2 < ymax <= 0.5:
            plt.yticks(np.linspace(-ymax,ymax,int(ymax*40)+1))
        if ymax <= 0.2:
            plt.yticks(np.linspace(-ymax,ymax,int(ymax*80)+1))
        
        plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax)*4)+1))
        
        plt.xlim([xmin,xmax])
        plt.ylim([-ymax,ymax])

        plt.legend(markerscale=10.0, fontsize=20)

        plt.savefig(f'{elements[n-1]} Density of States')

        plt.show()

        n=n+1

    
    plt.figure(figsize=(25,10))
    plt.title("Total Density of States", size=30) 
    plt.xlabel('Energy (eV)', size=20) 
    plt.ylabel('Denisty of States (arb.)', size=20)

    xmin=-3
    xmax=6

    ymax = max(max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_beta_tot']) if xmin <= x <= xmax], key=abs))
    
    j=1
    while j <= len(elements):
        ymax = max(ymax,max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot']) if xmin <= x <= xmax], key=abs))
        plt.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot'], marker = '',c=ELEMENT_COLOURS[elements[j-1]], label=f'{elements[j-1]}')
        plt.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot'],c=ELEMENT_COLOURS[elements[j-1]], marker = '')
        j=j+1

    plt.axvline(x=0, color='k', linestyle='--')

    plt.grid(visible=1)

    ymax = np.ceil(ymax*10)/10
    
    if ymax > 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(ymax*10)+1))
    if 0.5 < ymax <= 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(ymax*20)+1))
    if 0.2 < ymax <= 0.5:
        plt.yticks(np.linspace(-ymax,ymax,int(ymax*40)+1))
    if ymax <= 0.2:
        plt.yticks(np.linspace(-ymax,ymax,int(ymax*80)+1))
        
    plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax)*4)+1))

    plt.xlim([xmin,xmax])
    plt.ylim([-ymax,ymax])

    plt.legend(markerscale=10.0, fontsize=20)

    plt.savefig('Total Density of States')

    plt.show()
