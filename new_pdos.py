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



class ElementPDOS:

    def __init__(self,name,index,spin,sigma):
        self.name = name
        self.index = index
        self.spin = spin
        self.sigma = sigma
        self.alpha_file = None
        self.beta_file = None
        self.data = self.load_pdos_data()
    
    def find_files(self):
        files = os.listdir('.')
        k = self.index

        if self.spin:
            self.alpha_file = next(f for f in files if fnmatch.fnmatch(f, f'*ALPHA_k{k}*'))
            self.beta_file = next(f for f in files if fnmatch.fnmatch(f, f'*BETA_k{k}*'))
        else:
            matches = [f for f in files if fnmatch.fnmatch(f, f'*-k{k}*.pdos') or fnmatch.fnmatch(f, f'*k{k}*.pdos')]
            if matches:
                self.alpha_file = matches[0]

    def load_pdos_data(self):
        self.find_files()
        if self.spin:
            return smeared_pdos([self.alpha_file, self.beta_file], self.sigma)
        else:
            return smeared_pdos([self.alpha_file], self.sigma)
    
    def smear_data(self):

        self.find_files()

        with open(self.alpha_file, 'r') as alpha:
            data = np.loadtxt(alpha)
            header = alpha.readline().rstrip()

        f = 1 / (self.sigma * np.sqrt(2*np.pi))

        f * np.exp( - (E - E_i)**2 / (2 * self.sigma**2) )
    



for pdosfilename in pdosfilenames:
        with open(pdosfilename, 'r') as fhandle:
            match = HEADER_MATCH.match(fhandle.readline().rstrip())
            if not match:
                raise ValueError(f"The file '{pdosfilename}' does not look like a CP2K PDOS output.")
            fermi = float(match.group('Efermi'))
            header = fhandle.readline().rstrip().split()[1:]
            header[1:3] = [' '.join(header[1:3])]
            data = np.loadtxt(fhandle)
        alldata.append(data)
        orb_headers += header[DENSITY_COLUMN:]
        efermi = max(efermi, fermi)





    def extract_orbitals(self):
        x_vals = []
        alpha_orbitals = []
        beta_orbitals = [] if self.spin else None

        n_orbitals = (len(self.data[0]) - 1) // 2 if self.spin else len(self.data[0]) - 1

        for _ in range(n_orbitals):
            alpha_orbitals.append([])
            if self.spin:
                beta_orbitals.append([])

        for row in self.data[1:]:
            x_vals.append(row[0])
            for i in range(n_orbitals):
                alpha_orbitals[i].append(row[i + 1])
                if self.spin:
                    beta_orbitals[i].append(-row[i + 1 + n_orbitals])

        return x_vals, alpha_orbitals, beta_orbitals
    
    def get_ymax(self, xmin, xmax):
        x_vals, alpha_orbitals, beta_orbitals = self.extract_orbitals()

        ymax = 0

        def max_in_range(data):
            return max([abs(y) for x, y in zip(x_vals, data) if xmin <= x <= xmax], default=0)

        for alpha in alpha_orbitals:
            ymax = max(ymax, max_in_range(alpha))

        if self.spin and beta_orbitals:
            for beta in beta_orbitals:
                ymax = max(ymax, max_in_range(beta))

        return ymax




class PDOSPlotter:

    def __init__(self, elements, spin=True, grid=False, sigma=0.003, vbm=0, xmin=-3, xmax=6, figw=25, figh=10):
        self.elements = [ElementPDOS(name, i + 1, spin, sigma) for i, name in enumerate(elements)]
        self.spin = spin
        self.grid = grid
        self.vbm = vbm
        self.xmin = xmin
        self.xmax = xmax
        self.figw = figw
        self.figh = figh


    def plot_individual(self):

        for element in self.elements:
            x, alpha, beta = element.extract_orbitals()
            plt.figure(figsize=(self.figw, self.figh))
            plt.title(f"{element.name} Density of States", size=40)
            plt.xlabel('Energy (eV)', size=30) 
            plt.ylabel('Denisty of States (arb.)', size=30)

            plt.axvline(x=self.vbm, color='k', linestyle='--', alpha=0.5)
            plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
            plt.axhline(y=0, color='k', linestyle='-', alpha=0.5)

            orbital_labels = ['s', 'p', 'd', 'f']
            colors = [ELEMENT_COLOURS[element.name],  # s
                      lighten_hex_colour(ELEMENT_COLOURS[element.name], 0.5),  # p
                      lighten_hex_colour(ELEMENT_COLOURS[element.name], -0.2), # d
                      lighten_hex_colour(ELEMENT_COLOURS[element.name], -0.5)]  # f

            for i, (label, alpha_vals) in enumerate(zip(orbital_labels, alpha)):
                plt.plot(x, alpha_vals, label=f"{element.name}_{label}", color=colors[i])
                if self.spin:
                    plt.plot(x, beta[i], color=colors[i])

            
            plt.xlim([self.xmin, self.xmax])
            plt.xticks(np.linspace(self.xmin,self.xmax,int((abs(self.xmin)+self.xmax))+1),fontsize=20)

            ymax = element.get_ymax(self.xmin, self.xmax)
            ymax = np.ceil(ymax * 10) / 10 if ymax > 0.1 else np.ceil(ymax * 100) / 100  
            plt.ylim([-ymax, ymax] if self.spin else [0, ymax])         
            plt.yticks([])

            if self.grid:
                plt.grid(visible=1)

            plt.legend(markerscale=10.0, fontsize=20)
            plt.savefig(f"{element.name}_DOS.png", bbox_inches='tight')

            plt.show()

    def plot_total(self):
        
        plt.figure(figsize=(self.figw, self.figh))
        plt.title("Density of States", size=40)
        plt.xlabel('Energy (eV)', size=30) 
        plt.ylabel('Denisty of States (arb.)', size=30)

        plt.axvline(x=self.vbm, color='k', linestyle='--', alpha=0.5)
        plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.5)

        ymax = []

        for element in self.elements:
            x, alpha, beta = element.extract_orbitals()

            color = ELEMENT_COLOURS[element.name]

            alpha_vals = np.sum(alpha, axis=0)         
            beta_vals = np.sum(beta, axis=0)

            plt.plot(x, alpha_vals, label=f"{element.name}", color=color)
            if self.spin:
                plt.plot(x, beta_vals, color=color)
            
            ymax = np.append(ymax, element.get_ymax(self.xmin, self.xmax))

        plt.xlim([self.xmin, self.xmax])
        plt.xticks(np.linspace(self.xmin,self.xmax,int((abs(self.xmin)+self.xmax))+1),fontsize=20)

        ymax = max(ymax)
        ymax = np.ceil(ymax * 10) / 10 if ymax > 0.1 else np.ceil(ymax * 100) / 100  
        plt.ylim([-ymax, ymax] if self.spin else [0, ymax])         
        plt.yticks([])

        if self.grid:
            plt.grid(visible=1)

        plt.legend(markerscale=10.0, fontsize=20)
        plt.savefig(f"DOS.png", bbox_inches='tight')

        plt.show()










HEADER_MATCH = re.compile(
    r'\# Projected DOS for atomic kind (?P<element>\w+) at iteration step i = \d+, E\(Fermi\) = [ \t]* (?P<Efermi>[^\t ]+) a\.u\.')

EIGENVALUE_COLUMN = 1
DENSITY_COLUMN = 3


def smeared_pdos(pdosfilenames, sigma=0.02, de=0.001, scale=1, total_sum=False, no_header=False, output=None):
    alldata = []
    orb_headers = []
    efermi = -1000

    for pdosfilename in pdosfilenames:
        with open(pdosfilename, 'r') as fhandle:
            match = HEADER_MATCH.match(fhandle.readline().rstrip())
            if not match:
                raise ValueError(f"The file '{pdosfilename}' does not look like a CP2K PDOS output.")
            fermi = float(match.group('Efermi'))
            header = fhandle.readline().rstrip().split()[1:]
            header[1:3] = [' '.join(header[1:3])]
            data = np.loadtxt(fhandle)
        alldata.append(data)
        orb_headers += header[DENSITY_COLUMN:]
        efermi = max(efermi, fermi)

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
            







def pdos_plot(elements,spin=True,sigma=0.003,vbm=0,grid=True,xmin=-3,xmax=6,figw=25,figh=10):
    """
    
    """

        
    k=0
    names={}
    plots={}

    if spin==True:
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
    elif spin==False:
        for element in elements:
            k=k+1
            for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*-k'+str(k)+'*'+'.pdos'):
                    alpha=file
                elif fnmatch.fnmatch(file, '*k'+str(k)+'*'+'.pdos'):
                    alpha=file

            element_no=f"element{k}"
            
            names[element_no] = smeared_pdos([alpha],sigma)


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

        if spin==True:
            orbs = (len(names[f'element{n}'][0])-1)/2
        elif spin==False:
            orbs = (len(names[f'element{n}'][0])-1)

        while i < len(names[f'element{n}']):

            x.append(names[f'element{n}'][i][0])
            
            s.append(names[f'element{n}'][i][1])

            if orbs == 1:
                y = [s]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][2]*-1)
                    y_2 = [s_2]

            if orbs == 2:
                p.append(names[f'element{n}'][i][2])
                y = [s, p]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][3]*-1)
                    p_2.append(names[f'element{n}'][i][4]*-1)
                    y_2 = [s_2, p_2]

            if orbs == 3:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                y = [s, p, d]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][4]*-1)
                    p_2.append(names[f'element{n}'][i][5]*-1)
                    d_2.append(names[f'element{n}'][i][6]*-1)
                    y_2 = [s_2, p_2, d_2]

            if orbs == 4:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                f.append(names[f'element{n}'][i][4])
                y = [s, p, d, f]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][5]*-1)
                    p_2.append(names[f'element{n}'][i][6]*-1)
                    d_2.append(names[f'element{n}'][i][7]*-1)
                    f_2.append(names[f'element{n}'][i][8]*-1)
                    y_2 = [s_2, p_2, d_2, f_2]
                        

            i=i+1

        element_name = f'{elements[n-1]}'
        plots[f'{element_name}_x'] = x
        plots[f'{element_name}_alpha_tot'] = sum(map(np.array, y))
        if spin==True:
            plots[f'{element_name}_beta_tot']  = sum(map(np.array, y_2))


        plt.figure(figsize=(figw,figh))#.set_facecolor('#BCC2C3')
        plt.title(f"{elements[n-1]} Density of States", size=30) 
        plt.xlabel('Energy (eV)', size=20) 
        plt.ylabel('Denisty of States (arb.)', size=20)


        plt.axvline(x=vbm, color='k', linestyle='--', alpha=0.5)

        plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
        plt.axhline(y=0, color='k', linestyle='-', alpha=0.5)
        

        plt.plot(x, s, marker = '',c=ELEMENT_COLOURS[elements[n-1]],linestyle='solid', label=f'{elements[n-1]}_s')
        if spin==True:
            plt.plot(x, s_2, marker = '',c=ELEMENT_COLOURS[elements[n-1]],linestyle='solid')
        
        if spin==True:
            ymax = max(max([y for x, y in zip(x, s) if xmin <= x <= xmax]), max([y for x, y in zip(x, s_2) if xmin <= x <= xmax], key=abs))
        elif spin==False:
            ymax = max([y for x, y in zip(x, s) if xmin <= x <= xmax])
        if orbs >= 2:
            plt.plot(x, p, marker = '', label=f'{elements[n-1]}_p',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=0.5),linestyle='solid')
            if spin==True:
                plt.plot(x, p_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=0.5),linestyle='solid')
                ymax = max(ymax, max([y for x, y in zip(x, p) if xmin <= x <= xmax]), max([y for x, y in zip(x, p_2) if xmin <= x <= xmax], key=abs))
            elif spin==False:
                ymax = max([y for x, y in zip(x, s) if xmin <= x <= xmax])
            if orbs >= 3:
                plt.plot(x, d, marker = '', label=f'{elements[n-1]}_d',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.2),linestyle='solid')
                if spin==True:
                    plt.plot(x, d_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.2),linestyle='solid')
                    ymax = max(ymax, max([y for x, y in zip(x, d) if xmin <= x <= xmax]), max([y for x, y in zip(x, d_2) if xmin <= x <= xmax], key=abs))
                elif spin==False:
                    ymax = max([y for x, y in zip(x, s) if xmin <= x <= xmax])
                if orbs >= 4:
                    plt.plot(x, f, marker = '', label=f'{elements[n-1]}_f',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.5),linestyle='solid')
                    if spin==True:
                        plt.plot(x, f_2, marker = '',c=lighten_hex_colour(ELEMENT_COLOURS[elements[n-1]],factor=-0.5),linestyle='solid')
                        ymax = max(ymax, max([y for x, y in zip(x, f) if xmin <= x <= xmax]), max([y for x, y in zip(x, f_2) if xmin <= x <= xmax], key=abs))
                    elif spin==False:
                        ymax = max([y for x, y in zip(x, s) if xmin <= x <= xmax])

        if grid==True:
            plt.grid(visible=1)
        else:
            plt.grid(visible=0)

        if ymax >= 0.1:
            ymax = np.ceil(ymax*10)/10
        else:
            ymax = np.ceil(ymax*100)/100


        if ymax > 1.0:
            plt.yticks(np.linspace(-ymax,ymax,int(((ymax*5) + 1)/2)*2 + 1),labels=[])
        elif 0.5 < ymax <= 1.0:
            plt.yticks(np.linspace(-ymax,ymax,int(((ymax*20) + 1)/2)*2 + 1),labels=[])
        elif 0.2 < ymax <= 0.5:
            plt.yticks(np.linspace(-ymax,ymax,int(((ymax*40) + 1)/2)*2 + 1),labels=[])
        elif 0.05 < ymax <= 0.2:
            plt.yticks(np.linspace(-ymax,ymax,int(((ymax*80) + 1)/2)*2 + 1),labels=[])
        elif 0.02 < ymax <= 0.05:
            plt.yticks(np.linspace(-ymax,ymax,int(((ymax*400) + 1)/2)*2 + 1),labels=[])
        else:
           plt.yticks(np.linspace(-ymax,ymax,int(((ymax*1000) + 1)/2)*2 + 1),labels=[])
        #plt.yticks([])
        if grid==True:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax)*4)+1))
        else:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax))+1))
        
        plt.xlim([xmin,xmax])
        if spin==True:
            plt.ylim([-ymax,ymax])
        elif spin==False:
            plt.ylim([0,ymax])

        plt.legend(markerscale=10.0, fontsize=20)

        plt.savefig(f'{elements[n-1]} Density of States', bbox_inches='tight')

        plt.show()

        n=n+1

    
    plt.figure(figsize=(figw,figh))
    plt.title("Total Density of States", size=40) 
    plt.xlabel('Energy (eV)', size=30) 
    plt.ylabel('Denisty of States (arb.)', size=30)

    if spin==True:
        ymax = max(max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_beta_tot']) if xmin <= x <= xmax], key=abs))
    elif spin==False:
        ymax = max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_alpha_tot']) if xmin <= x <= xmax])
    plt.axvline(x=vbm, color='k', linestyle='--', alpha=0.5)

    plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    plt.axhline(y=0, color='k', linestyle='-', alpha=0.5)

    j=1
    while j <= len(elements):
        if spin==True:
            ymax = max(ymax,max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot']) if xmin <= x <= xmax], key=abs))
        if spin==False:
            ymax = max(ymax,max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot']) if xmin <= x <= xmax]), key=abs)
        plt.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot'], marker = '',lw=2,c=ELEMENT_COLOURS[elements[j-1]], label=f'{elements[j-1]}')
        if spin==True:
            plt.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot'],lw=2,c=ELEMENT_COLOURS[elements[j-1]], marker = '')
        j=j+1

    if grid==True:
        plt.grid(visible=1)
    else:
        plt.grid(visible=0)

    ymax = np.ceil(ymax*10)/10
    
    if ymax > 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*5) + 1)/2)*2 + 1),labels=[])
    if 0.5 < ymax <= 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*20) + 1)/2)*2 + 1),labels=[])
    if 0.2 < ymax <= 0.5:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*40) + 1)/2)*2 + 1),labels=[])
    if ymax <= 0.2:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*80) + 1)/2)*2 + 1),labels=[])
        
    if grid==True:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax)*4)+1),fontsize=20)
    else:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax))+1),fontsize=20)

    plt.xlim([xmin,xmax])
    if spin==True:
        plt.ylim([-ymax,ymax])
    if spin==False:
        plt.ylim([0,ymax])

    plt.legend(markerscale=10.0, fontsize=20)

    plt.savefig('Total Density of States', bbox_inches='tight')

    plt.show()



def ipr_pdos_plot(elements,spin=True,sigma=0.003,vbm=0,grid=True,xmin=-3,xmax=6,figw=25,figh=10):
    """
    
    """

        
    k=0
    names={}
    plots={}

    if spin==True:
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
    elif spin==False:
        for element in elements:
            k=k+1
            for file in os.listdir('.'):
                if fnmatch.fnmatch(file, '*-k'+str(k)+'*'+'.pdos'):
                    alpha=file
                elif fnmatch.fnmatch(file, '*k'+str(k)+'*'+'.pdos'):
                    alpha=file

            element_no=f"element{k}"
            
            names[element_no] = smeared_pdos([alpha],sigma)


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

        if spin==True:
            orbs = (len(names[f'element{n}'][0])-1)/2
        elif spin==False:
            orbs = (len(names[f'element{n}'][0])-1)

        while i < len(names[f'element{n}']):

            x.append(names[f'element{n}'][i][0])
            
            s.append(names[f'element{n}'][i][1])

            if orbs == 1:
                y = [s]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][2]*-1)
                    y_2 = [s_2]

            if orbs == 2:
                p.append(names[f'element{n}'][i][2])
                y = [s, p]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][3]*-1)
                    p_2.append(names[f'element{n}'][i][4]*-1)
                    y_2 = [s_2, p_2]

            if orbs == 3:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                y = [s, p, d]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][4]*-1)
                    p_2.append(names[f'element{n}'][i][5]*-1)
                    d_2.append(names[f'element{n}'][i][6]*-1)
                    y_2 = [s_2, p_2, d_2]

            if orbs == 4:
                p.append(names[f'element{n}'][i][2])
                d.append(names[f'element{n}'][i][3])
                f.append(names[f'element{n}'][i][4])
                y = [s, p, d, f]
                if spin==True:
                    s_2.append(names[f'element{n}'][i][5]*-1)
                    p_2.append(names[f'element{n}'][i][6]*-1)
                    d_2.append(names[f'element{n}'][i][7]*-1)
                    f_2.append(names[f'element{n}'][i][8]*-1)
                    y_2 = [s_2, p_2, d_2, f_2]
                        

            i=i+1

        element_name = f'{elements[n-1]}'
        plots[f'{element_name}_x'] = x
        plots[f'{element_name}_alpha_tot'] = sum(map(np.array, y))
        if spin==True:
            plots[f'{element_name}_beta_tot']  = sum(map(np.array, y_2))

    
    plt.figure(figsize=(figw,figh))
    plt.title("Total Density of States", size=40)
    
    ax = plt.subplot()
    ax.set_xlabel('Energy (eV)', size=30) 
    ax.set_ylabel('Denisty of States (arb.)', size=30)

    if spin==True:
        ymax = max(max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_beta_tot']) if xmin <= x <= xmax], key=abs))
    elif spin==False:
        ymax = max([y for x, y in zip(plots[f'{elements[0]}_x'], plots[f'{elements[0]}_alpha_tot']) if xmin <= x <= xmax])
    #plt.axvline(x=vbm, color='k', linestyle='--', alpha=0.5)

    #plt.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax.axhline(y=0, color='k', linestyle='-', alpha=0.5)

    j=1
    while j <= len(elements):
        if spin==True:
            ymax = max(ymax,max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot']) if xmin <= x <= xmax]), max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot']) if xmin <= x <= xmax], key=abs))
        if spin==False:
            ymax = max(ymax,max([y for x, y in zip(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot']) if xmin <= x <= xmax]), key=abs)
        ax.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_alpha_tot'], marker = '',lw=2,c=ELEMENT_COLOURS[elements[j-1]], label=f'{elements[j-1]}')
        if spin==True:
            ax.plot(plots[f'{elements[j-1]}_x'], plots[f'{elements[j-1]}_beta_tot'],lw=2,c=ELEMENT_COLOURS[elements[j-1]], marker = '')
        j=j+1

    if grid==True:
        plt.grid(visible=1)
    else:
        plt.grid(visible=0)

    ymax = np.ceil(ymax*10)/10
    
    if ymax > 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*5) + 1)/2)*2 + 1),labels=[])
    if 0.5 < ymax <= 1.0:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*20) + 1)/2)*2 + 1),labels=[])
    if 0.2 < ymax <= 0.5:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*40) + 1)/2)*2 + 1),labels=[])
    if ymax <= 0.2:
        plt.yticks(np.linspace(-ymax,ymax,int(((ymax*80) + 1)/2)*2 + 1),labels=[])
        
    if grid==True:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax)*4)+1),fontsize=20)
    else:
            plt.xticks(np.linspace(xmin,xmax,int((abs(xmin)+xmax))+1),fontsize=20)

    plt.xlim([xmin,xmax])
    if spin==True:
        plt.ylim([-ymax,ymax])
    if spin==False:
        plt.ylim([0,ymax])

    ax2 = ax.twinx()
    ax2.set_ylabel('IPR', size=30)

    ipr_alpha_x = []
    ipr_alpha_y = []
    ipr_beta_x = []
    ipr_beta_y = []

    with open('./IPR_output-Alpha.dat', 'r') as file:
        for line in file:
            columns = line.split()
            ipr_x_alpha.append(float(columns[0])*Ha_to_eV-vbm)
            ipr_y_alpha.append(float(columns[1]))
    with open('./IPR_output-Beta.dat', 'r') as file:
        for line in file:
            columns = line.split()
            ipr_x_beta.append(float(columns[0])*Ha_to_eV-vbm)
            ipr_y_beta.append(float(columns[1]))
   
    for xi, yi in zip(ipr_x_alpha, ipr_y_alpha):
        ax2.vlines(xi, 0, yi, colors='blue', linewidth=1)
    for xi, yi in zip(ipr_x_beta, ipr_y_beta):
        ax2.vlines(xi, 0, yi, colors='orange', linewidth=1)

    
    plt.legend(markerscale=10.0, fontsize=20)

    plt.savefig('Total Density of States', bbox_inches='tight')

    plt.show()
