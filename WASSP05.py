# Wannier & vASP Postprocessing module
# Utilities for VASP i coded or curated during my PhD
# This code is heavily based in the one by Martín Gutierrez for QEspresso
# PyBand by Qijing Zheng & Chengcheng Xiao and uses pymatgen.
# version =  1.0    date = 1/04/2022
# Author  = Irián Sánchez-Ramírez mail = irian.sanchez@dipc.org


from pymatgen.io.vasp.outputs import Vasprun
import pymatgen.core.structure as pst
import pymatgen.symmetry.analyzer as psa
from pymatgen.electronic_structure.plotter import BSPlotter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
from scipy.interpolate import make_interp_spline
import itertools as it
import re

rc('text', usetex=True)
rc('font', size=14)
rc('legend', fontsize=13)
rc('text.latex', preamble=r'\usepackage{cmbright}')

### Before Wannierization 

def tensorify(file,spin = True):
    
    """ Parse vasprun and obtain everything needed for plotting,
        this function generates a tensor which is structured this way:
        tensor[atom,orbital,spin] = DOS(atom, orbital, spin).
        Aditionally it throws also a list of atoms, orbitals, spins, energies
        and Fermi energy """

    # Use pymatgen to handle vasprun.xml

    output = Vasprun(file,parse_dos = True)
    pdos = output.complete_dos.pdos
    efermi = output.complete_dos.efermi
    energies = output.complete_dos.energies

    # Initialize the tensor

    if spin == True:
        tensor = np.zeros((len(energies),
                            len(pdos.keys()),
                            len(pdos[list(pdos.keys())[0]].keys()),
                            2),
                            dtype = np.float64)
    else:
        tensor = np.zeros((len(energies),
                            len(pdos.keys()),
                            len(pdos[list(pdos.keys())[0]].keys()),
                            1),
                            dtype = np.float64)
    
    
    # Fill it
    
    atoms = [key for key in pdos.keys()]
    orbitals = [key for key in pdos[atoms[0]]]
    
    for (ii,atom) in enumerate(atoms):
        for (i,orbital) in enumerate(pdos[atom].keys()):
            for (j,sp) in enumerate(pdos[atom][orbital].keys()):
                tensor[:,ii,i,j] = pdos[atom][orbital][sp]

    atoms = [str(atom)[-2::].strip() for atom in atoms]
    orbitals = [str(orbital) for orbital in orbitals]
    
    if spin == True:
        spins = ['up','down']

    return tensor, atoms, orbitals, spins, energies, efermi

def plot_pdos(file,_atoms = None,_orbitals = None,fig_size = (18,12),e_window = (-10,10)):
    
    """ Main function for plotting the partial density of states, for desired atoms and orbitals
    at all the different wyckoff positions.
    The usage is straightforward:
    _atoms = the atoms whose DOS you wanna plot as a list of strings e.g. ["I","S","Ra"]
    _orbitals = the orbitals of these atoms e.g. ["s","px","d","d_xy"]
    e_window = range of energies to plot.
    Everything is smoothed out :3
    """
    # Colors & linestiles
    ls_borb = ["--","-.",":",'densely dotted']
    ls_sorb = []
    # Get info
    tensor, atoms, orbitals, spins, energies, efermi = tensorify(file)
    structure = pst.IStructure.from_file(file)
    sga = psa.SpacegroupAnalyzer(structure)
    wyckoffs = sga.get_symmetry_dataset()['wyckoffs']
    
    smooth_energy = np.linspace(np.min(energies),np.max(energies),5000)
    smooth_array = smooth_energy    
    
    groups = {}
    unique_atoms = list(set(atoms))
    unique_wyckoffs = list(set(wyckoffs))

    # Generate an empty dictionary with combinations of atoms and wyckoffs
    for atom in unique_atoms:
        for wyckoff in unique_wyckoffs:
            groups[(atom,wyckoff)] = []
    # Fill only the ones matching --> later using a if != [] condition
    for (ii,atom) in enumerate(atoms):
        groups[(atom,wyckoffs[ii])].append(ii)
    # Parse wanted atoms & orbitals 
    if _orbitals != None:
        bigorbs = list(set([orbital[0] for orbital in _orbitals]))
    else:
        bigorbs = list(set([orbital[0] for orbital in orbitals]))
    comp_orbitals = []
    comp_atoms = []
    comp_bigorbs = []

    # Compare given input and existing atoms and orbitals.

    if _orbitals != None:
        for orbital in _orbitals:
            if orbital in orbitals:
                comp_orbitals.append(int(np.where(np.asarray(orbitals) == orbital)[0]))
            if orbital in bigorbs:
                comp_bigorbs.append(orbital)
    else: comp_bigorbs = bigorbs

    if _atoms != None:
        for (jj,atom) in enumerate(_atoms):
            if atom in atoms:
                comp_atoms.append(list(np.where(np.asarray(atoms) == atom)[0]))
    else: 
        _atoms = atoms
        for (jj,atom) in enumerate(_atoms):
            if atom in atoms:
                comp_atoms = list(np.where(np.asarray(atoms) == atom)[0])               
    
    comp_atoms = [j for i in comp_atoms for j in i]

    
    # Plot the small orbitals (aka: s,px,py,pz,dxy ...)

    fig = plt.figure(figsize = fig_size)
    if _orbitals != None:
        for atom_pos in comp_atoms:
            for wyckoff in unique_wyckoffs:
                if groups[(atoms[atom_pos],wyckoff)] != []:
                    for orbital_pos in comp_orbitals:
                        spl_tmp = make_interp_spline(energies-efermi,sum([tensor[:,val,orbital_pos,0] for val in groups[(atoms[atom_pos],wyckoff)]]) +
                                                                        sum([tensor[:,val,orbital_pos,1] for val in groups[(atoms[atom_pos],wyckoff)]]), k = 3)
                        plt.plot(smooth_energy - efermi, spl_tmp(smooth_energy-efermi),label =atoms[atom_pos]+" "+wyckoff+" "+orbitals[orbital_pos])
                    groups[(atoms[atom_pos],wyckoff)] = []
        groups = {}
        unique_atoms = list(set(atoms))
        unique_wyckoffs = list(set(wyckoffs))

    # Plot the big orbitals (aka: p = px+ py +pz ...)

        # Generate an empty dictionary with combinations of atoms and wyckoffs
        for atom in unique_atoms:
            for wyckoff in unique_wyckoffs:
                groups[(atom,wyckoff)] = []
        # Fill only the ones matching --> later using a if != [] condition
        for (ii,atom) in enumerate(atoms):
            groups[(atom,wyckoffs[ii])].append(ii) 
    
    for atom in unique_atoms:
        if atom in _atoms:
            for wyckoff in unique_wyckoffs:
                if groups[(atom,wyckoff)] != []:
                    
                    if "s" in comp_bigorbs:
                        temp_s = np.zeros(len(energies),dtype = np.float64)
                    if "p" in comp_bigorbs:
                        temp_p = np.zeros(len(energies),dtype = np.float64)
                    if "d" in comp_bigorbs:
                        temp_d = np.zeros(len(energies),dtype = np.float64)
                    if "f" in comp_bigorbs:
                        temp_f = np.zeros(len(energies),dtype = np.float64)

                    for bigorb in comp_bigorbs:
                        for (jj,orbital) in enumerate(orbitals):
                            if orbital[0] == "s" and orbital[0] == bigorb:
                                temp_s +=  sum([tensor[:,val,jj,0] for val in groups[(atom,wyckoff)]])
                            elif orbital[0] == "p" and orbital[0] == bigorb:
                                temp_p +=  sum([tensor[:,val,jj,0] for val in groups[(atom,wyckoff)]])
                            elif orbital[0] == "d" and orbital[0] == bigorb:
                                temp_d +=  sum([tensor[:,val,jj,0] for val in groups[(atom,wyckoff)]])
                            elif orbital[0] == "f" and orbital[0] == bigorb:
                                temp_f +=  sum([tensor[:,val,jj,0] for val in groups[(atom,wyckoff)]])
                    
                    groups[(atom,wyckoff)] = []

                    if "s" in comp_bigorbs:
                        if _orbitals != None and "s" not in _orbitals:
                            spl_s = make_interp_spline(energies-efermi,temp_s,k=3)
                            plt.plot(smooth_array-efermi,spl_s(smooth_array-efermi),label = atom+" "+wyckoff+" "+"s",ls = "--")
                        else:
                            spl_s = make_interp_spline(energies-efermi,temp_s,k=3)
                            plt.plot(smooth_array-efermi,spl_s(smooth_array-efermi),label = atom+" "+wyckoff+" "+"s", ls = "--")
                    if "p" in comp_bigorbs:
                        spl_p = make_interp_spline(energies-efermi,temp_p,k=3)
                        plt.plot(smooth_array-efermi,spl_p(smooth_array-efermi),label = atom+" "+wyckoff+" "+"p",ls = "-.")
                    if "d" in comp_bigorbs:
                        spl_d = make_interp_spline(energies-efermi,temp_d,k=3)
                        plt.plot(smooth_array-efermi,spl_d(smooth_array-efermi),label = atom+" "+wyckoff+" "+"d", ls = ":")
                    if "f" in comp_bigorbs:
                        spl_f = make_interp_spline(energies-efermi,temp_f,k=3)         
                        plt.plot(smooth_array-efermi,spl_f(smooth_array-efermi),label = atom+" "+wyckoff+" "+"f", ls = ".")

    # Plot the total DOS

    total_dos = np.zeros(len(energies),np.float64)
        
    for (ii,atom) in enumerate(atoms):
        for (jj,orbital) in enumerate(orbitals):
            for (kk, spin) in enumerate(spins):
                total_dos += tensor[:,ii,jj,kk]


    spl_tot = make_interp_spline(energies-efermi,total_dos,k=3)

    plt.plot(smooth_array-efermi,spl_tot(smooth_array-efermi),label = "Total")

    idx_min = (np.abs(energies - e_window[0])).argmin()
    idx_max = (np.abs(energies - e_window[1])).argmin()
    max_dos = np.max(total_dos[idx_min:idx_max])
    plt.xlim(e_window[0]-0.01*(e_window[1]-e_window[0]),e_window[1]+0.01*(e_window[1]-e_window[0]))

    plt.axvline(x = 0.0,linestyle = '--', c = "grey")
    plt.legend(prop={'size': 25})
    plt.xlim(e_window[0],e_window[1])
    plt.ylim(0,max_dos + max_dos*0.05)
    plt.xlabel("Enegy [$eV$]",fontsize = 25)
    plt.ylabel("DOS a.u.",fontsize = 25)
    plt.show()


def band_counter_gamma(file = "OUTCAR", emin = 0.0, emax =0.0):

    """ Counts the number of bands in a certain energy window at Gamma point.
        Its fast but only gives a first guess for the energy windows.
        Energy window must be suplied in nergies regulated by fermi energy:
            --> Regardless Ef = 5.0 eV, think of it being Ef = 0.0 eV
        Recomended upplied window values for .win will be indicated as output."""

    outcar = [line for line in open(file)]
    band_lines = np.zeros(2,dtype = int)
    counter = 0
    my_formatter = "{0:4.2f}"

    for (ii,line) in enumerate(outcar):
        if "E-fermi :" in line:
            efermi = np.float64(line.split()[2])
        if "k-point     1" in line:
            band_lines[0] = ii
        if "k-point     2" in line:
            band_lines[1] = ii
            break
    
    for (ii,line) in enumerate(outcar[band_lines[0]+2:band_lines[1]-1]):
        if np.float64(line.split()[1]) >= emin+efermi and np.float64(line.split()[1]) <= emax+efermi:
            counter += 1
    print("Efermi="+str(efermi))
    print("The number of bands between "+my_formatter.format(emin)+" eV ("+my_formatter.format(emin+efermi)+" eV) and "+my_formatter.format(emax)+" eV ("+my_formatter.format(emax+efermi)+" eV) is "+str(counter)+" at Gamma point.")
    

def band_counter(file = "vasprun.xml",emin =0.0, emax= 0.0):

    """ Counts the minimun number of bands for a energy window in the 1bz.
        Its slower but gives a good guess for the energy windows.
        Energy window must be suplied in nergies regulated by fermi energy:
            --> Regardless Ef = 5.0 eV, think of it being Ef = 0.0 eV
        Recomended upplied window values for .win will be indicated as output. """

    output = Vasprun(file,parse_eigen = True)    
    eigdic = output.eigenvalues
    nkpts = len(eigdic[list(eigdic.keys())[0]])
    bands = len(eigdic[list(eigdic.keys())[0]][0])
    efermi = output.efermi 
    bandata = np.zeros(nkpts,dtype = int)
    iterator = it.product(np.linspace(0,nkpts-1,nkpts,dtype = int),np.linspace(0,bands-1,bands,dtype = int))
    my_formatter = "{0:4.2f}"

    for comb in iterator:
        if emin+efermi <= eigdic[list(eigdic.keys())[0]][comb[0]][comb[1]][0] <= emax+efermi:
            bandata[comb[0]] += 1

    num_bands = np.min(bandata)
   
    print("Efermi = "+str(efermi)+".")
    print("Total bands = "+str(bands)+".")
    print("The number of bands between "+my_formatter.format(emin)+" eV ("+my_formatter.format(emin+efermi)+" eV) and "+my_formatter.format(emax)+" eV ("+my_formatter.format(emax+efermi)+" eV) is "+str(num_bands)+".")

### After Wannierization

def plot_wannierbands(file_dat = "wannier90_band.dat", gnu = "wannier90_band.gnu",efermi = 0.0, e_window = None, fig_size = (15,8),savename = "wannierbands.png"):
    
    """ Plotting the wannier bands using the .dat and .gnu outputs from a wannier90 run 
        Mainly this code just interprets how the .dat file is written and translates it to 
        numpy and matplotlib ploteable data """

    data = np.loadtxt(file_dat)
    block_length = np.where(data[:,0]==0.0)[0][1]                        # Number of lines in each block of .dat
    block_number = np.shape(np.where(data[:,0]==0.0))[1]                 # Number of blocks 
    data = np.reshape(data,(block_length,2*block_number),order = "F")    # Reshape array
    num_wann = block_number-1
    data = np.delete(data,np.s_[0:num_wann],1)

    fig, ax = plt.subplots(figsize = fig_size)
    ax.plot(data[:,0],data[:,1],color = "red")
    for i in np.linspace(3,num_wann,num_wann-2):
        if np.max(data[:,int(i-1):int(i)]) >= efermi+0.5:
            ax.plot(data[:,0],data[:,int(i-1):int(i)]-efermi,color = "red")
        else:
            ax.plot(data[:,0],data[:,int(i-1):int(i)]-efermi,color = "blue")

    ax.set_xticks([])
    
    kpts = [line for line in open(gnu) if line.strip()]
    x_ticks =  []
    x_labels = []
    
    for line in kpts:
        if "set xtics" in line:
            label_info = line
    
    for word in [words for words in label_info[11:-2].split(",")]:
        x_ticks.append(np.float64(word.split()[1]))
        if word.split()[0][1:-1] == 'G':
            x_labels.append("$\Gamma$")
        else:
            x_labels.append(word.split()[0][1:-1])

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels)

    for line in x_ticks:
        ax.axvline(x = line, ls  ="--", color = "grey", lw = 0.5)
    
    ax.axhline(y = 0.0, ls = "-",color ="grey", lw=  0.6)
    ax.set_xlim(x_ticks[0],x_ticks[-1]) 
    if e_window != None:
        ax.set_ylim(e_window[0],e_window[1])
    
    if efermi != 0.0:
        ax.set_ylabel("$E-E_f\;(eV)$")
    else: 
        ax.set_ylabel("$E;(eV)$")

    plt.savefig(savename)
    plt.show()

# Functions for plotting VASP output:

def kpath_name_parse(KPATH_STR):

    '''
    Parse the kpath string
    '''

    KPATH_STR = KPATH_STR.upper()
    # legal kpath separators: blank space, comma, hypen or semicolon
    KPATH_SEPARATORS = ' ,-;'
    # Greek Letters Dictionaries
    GREEK_KPTS = {
        'G':      r'$\mathrm{\mathsf{\Gamma}}$',
        'GAMMA':  r'$\mathrm{\mathsf{\Gamma}}$',
        'Gamma':  r'$\mathrm{\mathsf{\Gamma}}$',
        'DELTA':  r'$\mathrm{\mathsf{\Delta}}$',
        'LAMBDA': r'$\mathrm{\mathsf{\Lambda}}$',
        'SIGMA':  r'$\mathrm{\mathsf{\Sigma}}$',
    }

    # If any of the kpath separators is in the kpath string
    if any([s in KPATH_STR for s in KPATH_SEPARATORS]):
        kname = [
            GREEK_KPTS[x] if x in GREEK_KPTS else 
        r'$\mathrm{{\mathsf{{{}}}}}$'.format(x)
            for x in re.sub('['+KPATH_SEPARATORS+']', ' ', KPATH_STR).split()
        ]
    else:
        kname = [
            GREEK_KPTS[x] if x in GREEK_KPTS else 
        r'$\mathrm{{\mathsf{{{}}}}}$'.format(x)
            for x in KPATH_STR
        ]

    return kname

def get_bandInfo1(inFile='OUTCAR',kpointfile = "KPOINTS"):

    """
    extract band energies from OUTCAR
    """

    outcar = [line for line in open(inFile) if line.strip()]

    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

        if 'NELECT' in line:
            nelect = float(line.split()[2])
            # break

    # basis vector of reciprocal lattice
    # B = np.array([line.split()[3:] for line in outcar[ibasis:ibasis+3]],

    # When the supercell is too large, spaces are missing between real space
    # lattice constants. A bug found out by Wei Xie (weixie4@gmail.com).
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]],
                 dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:, :3]
    wkpts = tmp[:, -1]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi+1:LineEfermi + N+1]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))


    kp = open(kpointfile).readlines()

    if kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0])
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)

        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end, :] - vkpts[start, :]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds = np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1], ]))

    return kpt_path, bands, efermi, kpt_bounds, wkpts, nelect

def get_bandInfo2(inFile='OUTCAR',kpointfile = "KPOINTS"):

    """
    extract band energies from OUTCAR
    """

    outcar = [line for line in open(inFile) if line.strip()]

    for ii, line in enumerate(outcar):
        if 'NKPTS =' in line:
            nkpts = int(line.split()[3])
            nband = int(line.split()[-1])

        if 'ISPIN  =' in line:
            ispin = int(line.split()[2])

        if "k-points in reciprocal lattice and weights" in line:
            Lvkpts = ii + 1

        if 'reciprocal lattice vectors' in line:
            ibasis = ii + 1

        if 'E-fermi' in line:
            efermi = float(line.split()[2])
            LineEfermi = ii + 1
            # break

        if 'NELECT' in line:
            nelect = float(line.split()[2])
            # break

    # basis vector of reciprocal lattice
    # B = np.array([line.split()[3:] for line in outcar[ibasis:ibasis+3]],

    # When the supercell is too large, spaces are missing between real space
    # lattice constants. A bug found out by Wei Xie (weixie4@gmail.com).
    B = np.array([line.split()[-3:] for line in outcar[ibasis:ibasis+3]],
                 dtype=float)
    # k-points vectors and weights
    tmp = np.array([line.split() for line in outcar[Lvkpts:Lvkpts+nkpts]],
                   dtype=float)
    vkpts = tmp[:, :3]
    wkpts = tmp[:, -1]

    # for ispin = 2, there are two extra lines "spin component..."
    N = (nband + 2) * nkpts * ispin + (ispin - 1) * 2
    bands = []
    # vkpts = []
    for line in outcar[LineEfermi:LineEfermi + N]:
        if 'spin component' in line or 'band No.' in line:
            continue
        if 'k-point' in line:
            # vkpts += [line.split()[3:]]
            continue
        bands.append(float(line.split()[1]))

    bands = np.array(bands, dtype=float).reshape((ispin, nkpts, nband))


    kp = open(kpointfile).readlines()

    if kp[2][0].upper() == 'L':
        Nk_in_seg = int(kp[1].split()[0])
        Nseg = nkpts // Nk_in_seg
        vkpt_diff = np.zeros_like(vkpts, dtype=float)

        for ii in range(Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            vkpt_diff[start:end, :] = vkpts[start:end, :] - vkpts[start, :]

        kpt_path = np.linalg.norm(np.dot(vkpt_diff, B), axis=1)
        # kpt_path = np.sqrt(np.sum(np.dot(vkpt_diff, B)**2, axis=1))
        for ii in range(1, Nseg):
            start = ii * Nk_in_seg
            end = (ii + 1) * Nk_in_seg
            kpt_path[start:end] += kpt_path[start-1]

        # kpt_path /= kpt_path[-1]
        kpt_bounds = np.concatenate((kpt_path[0::Nk_in_seg], [kpt_path[-1], ]))

    return kpt_path, bands, efermi, kpt_bounds, wkpts, nelect


def bandplot(kpath, bands, efermi, kpt_bounds, nelect, kpointfile =  "KPOINTS",fig_size = (15,8),e_window = None,line_w = 0.5,title = "bands.png"):
    '''
    Use matplotlib to plot band structure
    '''

    width, height = fig_size
    if e_window != None:
        ymin, ymax = e_window[0], e_window[1]
    else:
        ymin, ymax = -4.0,4.0

    fig = plt.figure()
    fig.set_size_inches(width, height)
    ax = plt.subplot(111)

    nspin, nkpts, nbands = bands.shape

    clrs = ['r', 'b']

    for Ispin in range(nspin):
        for Iband in range(nbands):
            
            # if Iband == 0 else line.get_color()
            lc = None if Iband == 0 else line.get_color()
            
            #if nspin == 1:
            #    new_nelect = nelect/2
            if float(Iband) >= nelect:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], color='red',lw=line_w, zorder=0,
                                alpha=0.8)
            else:
                line, = ax.plot(kpath, bands[Ispin, :, Iband], lw=line_w, zorder=0,
                                alpha=0.8,
                                color='blue')

    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)

    # add extra horizontal/vertical lines

    ax.set_ylabel('$E - E_f$ [eV]',  # fontsize='small',
                  labelpad=5)
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(kpath.min(), kpath.max())

    ax.set_xticks(kpt_bounds)

    # Read and use kpoint file

    with open(kpointfile,'r') as KPointsFile:
        TmpFlag = 0;
        TmpLabels = [];
        for TmpLine in KPointsFile:
            TmpLine = TmpLine.strip()
            if TmpFlag == 1:                    
                TmpLine = re.sub(r'^.*\!\s?', '', TmpLine)
                TmpLine.strip()
                if TmpLine != "":
                    if TmpLine == "G":
                        TmpLabels.append(r'$\mathrm{{\mathsf{\Gamma}}}$')
                    else:
                        TmpLabels.append(r'$\mathrm{\mathsf{'+TmpLine+'}}$')
            if (TmpLine == "reciprocal") | (TmpLine == "rec"):
                TmpFlag = 1
        TmpLabels2 = [TmpLabels[0]]
        TmpIndex = 1
        while TmpIndex < (len(TmpLabels) - 1):
            if TmpLabels[TmpIndex + 1] == TmpLabels[TmpIndex]:
                TmpLabels2.append(TmpLabels[TmpIndex])
            else:
                TmpLabels2.append(TmpLabels[TmpIndex]+'|'+TmpLabels[TmpIndex + 1])
            TmpIndex += 2
        TmpLabels2.append(TmpLabels[len(TmpLabels) - 1])
        ax.set_xticklabels(TmpLabels2)

    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.axhline(y=0, xmin=0, xmax=1, linestyle='dotted', color='black', linewidth=0.5)
    plt.tight_layout(pad=0.20)
    plt.show()
    plt.savefig(title)
    return kpath, bands

def plot_vaspbands(outcar = "OUTCAR", kpoints = "KPOINTS"):

    """ Plot vasp bands using the aforedefined functions """

    # Use non-interactive backend in case there is no display
    mpl.use('agg')
    mpl.rcParams['axes.unicode_minus'] = False
    try:    
        kpath, bands, efermi, kpt_bounds, wkpts, nelect = get_bandInfo1(outcar,kpoints)
    except IndexError:
        kpath, bands, efermi, kpt_bounds, wkpts, nelect = get_bandInfo2(outcar,kpoints)

    bandplot(kpath,bands-efermi,efermi,kpt_bounds,nelect,kpointfile=kpoints)
    
def plot_comparison(outcar = "OUTCAR", kpoints = "KPOINTS",file_dat = "wannier90_band.dat", gnu = "wannier90_band.gnu",efermi = 0.0, fig_size = (12,8), e_window = (-4,4),savename = "comparison.png"):
     
    """ Plots the comparison by plotting the values into the same subplot but with different axis
        in order to avoid problems with array sizes """

    # Get the info for VASP
    try:
        kpath, bands, efermi, _, _, _ = get_bandInfo1(outcar,kpoints)
    except IndexError:
        kpath, bands, efermi, _, _, _ = get_bandInfo2(outcar,kpoints)
        
    # We plot both bandstructures in the same axis and then using k-labels from gnu wannier file
    fig=plt.figure(figsize=fig_size)
    ax2=fig.add_subplot(111, label="Wannier")
    ax1=fig.add_subplot(111, label="VASP", frame_on=False)

    # VASP
    for i in range(np.shape(bands)[2]):
        ax1.plot(kpath,bands[0,:,i]-efermi, color = "grey",lw = 0.65)
    
    # Wannier90
    
    data = np.loadtxt(file_dat)
    block_length = np.where(data[:,0]==0.0)[0][1]                        # Number of lines in each block of .dat
    block_number = np.shape(np.where(data[:,0]==0.0))[1]                 # Number of blocks 
    data = np.reshape(data,(block_length,2*block_number),order = "F")    # Reshape array
    num_wann = block_number-1
    data = np.delete(data,np.s_[0:num_wann],1)

    for i in np.linspace(3,num_wann,num_wann-2):
        ax2.plot(data[:,0],data[:,int(i-1):int(i)]-efermi,color = "red",ls="-.",lw = 0.8)


    ax2.set_xticks([])
    ax1.set_xticks([])
    
    kpts = [line for line in open(gnu) if line.strip()]
    x_ticks =  []
    x_labels = []
    
    for line in kpts:
        if "set xtics" in line:
            label_info = line
    
    for word in [words for words in label_info[11:-2].split(",")]:
        x_ticks.append(np.float64(word.split()[1]))
        if word.split()[0][1:-1] == 'G':
            x_labels.append("$\Gamma$")
        else:
            x_labels.append(word.split()[0][1:-1])

    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_labels)

    for line in x_ticks:
        ax2.axvline(x = line, ls  ="--", color = "grey", lw = 0.5)
    
    ax2.axhline(y = 0.0, ls = "-",color ="grey", lw=  0.6)
    ax2.set_xlim(x_ticks[0],x_ticks[-1]) 
    ax1.set_xlim(kpath[0],kpath[-1])

    if efermi != 0.0:
        ax1.set_ylabel("$E-E_f\;(eV)$")
    else: 
        ax1.set_ylabel("$E;(eV)$")

    ax2.set_ylim(e_window[0]+0.01,e_window[1]-0.01)
    ax1.set_ylim(e_window[0]+0.01,e_window[1]-0.01)
    
    vasp = mpl.lines.Line2D([0], [0],color='grey', label='VASP',ls = "-",lw = 0.65)
    wanni = mpl.lines.Line2D([0], [0],color='red', label='Wannier90',ls = "-.",lw = 0.8)
    plt.legend(handles=[vasp,wanni],loc = 1,fontsize = 15)

    plt.savefig(savename,dpi = 500)
    plt.show()

## Utilities for plotting MBJ functional vasp runs.

def scrap_data(vasprun_file,kpoint_file):

    """ Scrap data from vasprun file for comparison plot """

    vr_dic = Vasprun(vasprun_file)
    data_dict = BSPlotter(vr_dic.get_band_structure(kpoints_filename = kpoint_file)).bs_plot_data()
    energies = data_dict["energy"]["1"]
    distances = data_dict["distances"]
    tick_place = data_dict["ticks"]["distance"]
    tick_label = data_dict["ticks"]["label"]

    return energies, distances, tick_place, tick_label


def compare_MBJ(vasprun_pbe = "vasprun1.xml",
            vasprun_mbj = "vasprun2.xml",
            kpoint_file = "KPOINTS",
            e_window = (-4,4),
            fig_title = None,
            fig_name = "comparison.png"):

    """ Plot a comparison between MBJ and PBE potentials"""
    
    energies1, distances1, tick_place1, tick_label1 = scrap_data(vasprun_pbe,kpoint_file)
    energies2, distances2, _, _ = scrap_data(vasprun_mbj,kpoint_file)

    for ii, tick in enumerate(tick_label1): 
        if tick == "GAMMA" or tick =="Gamma" or tick=="G":
            tick_label1[ii] = "$\Gamma$"

    energies2 = energies2[0:-1]
    distances2 = distances2[0:-1]

    fig, ax = plt.subplots(figsize = (12,8))
    ax.set_ylim(e_window[0],e_window[1])
    ax.set_xlim(distances2[0][0],distances2[-1][-1])

    for dist, ene in zip(distances1,energies1):
        ax.plot(dist,ene.T,c = "k",lw = 0.7)

    for dist, ene in zip(distances2,energies2):
        ax.plot(dist,ene.T,c = "r",lw = 0.7, ls = ":")

    ax.set_xticks([])
    ax.set_xticks(tick_place1)
    ax.set_xticklabels(tick_label1)
    ax.axhline(y = 0, ls = "--",c="grey",lw = 0.5)

    for pos in tick_place1[0:-1]:
        ax.axvline(pos,ls = "-.",c="grey",lw = 0.4)

    ax.set_ylabel("$E-E_f(eV)$")

    pbe = mpl.lines.Line2D([0], [0],color='k', label='PBE',ls = "-",lw = 0.7)
    mbj = mpl.lines.Line2D([0], [0],color='red', label='MBJ',ls = ":",lw = 0.7)
    plt.legend(handles=[pbe,mbj],loc = 1,fontsize = 15)
    plt.title(fig_title)
    plt.savefig(fig_name,dpi = 500)

    plt.show()
    
    def wann_kpoints(file = "KPOINTS"):

    """ Generates a kpoint path from a KPOINT file of VASP nsc calculation suitable for seedname.win"""

    kpoints = [line for line in open(file) if line.strip() and line != ""][4::]
    new_kpoints = []
    temp_line = []
    for ii, line in enumerate(kpoints):
        if ii%2 == 0 and ii < len(kpoints)-1:
            if line.split()[-1] == "GAMMA":
                temp_line.append("G"+" "+" ".join(line.split()[0:3]))
            else: temp_line.append(line.split()[-1]+" "+" ".join(line.split()[0:3]))
            if kpoints[ii+1].split()[-1] == "GAMMA":
                temp_line.append("G"+" "+" ".join(kpoints[ii+1].split()[0:3]))
            else: temp_line.append(kpoints[ii+1].split()[-1]+" "+" ".join(kpoints[ii+1].split()[0:3]))
            temp_line.append("\n")
            new_kpoints.append(" ".join(temp_line))
            temp_line = []

    with open("WKPTS.txt","w") as f:
        for line in new_kpoints:
            f.write(line)

def bandplot2(kpt_path, bands, efermi, kpt_bounds, nelect, W,L,ymin, ymax,dpi,linewidth,kp_i=None,kp_f=None,title=None,fname=None):                            # y limits

    fig = plt.figure(figsize=(W,L))
    ax = plt.subplot(111)
    nspin, nkpts, nbands = bands.shape                # get spin, number of kpoints and number of bands.

    if kp_i is not None and kp_f is not None:
        km_i = np.where(kpt_path == kpt_bounds[kp_i])[0][0]
        km_f = np.where(kpt_path == kpt_bounds[kp_f])[0][0]+1
        if km_f >= len(kpt_path): km_f = km_f-1

        ax.set_xlim(kpt_path[km_i],kpt_path[km_f])
    else:
        km_i = 0
        km_f = len(kpt_path)-1
        ax.set_xlim(kpt_path[km_i],kpt_path[km_f])

    for Ispin in range(nspin):
        for Iband in range(nbands):
            lc = None if Iband == 0 else line.get_color()
            if Iband >= nelect:
                line, = ax.plot(kpt_path[km_i:km_f], bands[Ispin,:, Iband][km_i:km_f]-efermi, lw= linewidth, zorder=0,
                            alpha=0.8,
                            color='red',
                            )
            else:
                line, = ax.plot(kpt_path[km_i:km_f], bands[Ispin, :, Iband][km_i:km_f]-efermi, lw= linewidth, zorder=0,
                            alpha=0.8,
                            color='blue',
                            )

    for bd in kpt_bounds:
        ax.axvline(x=bd, ls='-', color='k', lw=0.5, alpha=0.5)

    ax.set_ylabel('$E - E_f$ [eV]',fontsize='x-large',labelpad=5)
    ax.set_ylim(ymin, ymax)



    with open('KPOINTS','r') as KPointsFile:
        TmpFlag = 0;
        TmpLabels = [];
        for TmpLine in KPointsFile:
            TmpLine = TmpLine.strip()
            if TmpFlag == 1:
                TmpLine = re.sub(r'^.*\!\s?', '', TmpLine)
                TmpLine.strip()
                if TmpLine != "":
                    if TmpLine == "G" or TmpLine == "Gamma" or TmpLine == "GAMMA":
                        TmpLabels.append(r'$\mathrm{{\mathsf{\Gamma}}}$')
                    else:
                        TmpLabels.append(r'$\mathrm{\mathsf{'+TmpLine+'}}$')
            if (TmpLine == "reciprocal") | (TmpLine == "rec"):
                TmpFlag = 1
        TmpLabels2 = [TmpLabels[0]]
        TmpIndex = 1
        while TmpIndex < (len(TmpLabels) - 1):
            if TmpLabels[TmpIndex + 1] == TmpLabels[TmpIndex]:
                TmpLabels2.append(TmpLabels[TmpIndex])
            else:
                TmpLabels2.append(TmpLabels[TmpIndex]+'|'+TmpLabels[TmpIndex + 1])
            TmpIndex += 2
        TmpLabels2.append(TmpLabels[len(TmpLabels) - 1])

        if kp_i is not None and kp_f is not None:
            ax.set_xlim(kpt_bounds[kp_i],kpt_bounds[kp_f])
            if kp_f >= len(TmpLabels2):
                ax.set_xticks(kpt_bounds[kp_i:kp_f])
                ax.set_xticklabels(TmpLabels2[kp_i:kp_f],Fontsize= 12)
            else:
                ax.set_xticks(kpt_bounds[kp_i:kp_f+1])
                ax.set_xticklabels(TmpLabels2[kp_i:kp_f+1],Fontsize= 12)
        else:
            ax.set_xticks(kpt_bounds)



    ax.yaxis.set_minor_locator(mpl.ticker.AutoMinorLocator(2))
    ax.axhline(y=0, xmax=1, linestyle='dotted', color='black', linewidth=0.5)
    plt.tight_layout(pad=1.20)
    if title != None:
        plt.title(title)
    plt.plot(dpi=dpi)
    if fname != None:
        plt.savefig(fname, dpi=dpi)

def plot_custom_vaspbands(outcar = "OUTCAR",kpoints = "KPOINTS",figsize = (10,7.5),ewindow = (-3,3),dpi = 500,linewidth = 0.5,kp_i=None,kp_f=None,title=None,fname=None):

    try:    
        kpt_path, bands, Efermi, kpt_bounds, _ ,nelect = get_bandInfo1(outcar,kpoints)
    except IndexError:
        kpt_path, bands, Efermi, kpt_bounds, _, nelect = get_bandInfo2(outcar,kpoints)
    
    bandplot2(kpt_path, bands, Efermi, kpt_bounds, nelect, figsize[0],figsize[1],ewindow[0],ewindow[1],dpi,linewidth,kp_i,kp_f,title,fname)
    
def tb_to_hr(file = "wannier90_tb.dat"):
    
    """ Generates a seedname_hr.dat file from seedname_tb.dat
        Intended to generate only _tb.dat and later using it for
        WanierBerri & processing w/o restarting wannier90.x calc """
    
    import datetime as dt 
    
    seedname = file[0:-7]
    file_tb = [line for line in open(file)]
    num_wann = int(file_tb[4])
    num_rpts = int(file_tb[5])

    for (ii,line) in enumerate(file_tb):
        if ii>1 and len(line.split()) == 3 and file_tb[ii-1] == '\n':
            loc = ii
            break

    newfile = []; newfile.append("Generated at "+dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')+" from "+seedname+"_tb.dat\n")
    newfile = newfile + file_tb[4:loc+(num_wann**2+2)*num_rpts-1]
    filename = seedname + "_hr.dat"
    
    with open(filename,"w") as f:
        for line in newfile:
            f.write(line) 
       
