# WASPP module
Wannier &amp; vASP Postprocessing module

## Main functions:
### Wannier90 - VASP interface preparation and evaluation.

``` plot_pdos ```

Main function for plotting the partial density of states, for desired atoms and orbitals at all the different Wyckoff positions. Usage is straightforward:

* `file = "vasprun.xml"` of your VASP run.
* `_atoms` = Atoms whose pDOS you wanna plot as a list of strings e.g. ["P","I","O","Rn","Al"].
* `_orbitals` = Orbitals whose pDOS you wanna know from _atoms as a list of strings. They can be "big" orbitals (s,p,d,f) or "small" (px,py,pz,dxy,...), e.g ["s","px","d"].
* `e_window` = Range of energies of interest.

If `_atoms` and `_orbitals` are `None`, the function will plot every atom and big orbital at every different Wyckoff position. For a better visualization each kind of big orbital is displayed with a different linestyle.

For example (Ta6Se24I2): 
```import WASPP_0.5 as wasp
wasp.plot_dos("vasprun.xml",e_window = (-7.5,4))
```

Returns:

![TaSeI pDOS](TaSeI.png)

The first tag in the legend is the atom, the second is the Wyckoff position and the third the orbital.


``` band_counter ```

Counts the numebr of bands in a energy window in the whole FBZ. It gives a good clue of how to choose the energy window. Usage is as following:

``` 
band_counter(file = "vasprun.xml", emin = 0.0, emax = 0.0)
```
* `file = "vasprun.xml"` of your VASP run.
* `emin` and `emax` are the lower and upper part of the energy window given in reference to Fermi energy.

For example:
```import WASPP_0.5 as wasp
wasp.band_counter(file = "vasprun.xml", emin = -7.0, emax = 4.0)
```
Returns:
```
Efermi = 3.01610496.
Total bands = 544.
The number of bands between -7.00 eV (-3.98 eV) and 4.00 eV (7.02 eV) is 368.
```
Which is, the Fermi energy in eV, the total number of bands of the vasp run and the number of bands in the energy window (with real energies in parenthesis for wannier90.win)

` plot_wannierbands `

Function for plotting wannier bands from `.dat` and `.gnu` files.

Usage:

```plot_wannierbands(file_dat = "wannier90_band.dat", gnu = "wannier90_band.gnu",efermi = 0.0, e_window = None, fig_size = (15,8),savename = "wannierbands.png")
```
* `file_dat`: `*_band.dat` output file from a wannier90.x run.
* `file_dat`: `*_band.gnu` output file from a wannier90.x run.
* `efermi`: Fermi energy.
* `e_window` = Energy window for the plot

It generates a `"wannierbands.png"` file.

`plot_vaspbands`

Function for plotting VASP bands from a non self-consistent calculation in a KPATH. Usage:

`plot_vaspbands(outcar = "OUTCAR", kpoints = "KPOINTS")`

* `outcar`: OUTCAR file from VASP run.
* `kpoints`: KPOINTS file from nsc VASP run (linemode expected).

`plot_comparison`

Function for comparing VASP and Wannier90 bandstructures combining the previous functions and tags. Usage:

``` plot_comparison(outcar = "OUTCAR", kpoints = "KPOINTS",file_dat = "wannier90_band.dat", gnu = "wannier90_band.gnu",efermi = 0.0, fig_size = (12,8), e_window = (-4,4),savename = "comparison.png"):
```

* `outcar`: OUTCAR file from VASP run.
* `kpoints`: KPOINTS file from nsc VASP run (linemode expected).
* `file_dat`: `*_band.dat` output file from a wannier90.x run.
* `file_dat`: `*_band.gnu` output file from a wannier90.x run.
* `efermi`: Fermi energy.
* `e_window` = Energy window for the plot

Example (badly frozen e_win chosen in NbGe2):
 ``` plot_comparison(outcar = "OUTCAR", kpoints = "KPOINTS",file_dat = "wannier90_band.dat", gnu = "plottt/wannier90_band.gnu",efermi = 0.0, fig_size = (12,8), e_window = (-4,4),savename = "comparison.png") ```

![NbGe2 VASP vs Wannier90](NbGe2.png)

### MBJ and PBE potentials bandstructure comparison.


