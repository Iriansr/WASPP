# WASPP module
Wannier &amp; vASP Postprocessing module

## Main functions:

``` plot_pdos ```

Main function for plotting the partial density of states, for desired atoms and orbitals at all the different Wyckoff positions. Usage is straightforward:

* `file = "vasprun.xml"` of your VASP run.
* `_atoms` = Atoms whose pDOS you wanna plot as a list of strings e.g. ["P","I","O","Rn","Al"].
* `_orbitals` = Orbitals whose pDOS you wanna know from _atoms as a list of strings. They can be "big" orbitals (s,p,d,f) or "small" (px,py,pz,dxy,...), e.g ["s","px","d"].
* `e_window` = Range of energies of interest.

If `_atoms` and `_orbitals` are `None`, the function will plot every atom and big orbital at every different Wyckoff position. For a better visualization each kind of big orbital is displayed with a different linestyle.

For example: 
```import WASPP_0.5 as wasp
wasp.plot_dos("vasprun.xml",e_window = (-7.5,4))
```

Returns:

![TaSeI pDOS](TaSeI.png)
