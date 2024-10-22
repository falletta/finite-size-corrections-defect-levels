# Finite-Size Corrections of Defect Energy Levels Involving Ionic Polarization


<img align="right" src="logo.png" width="20%">

This script calculates the finite-size corrections of total energies and single-particle energy levels involving defect states with built-in ionic polarization in supercell calculations. The method accounts on an equal footing for the screening of the electrons and of the ionic polarization charge arising from the lattice distortions. These corrections allow one to achieve accurate optical transition energies and single-particle defect levels without requiring computationally prohibitive system-size scalings. An extensive documentation of the method can be found in the paper of [Falletta *et al.*](https://journals.aps.org/prb/accepted/a307bYebYaa1f267498c8912422be5af7ddfad0fc)

## Use

To run the script, type: ``python3 finite-size-corrections-defect-levels.py infile.dat``, where ``infile.dat`` takes the following inputs:
* ``system``: name of the system
* ``latt_param[A,deg]``: lattice parameters in Angstrom and angles in degrees
* ``r_defect[A]``: coordinates of the defect in Angstrom
* ``E_cutoff[Ry]``: energy cutoff in Rydberg
* ``direction``: direction for the average of the electrostatic potential (``x``, ``y``, ``z``)
* ``sigma[Bohr]``: variance Gaussian for calculating the alignment term in Bohr
* ``eps_0``: static dielectric constant
* ``eps_inf``: high-frequency dielectric constant
* ``alignment``: flag for calculating the alignment term (``yes``, ``no``)
* ``N_states``: number of states

The states are denoted as ``(qC,qR)``, where ``qC`` is the charge state and ``qR`` the charge generating the lattice distortions. The script gives in output the corrections for the total energy and single-particle defect level of each state. If ``alignment`` is enabled, the cubefiles of the electrostatic potentials of the states ``(qC,qR)``, ``(qR,qR)`` and ``(0,0)`` are required in input. The calculation of the alignment terms is plotted in a pdf file. A visual inspection of this file is recommended to check the correct functioning of the code. In particular, the DFT defect potential is smoothened with a window function to reduce the noise. One might change the parameters of the window function to ensure a correct noise reduction.

Considering that the cubefiles could be written in different units depending on the code used, the constant ``Eunit2eV`` must be set to allow for the correct conversion to eV. For instance, if one uses the cubefiles printed by the code CP2K, the potentials are expressed in Hartree and therefore ``Eunit2eV`` must convert Hartree to eV. At variance, if one uses the code Quantum Espresso, the potentials are expressed in Rydberg and therefore ``Eunit2eV`` must convert Rydberg to eV.   

An example of input file is ``input-MgO.dat``, where the finite-size corrections are calculated in the cases of the charged and neutral hole polaron in MgO (supercell 64 atoms). The complete set of data used in the paper of [Falletta *et al.*](https://journals.aps.org/prb/accepted/a307bYebYaa1f267498c8912422be5af7ddfad0fc) can be found in the [Materials Cloud repository](https://archive.materialscloud.org/record/446).

## Cite
Please cite [our paper](https://journals.aps.org/prb/accepted/a307bYebYaa1f267498c8912422be5af7ddfad0fc) if you use this code in your own work.
```
@article{PhysRevB.102.041115,
  title = {Finite-size corrections of defect energy levels involving ionic polarization},
  author = {Falletta, Stefano and Wiktor, Julia and Pasquarello, Alfredo},
  journal = {Phys. Rev. B},
  volume = {102},
  issue = {4},
  pages = {041115},
  numpages = {5},
  year = {2020},
  month = {Jul},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.102.041115},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.102.041115}
}
```

