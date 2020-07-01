# Finite-Size Corrections of Defect Energy Levels Involving Ionic Polarization

This script calculates the finite-size corrections of total energies and single-particle energy levels involving defect states with built-in ionic polarization in supercell calculations. The method accounts on an equal footing for the screening of the electrons and of the ionic polarization charge arising from the lattice distortions. These corrections allow one to achieve accurate optical transition energies and single-particle defect levels without requiring computationally prohibitive system-size scalings. An extensive documentation of the method can be found in the paper of [Falletta *et al.*](https://journals.aps.org/prb/accepted/a307bYebYaa1f267498c8912422be5af7ddfad0fc)

## Use

To run the script, type: ``python3 finite-size-corrections-defect-levels.py infile.dat``, where ``infile.dat`` in the input file, which requires the following inputs:
* ``system``: name of the system
* ``latt_param[A,deg]``: lattice parameters in Angstrom and angles in degrees
* ``r_defect[A]``: coordinates of the defect in Angstrom
* ``E_cutoff[Ry]``: energy cutoff in Rydberg
* ``direction``: direction for the average of the electrostatic potential (``x``, ``y``, ``z``)
* ``sigma[Bohr]``: variance Gaussian for calculating the alignment term
* ``eps_0``: static dielectric constant
* ``eps_inf``: high-frequency dielectric constant
* ``alignment``: flag for calculating the alignment term (``yes``, ``no``)
* ``N_states``: number of states

The states are denoted as ``(qC,qR)``, where ``qC`` is the charge state and ``qR`` the charge generating the lattice distortions. The script gives in output the corrections for the total energy and single-particle defect level of each state. If ``alignment`` is enabled, theb the cubefiles of the electrostatic potentials associated with each defect states are required in input. A plot file is also printed to show the calculation of the alignments.

