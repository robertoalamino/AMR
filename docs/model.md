# Model and output conventions

Both programs use a 50 by 50 periodic lattice and genotype dimension 50.
Corrected reproduction selects uniformly among currently empty neighbours.
Each initial parent is considered once per reproduction round. Culling
independently removes each cell present at the start of a pass with
probability 0.5; it does not impose an exact quota.

Lower antimicrobial death probability q indicates greater resistance under
the model exposure. The reference calibration q = 0.3 does not set the
living-population mean to 0.3. The reversal reference is calculated at
onset from stored genotypes over all sites, including empty sites; this
calculation is preserved and is not an occupied-cell reference.

Treatment acts on [start, end). Susceptibility is measured before culling;
extinction time records the first step ending empty. In extinction outputs,
0 means initially empty and -1 means no extinction observed by the run end.
After extinction, susceptibility carries the last valid observation; it is
not a living-population measurement. With no valid observation it is NaN.

Reversal population rows cover steps 0 through T, and susceptibility rows
cover start through T. Reversal extinction rows contain realisation and
time. Sweep outcome rows contain alpha, beta, realisation, final recorded
q and extinction time. Sweep extinction rows omit q. Sweep aggregate rows
contain alpha, beta and the strict mean over all realisations, including
carried values; unavailable values propagate.

The random generator is initialised once per execution. Realisations use
the continuing stream. Equal seeds across different algorithms do not
imply paired trajectories.
