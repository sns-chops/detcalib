# Overview

The procedure includes
* Use Si+V data to calibrate long packs
* Use C60+V data to calibrate short packs

The coordinate system
* z: beam
* y: vertical up

The basic idea is
* Use white beam powder data (Si/C60) to obtain difc
  - Compute initial difc from nominal geometrical info of det packs and save it
  - Obtain I(tof) spectrum for each pixel. see for example, Si-generate-I_tof_spectra.ipynb
  - For each pixel, fit I(tof) to peaks with a scale factor (difc) to convert to d-spacing.
* Use V powder data to obtain L2
  - For each pixel, compute tof of the elastic peak, then determine L2
* Use difc+L2 of all pixels in a pack (excluding bad fits etc) to optimize the
  position (x,z) and orientation of this pack (y).
* Check the calibration by
  - Visualize the det system in mantid
  - Compute I(d) spectrum of a pack


# Use vanadium powder data to obtain L2

See [V-L2 notebook](./V-L2.ipynb)

# Nominal difc

See [notebook](nominal_difc.ipynb) or [py module](lib/nominal_difc.py).

It is also available in [get_difc_from_Itof](lib/get_difc_from_Itof.py).

# I(tof) spectra
See [Si notebook](Si-generate-I_tof-spectra.ipynb)


# difc

For each pack, try to find the best dataset (either Si or C60) for the purpose
of calculating difc.

The [calibrate module](lib/calibrate.py) is used for calculate difc and optionally
fit position and orientation of a pack.

For each row, there is a master notebook. For example [C row](C-row/calibrate-C-row.ipynb).
Some packs need special attention because it could be hard to tell whether Si or C60 is
the better dataset to use.

C row: other than the short packs, most packs are fairly straight forward.
D row: most packs are good, but some packs on top of the forward beam are hard
B row: similar to D. packs near the forward beam are hard. Especially hard is B25 and B26.
See the Dec 29, 2017 update for more details.


## Check alignment
* using notebooks. See for example  [notebook for C25B](./check-C25B-only-C60.ipynb).


# Misc Notes

* The original positions of the short packs are very wrong. Have to make a guess and put them to reasonable positions
* Have to use all C60 data. Each one is too large to work with as a whole (out of memory). Have to slice each of
  them to ~10 slices and process.

