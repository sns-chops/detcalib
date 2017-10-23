# Procedure

* Use Si+V data to calibrate long packs
* Use C60+V data to calibrate short packs

Coordinate system
* z: beam
* y: vertical up

The basic idea is
* Use white beam single crystal data (Si/C60) to obtain difc
  - Compute initial difc from nominal geometrical info of det packs
  - Obtain I(d) spectrum for each pixel using the initial geo info of det packs
  - For each pixel, fit I(d) to peaks and compare them to desired d spacing values of the crystal,
    calculate a ratio of adjustment, multiply the ratio wit the initial difc to obtain
    the calibrated difc
* Use V powder data to obtain L2
  - For each pixel, compute tof of the elastic peak, then determine L2
* Use difc+L2 of all pixels in a pack (excluding bad fits etc) to optimize the
  position (x,z) and orientation of this pack (y).
* Check the calibration by
  - Visualize the det system in mantid
  - Compute I(d) spectrum of a pack


# Vanadium powder data: L2

See [V-L2 notebook](./V-L2.ipynb)


# Long packs

## Si: difc

See [Si-difc notebook](./Si-difc-2.ipynb)

New difc and mask are saved.

## Align

See [Alignment python script](./align_longpacks.py)

# Short packs

## C60: difc

* Compute I(d) spectrum of all pixels for the short packs

[python script](./C60-I_d_shortpacks.py)

** This takes ~5 hours at ndav2! **

* Compute difc

For each pack, modify and run [this python script](./difc_shortpacks.py)

* Calibrate

For each pack, modify and run [this python script](./align_shortpacks.py)

This script generates "new.xml". Copy the content to SEQ_new.xml.fit-to-difc-and-L2_Si+C60.



# Misc Notes

* The original positions of the short packs are very wrong. Have to make a guess and put them to reasonable positions
* Have to use all C60 data. Each one is too large to work with as a whole (out of memory). Have to slice each of
  them to ~10 slices and process.

