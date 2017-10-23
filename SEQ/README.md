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


# Long packs

## Si: difc

See ./Si-difc-2.ipynb

