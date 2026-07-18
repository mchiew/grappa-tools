Notes:

- kernel sizes (kx,ky) should both be odd
- if the input data is CAIPI shifted, then the calib data should be shifted to match.
- this recon only unaliases 1->N slices, it doesn't undo the FOV shifts, that has to be done separately afterwards
- split-slice GRAPPA is essentially the same thing as the "LeakBlock" option in CMRR language