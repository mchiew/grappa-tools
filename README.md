# GRAPPA Reconstruction Tools

These are MATLAB-based tools for performing 2D or 3D GRAPPA parallel imaging reconstructions, estimating GRAPPA g-factor maps, and simultaneous multi-slice unaliasing GRAPPA methods. These use and extend the same framework as the GRAPPA tutorial provided here, and should all be completely self-contained without any external dependencies. Included are:

```grappa/grappa.m```
: 2D or 3D regular GRAPPA reconstruction
    
```grappa/grappa_gfactor.m```
: Same as above, but performs unaliasing in image space, and can provide analytical g-factor noise amplification maps
    
```slice-grappa/sg.m```
: Slice-GRAPPA simultaneous multi-slice unaliasing
    
```slice-grappa/spsg.m```
: Split-Slice (LeakBlock) GRAPPA simultaneous multi-slice unaliasing
