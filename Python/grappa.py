"""
GRAPPA 

Mark Chiew (mark.chiew@utoronto.ca)
Python port of https://github.com/mchiew/grappa-tools/tree/main/grappa

"""
import numpy as np

def grappa(data, calib, R, kernel, delta=0):
    """
    inputs: 
        data    -   [nc, nx, ny, nz] shaped array of complex undersampled k-space data
        calib   -   [nc, cx, cy, cz] shaped array of complex fully-sampled calibration k-space data
        R       -   (1, Ry) or (1, Ry, Rz) acceleration factors tuple
                    requires that the first spatial dimension (x) be fully sampled
        kernel  -   (kx, ky) or (kx, ky, kz) kernel size tuple
        delta   -   optional, specifies ∆ CAIPI shift (i.e. ∆ky shift for successive kz planes)
                    if not specified, defaults to 0 (no CAIPI shift)

    output:
        recon   -   [nc, nx, ny, nz] shaped array of complex reconstructed k-space data


    example usage:
        # for a 2D dataset with shape [nc, nx, ny], Ry=2, kernel 3 points in kx, 2 points in ky
        recon = grappa.grappa(data, calib, (1,2), (3,2))
    """

    # If 2D inputs, add 3rd dimension
    if data.ndim == 3:
        data = data[:,:,:,None]
    if calib.ndim == 3:    
        calib = calib[:,:,:,None]
    if len(R) == 2:
        R = R + (1,)
    if len(kernel) == 2:
        kernel = kernel + (1,)

    # generate different kernel geometries
    kernels = build_kernels(R, kernel, delta)        

    # pad and initialize output
    recon, pad = __pad(data, kernels)

    # find sampling mask
    sampling_mask = recon[0] != 0

    # loop over possible kernel types
    for ker in kernels:

        # get src, trg points for training
        src_idx, trg_idx = get_indices(calib.shape, ker)
        src = calib[src_idx]
        src = src.reshape(-1, src.shape[-1])
        trg = calib[trg_idx]

        # fit/estimate weights
        # W = trg @ np.linalg.pinv(src) 
        W = np.linalg.lstsq(src.T, trg.T, rcond=None)[0].T

        # get src, trg points for reconstruction
        src_idx, trg_idx = get_indices(recon.shape, ker, mask=sampling_mask)
        src = recon[src_idx]
        src = src.reshape(-1, src.shape[-1])
        
        # apply weights to reconstruct missing data 
        recon[trg_idx] = W @ src

    # unpad and return recon
    return __unpad(recon, pad)
    

def __pad(data, kernels):
    """
    Zero pad data around edges
    """
    
    all_offsets = np.concatenate([k.source_offsets for k in kernels], axis=0)
    pad = tuple(int(np.abs(all_offsets[:, i]).max()) for i in range(3))
    return np.pad(data, ((0, 0), (pad[0], pad[0]), (pad[1], pad[1]), (pad[2], pad[2]))), pad

def __unpad(data, pad):
    """
    Remove padding applied by __pad()
    """

    if data.shape[3] == 1:
        return data[:, pad[0]:-pad[0], pad[1]:-pad[1], 0]
    else:
        return data[:, pad[0]:-pad[0], pad[1]:-pad[1], pad[2]:-pad[2]]
    
class GrappaKernel:
    """
    One kernel type corresponding to a particular missing point
    within the Ry x Rz acceleration cell.
    """

    def __init__(self, source_offsets, target_offset):
        self.source_offsets = np.asarray(source_offsets, dtype=int)
        self.target_offset = np.asarray(target_offset, dtype=int)

def build_kernels(R, kernel_size, delta):
    """
    Create all GRAPPA kernel geometries.
    """

    Rx, Ry, Rz = R
    kx, ky, kz = kernel_size

    if Rx != 1:
        raise ValueError("This implementation assumes x is fully sampled")

    source_offsets = []

    for dx in range(-((kx-1)//2), kx//2 + 1):
        for dy in range(-((ky-1)//2), ky//2 + 1):
            for i, dz in enumerate(range(-((kz-1)//2), kz//2 + 1)):
                source_offsets.append((dx * Rx, dy * Ry + (i*delta)%Ry, dz * Rz))

    kernels = []

    for ry in range(Ry):
        for rz in range(Rz):

            if ry == 0 and rz == 0:
                continue

            target_offset = np.array([0, ry, rz])
            kernels.append(GrappaKernel(source_offsets, target_offset))

    return kernels

def __valid_ranges(src_offsets, nx, ny, nz):
    """Given the (nOff, 3) offset array, return valid center index
    ranges [max(|offset|), n - max(|offset|)) per axis."""
    maxx = int(np.abs(src_offsets[:, 0]).max())
    maxy = int(np.abs(src_offsets[:, 1]).max())
    maxz = int(np.abs(src_offsets[:, 2]).max())
    xs = np.arange(maxx, nx - maxx)
    ys = np.arange(maxy, ny - maxy)
    zs = np.arange(maxz, nz - maxz)
    return xs, ys, zs


def get_indices(shape, kernel, mask=None):
    """
    Build fancy-index tuples for gathering source/target points for
    this kernel's geometry.

    shape : (nc, nx, ny, nz) shape of the array these indices will be
            used against.
    mask  : optional sampling mask (shape (nx,ny,nz). When given, points 
            are filtered down to those where the target is actually missing 
            When mask is None (e.g. training on fully-sampled calib
            data), no filtering is done.

    Returns
        src_idx : tuple usable as arr[src_idx] -> (n_off, nc, n_points)
        trg_idx : tuple usable as arr[trg_idx] -> (nc, n_points)
    """
    
    nc, nx, ny, nz = shape
    src_offsets = kernel.source_offsets
    trg_offset = kernel.target_offset
    n_off = src_offsets.shape[0]

    xs, ys, zs = __valid_ranges(src_offsets, nx, ny, nz)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing='ij')
    Xf, Yf, Zf = X.ravel(), Y.ravel(), Z.ravel()

    Tx = Xf + trg_offset[0]
    Ty = Yf + trg_offset[1]
    Tz = Zf + trg_offset[2]

    if mask is not None:
        keep = ~mask[Tx, Ty, Tz]

        Xf, Yf, Zf = Xf[keep], Yf[keep], Zf[keep]
        Tx, Ty, Tz = Tx[keep], Ty[keep], Tz[keep]

    n_points = Xf.size

    # one row of coordinates per source tap, shape (n_off, n_points)
    Xo = np.empty((n_off, n_points), dtype=np.intp)
    Yo = np.empty((n_off, n_points), dtype=np.intp)
    Zo = np.empty((n_off, n_points), dtype=np.intp)
    for i, (dx, dy, dz) in enumerate(src_offsets):
        Xo[i] = Xf + dx
        Yo[i] = Yf + dy
        Zo[i] = Zf + dz

    # explicit coil index (instead of a `:` slice) so the broadcast of
    # all four index arrays lands on (n_off, nc, n_points) directly
    coil = np.arange(nc)[None, :, None]
    src_idx = (coil, Xo[:, None, :], Yo[:, None, :], Zo[:, None, :])
    trg_idx = (slice(None), Tx, Ty, Tz)

    return src_idx, trg_idx