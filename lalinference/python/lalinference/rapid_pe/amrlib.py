import itertools
import copy

import numpy
import h5py

import lal

from . import lalsimutils

m1m2 = numpy.vectorize(lalsimutils.m1m2)

#
# Utility functions
#

OFFSET_VECTORS = (1, -1, 0)
def ndim_offsets(ndim):
    for off in itertools.product(*numpy.tile(OFFSET_VECTORS, (ndim, 1))):
        yield numpy.array(off)

#
# Refinement strategies
#

# TODO: Optionally return cells
def refine_regular_grid(grid, grid_spacing, return_cntr=False):
    """
    Given a regular grid in N dimensions with spacing given by grid_spacing, refine the grid by bisection along all dimensions. Returns the new grid and new spacing. The grid will not contain the center point by default (since it is assumed the point already exists in the unrefined grid, but this behavior can be changed by setting return_cntr to True.

    >>> region = Cell(numpy.array([(-1., 1.), (-2., 2.)]))
    >>> grid, spacing = create_regular_grid_from_cell(region, side_pts=2)
    >>> grid, spacing = refine_regular_grid(grid, spacing)
    >>> print(grid)
    [array([-2., -4.]), array([-2., -2.]), array([-1., -4.]), array([-1., -2.]), array([-2.,  0.]), array([-2.,  2.]), array([-1.,  0.]), array([-1.,  2.]), array([ 0., -4.]), array([ 0., -2.]), array([ 1., -4.]), array([ 1., -2.]), array([ 0.,  0.]), array([ 0.,  2.]), array([ 1.,  0.]), array([ 1.,  2.])]
    >>> print(spacing)
    [ 1.  2.]
    """

    # FIXME: We can probably allocate space ahead of time
    new_pts = []
    for cells in grid_to_cells(grid, grid_spacing):
        for cell in cells.refine_full(return_cntr=return_cntr):
            new_pts.append(cell._center)

    return new_pts, grid_spacing / 2

# Take midpoint between target point and neighbor
def midpoint(pt1, pt2):
    diff = pt2 - pt1
    mid = diff / 2
    return pt1 + mid

class Cell(object):
    def __init__(self, boundaries, center=None):
        self._bounds = boundaries.copy()
        if center is not None:
            self._center = center.copy()
            assert len(self._center) == len(self._bounds)
            assert all([b[0] < c < b[1] for c, b in zip(self._center, self._bounds)])
        else:
            #self._center = [(a+b)/2.0 for a, b in self._bounds]
            self._center = (self._bounds[:,0] + self._bounds[:,1])/2.0

    def area(self):
        return numpy.abs(numpy.diff(self._bounds)).prod()

    def __intersects_1d(self, s1, s2):
        return s2[1] <= s1[0] > s2[0] or s2[1] < s1[1] > s2[0]

    def intersects(self, other):
        assert self._bounds.shape == other._bounds.shape
        return all(self.__intersects_1d(s1, s2) for s1, s2 in zip(self._bounds, other._bounds))

    def divide_if(self, testfunc, depth=1):
        """
        Subdivide once along each dimension, recursively.
        """
        if depth < 1:
            raise ArgumentError("Depth value must be greater than 0")

        # Base case: We're at the finest resolution, divide and check
        if depth == 1:
            cells = self.divide()
            if testfunc(cells):
                return cells
            else:
                return self
        else:

            daughters, divided = [], False
            for d in self.divide():
                # Recurse downward and get result which is returned
                dcells = d.divide_if(testfunc, depth-1)

                # We got the cell back, it didn't divide. So we preserve it to
                # check if the parent cell should divide or not
                if dcells == d:
                    daughters.append(d)

                # We got a set of divided cells, so we keep to pass back up
                else:
                    daughters.extend(dcells)
                    divided = True

            # No child (or below) of this cell was divided, pass ourself back
            # up
            if not divided and not testfunc(daughters):
                return self

        return daughters

    def refine_full(self, return_cntr=True):
        """
        Refine each cell such that new cells are created along all permutations of a displacement vector with unit vectors in each dimension and the 0 vector.
        """
        # FIXME: needs "return center" argument
        cells = []
        extent = numpy.diff(self._bounds).flatten() / 2
        # Iterate through all possible offsets for this dimension -- creating a
        # new cell and translating it along this offset
        # Note the factor of two required in the offset is preincluded above.
        for offset in ndim_offsets(len(self._center)):
            if not return_cntr and numpy.abs(offset).sum() == 0.0:
                continue
            cell = Cell(numpy.copy(self._bounds))
            # Shrink cell
            cell._bounds -= cell._center[:,numpy.newaxis]
            cell._bounds /= 2
            cell._bounds += cell._center[:,numpy.newaxis]
            # translate to new offset
            cell.translate(offset * extent)
            cells.append(cell)

        return cells

    def translate(self, offset):
        """
        Translate this cell's position and extent by offset.
        """
        self._bounds += offset[:,numpy.newaxis]
        self._center += offset

    def refine(self):
        """
        Refine each cell such that 2*dim new cells are created. In contrast to divide, refine will place the center of each new cell at the midpoint along the boundary of the mother cell. Its dimensions will be the half width from the center to the edge on each side, except for the dimension along the splitting axis, where the dimension will be twice the half width of the original cell in that dimension (and on that side).
        """
        cells = []
        dim = len(self._bounds)

        cntrs = []
        # New centers are bounds of previous cell
        for i in range(dim):
            # Split left
            cntr = numpy.copy(self._center)
            cntr[i] = self._bounds[i][0]
            cntrs.append(cntr)

            # Split right
            cntr = numpy.copy(self._center)
            cntr[i] = self._bounds[i][1]
            cntrs.append(cntr)

        # Determine length of boundaries
        bounds = []
        for i, (lb, rb) in enumerate(self._bounds):
            half_left = (self._center[i] - lb) / 2.0
            half_right = (rb - self._center[i]) / 2.0
            bounds.append( (half_left, half_right) )
        bounds = numpy.array(bounds)

        # Create each new cell
        for i, ci in enumerate(cntrs):
            cbnds = list(copy.copy(cntr))
            for j in range(len(bounds)):
                if i == j:
                    cbnds[j] = (ci[j] - bounds[j,i%2], ci[j] + bounds[j,i%2])
                else:
                    cbnds[j] = (ci[j] - bounds[j,0], ci[j] + bounds[j,1])
            cells.append(Cell(numpy.array(cbnds), numpy.array(ci)))

        return cells

    def divide(self):
        """
        Subdivide once along each dimension, recursively.
        """
        return self.__recursive_divide(len(self._bounds)-1)

    def __recursive_divide(self, dim):
        """
        Do not call directly!
        """

        if dim > 0:
            cells = self.__recursive_divide(dim-1)
        else:
            cells = [self]

        divided_cells = []
        for cell in cells:
            d1, d2 = copy.deepcopy(cell), copy.deepcopy(cell)
            # Divide left
            d1._bounds[dim][1] = self._center[dim]
            d1._center[dim] = (self._bounds[dim][0] + self._center[dim]) / 2.0
            # Divide right
            d2._bounds[dim][0] = self._center[dim]
            d2._center[dim] = (self._bounds[dim][1] + self._center[dim]) / 2.0
            divided_cells.extend([d1, d2])

        return divided_cells

    @staticmethod
    def make_cell_from_boundaries(inpt_pt, pts, symmetric=True):
        """
        Construct a 'virtual cell' from the extent of the points, centered on inpt_pt. If symmetric is True, the boundaries will be symmetric about the input point, if False, the cell will follow the extent of the points on either side.
        """
        cell_bounds = []

        pts = numpy.array(pts)

        # Find the extent of points in each dimension
        ext_right = numpy.max(pts, axis=0)
        ext_left = numpy.min(pts, axis=0)
        if symmetric:
            for el, er, pt in zip(ext_left, ext_right, inpt_pt):
                max_ext = max(abs(pt - el), abs(pt - er))
                bound = (pt - max_ext, pt + max_ext)
                cell_bounds.append(bound)
        else:
            for el, er in zip(ext_left, ext_right):
                bound = (el, er)
                cell_bounds.append(bound)

        return Cell(numpy.array(cell_bounds), inpt_pt)

class GriddingException(Exception):
    def __init__(self, *args, **kwargs):
        super(Exception, self).__init__(*args, **kwargs)

#
# Gridding / cell utilities
#
def grid_to_indices(pts, region, grid_spacing, check=True):
    """
    Convert points in a grid to their 1-D indices according to the grid extent in region, and grid spacing. If 'check' is True, ensure that all points exist within the expanded region before indexing --- this tends to avoid more cryptic errors about indexing later on.
    """
    region = region.copy() # avoid referencing
    # FIXME: This is not really exact. Because the information about how much
    # the initial region is expanded by the refinement is lost, we'll assume
    # our region extends in 3x in all directions.
    # This shouldn't affect anything since the absolute index relative to the
    # original region doesn't matter so much as the relative position
    region[:,0] -= numpy.diff(region).flatten()
    region[:,1] += numpy.diff(region).flatten()
    if check:
        if not ((region[:,0] < pts).all() and (region[:,1] > pts).all()):
            raise GriddingException("Some or all of provided points are not within the region. Are your dimension labels swapped?")

    extent = numpy.diff(region)[:,0]
    pt_stride = numpy.round(extent / grid_spacing).astype(int)
    # Necessary for additional point on the right edge
    pt_stride += 1
    #print pt_stride
    idx = numpy.round((pts - region[:,0]) / grid_spacing).astype(int)
    indices = []
    for i in idx:
        indices.append(numpy.ravel_multi_index(i, pt_stride))
        #print i, indices[-1]
    return numpy.array(indices)

def prune_duplicate_pts(pts, region, grid_spacing):
    """
    Remove identical points from list.
    """
    ind = grid_to_indices(pts, region, grid_spacing)
    _, ind = numpy.unique(ind, return_index=True)
    return numpy.array(pts)[ind]

def create_regular_grid_from_cell(cell, side_pts=5, return_cells=False):
    """
    Given a region (amrlib.Cell), grid it with side_pts to a dimension. If return_cells is given, the divided cells will be returned rather than the grid points.

    >>> import numpy
    >>> region = Cell(numpy.array([(-1., 1.), (-2., 2.)]))
    >>> grid, spacing = create_regular_grid_from_cell(region, side_pts=2)
    >>> print(grid)
    [[-1. -2.]
     [-1.  2.]
     [ 1. -2.]
     [ 1.  2.]]
    >>> print(spacing)
    [ 2.  4.]
    """

    # TODO Recenter on trigger point

    # This creates the logic to create points along a given axis
    grid_pts = [slice(ls, rs, side_pts * 1j) for ls, rs in cell._bounds]
    # This produces a representation of the points, but not in a tuple form
    grid_pts = numpy.mgrid[grid_pts]
    # This tuplizes them
    grid_pts = numpy.array(zip(*[grid.flatten() for grid in grid_pts]))

    # useful knowledge for later
    grid_spacing = numpy.array([(rb - lb) / (side_pts - 1) for lb, rb in cell._bounds])
    return grid_pts, grid_spacing

# FIXME: This isn't currently used anywhere, unsure of it's usefulness
def create_new_grid(cell, pts, max_level=6):
    """
    Subdivide until all cells are self contained - e.g. every cell contains no more than one pt from pts. If the number of divisions exceeds max_level, stop dividing.
    """
    cells = [cell]
    i = 1
    while True:
        #print "Division level %d, ncells: %d" % (i, len(cells))
        new_cells = []
        for c in cells:
            new_cells.extend(c.divide())

        cell_centers = numpy.array([c._center for c in new_cells])
        cell_tree = BallTree(cell_centers)
        cells_occupied = set()
        for j, pt in enumerate(pts[idx]):
            cell_idx = cell_tree.query(pt, k=1, return_distance=False)[0]
            selected_cell = tuple(cell_centers[cell_idx][0])
            #print "Template falls in cell centered on %s. %d / %d pts checked" % (str(selected_cell), j+1, len(pts[idx]))
            if selected_cell not in cells_occupied:
                cells_occupied.add(selected_cell)
            else:
                break

        if len(cells_occupied) == len(pts):
            break

        cells = new_cells
        if i >= max_level:
            break
        i += 1
    return cells

# FIXME: Make into a generator
def grid_to_cells(grid, grid_spacing):
    return [Cell(
            numpy.array([(c-sp/2, c+sp/2) for c, sp in zip(pt, grid_spacing)]),
            pt) for pt in grid]

#
# Cell grid serialization and packing routines
#

def init_grid_hdf(init_region, h5file, overlap_thresh, crd_sys, intr_prms=None, base_grp="rapidpe_grids"):
    """
    Set up a new HDF5 file (h5file), truncating any existing file with this name. A new 'folder' called 'base_grp' is set up and the initial region with attribute 'overlap_thresh' is set up under it.
    """
    hfile = h5py.File(h5file, "w")
    if base_grp not in hfile:
        hfile.create_group(base_grp)
    hfile[base_grp].create_dataset("init_region", data=init_region._bounds)
    hfile[base_grp].attrs.create("overlap_thresh", overlap_thresh)
    hfile[base_grp].attrs.create("distance_coordinates", crd_sys)
    if intr_prms is not None:
        hfile[base_grp].attrs.create("labels", intr_prms)
    return hfile[base_grp]

def save_grid_cells_hdf(base_grp, cells, crd_sys, intr_prms=None, check=True):
    """
    Under the base_grp, a new level of grid points is saved. It will create a subgroup called "grids" if it does not already exist. Under this subgroup, levels are created sequentially, the function looks for the last level (e.g. the subsubgroup with the largest "level" attribute, and appends a new level with +1 to that number. If check is enabled (that is the default), a safety check against the new level versus the last level is made to ensure the resolution is a factor of two smaller.
    """
    if "grids" not in base_grp:
        grids = base_grp.create_group("grids")
        grids.attrs.create("distance_coordinates", crd_sys)
    else:
        grids = base_grp["grids"]
        assert grids.attrs["distance_coordinates"] == crd_sys

    levels = []
    for name, ds in grids.iteritems():
        levels.append((ds.attrs["level"], name))

    grid_res = numpy.diff(cells[0]._bounds).flatten()
    if len(levels) == 0:
        ds = grids.create_dataset("level_0", data=numpy.array([c._center for c in cells]))
        ds.attrs.create("level", 0)
        ds.attrs.create("resolution", grid_res)
        lvl = 0
    else:
        levels.sort()
        lvl, name = levels[-1]
        lvl += 1
        if check:
            assert numpy.allclose(grids[name].attrs["resolution"] / 2, grid_res)

        ds = grids.create_dataset("level_%d" % lvl, data=numpy.array([c._center for c in cells]))
        ds.attrs.create("level", lvl)
        ds.attrs.create("resolution", grid_res)

    if intr_prms is not None:
        assert len(intr_prms) == ds.shape[1]
        ds.attrs.create("labels", intr_prms)

    return lvl

def load_init_region(h5file, get_labels=False, base_grp="rapidpe_grids"):
    """
    Load the initial region for a set grid points (in the form of cells) from h5file.
    """
    hfile = h5py.File(h5file, "r")
    if get_labels:
        if "labels" in hfile[base_grp].attrs:
            labels = hfile[base_grp].attrs["labels"]
        else:
            labels = tuple()

    if get_labels:
        return Cell(hfile[base_grp]["init_region"][:]), labels
    else:
        return Cell(hfile[base_grp]["init_region"][:])

def load_grid_level(h5file, level, return_labels= False, base_grp="rapidpe_grids"):
    """
    Load a set grid points (in the form of cells) from h5file. If level is None, return the base_grp, if level is -1, return the highest resolution available, otherwise an ArgumentError is raised if 'level' is not represented in the attributes of one of the levels stored under base_grp/'grids'
    """
    hfile = h5py.File(h5file, "a")
    if level is None:
        return hfile[base_grp]

    grids = hfile[base_grp]["grids"]
    if level == -1:
        level = sorted(dat.attrs["level"] for name, dat in grids.iteritems())[-1]

    for name, dat in grids.iteritems():
        if dat.attrs["level"] == level:
            grid_res = dat.attrs["resolution"][numpy.newaxis,:]
            if return_labels:
                labels = dat.attrs["labels"] if "labels" in dat.attrs else tuple()
                return unpack_grid_cells(numpy.concatenate((grid_res, dat[:]))), level, labels
            else:
                return unpack_grid_cells(numpy.concatenate((grid_res, dat[:]))), level
    raise ArgumentError("No grid refinement level %d" % level)

def pack_grid_cells(cells, resolution):
    """
    Pack the cell centers (grid points) into a convienient numpy array. The resolution of the grid is prepended for later use.
    """
    return numpy.vstack((resolution, [c._center for c in cells]))

def unpack_grid_cells(npy_cells):
    """
    Unpack the cell centers (grid points), the cell collection is returned.
    """
    res = npy_cells[0]
    cells = []
    for cntr in npy_cells[1:]:
        cell = Cell(numpy.array((cntr - res / 2, cntr + res / 2)).T)
        cells.append(cell)
    return cells, res

def serialize_grid_cells(level_dict, fname):
    """
    Given a dictionary of refinement levels, pack it up and save it to a numpy file.
    FIXME: This can't yet distinguish the level of refinement. Need to teach it to pay attention to the resolution.
    """
    cell_block = numpy.array([level_dict[i] for i in sorted(level_dict)])
    numpy.save(fname, cell_block)

def deserialize_grid_cells(fname):
    """
    Unpacks a numpy file of level grid points and returns a cell collection.
    FIXME: This can't yet distinguish the level of refinement. Need to teach it to pay attention to the resolution.
    """
    cell_block = numpy.load(fname)
    cells = []
    for level in cell_block:
        level_cells, res = unpack_grid_cells(level)
        cells.extend(level_cells)
    return cells, res


#
# Coordinate transformations
#
def transform_m1m2_mceta(m1, m2):
    return lalsimutils.Mceta(m1, m2)

def transform_mceta_m1m2(mc, eta):
    return m1m2(mc, eta)

__prefac_0 = 5. / 256 / numpy.pi
__prefac_3 = 1. / 8
__dim_mass = lal.G_SI / lal.C_SI**3 * lal.MSUN_SI
def transform_m1m2_tau0tau3(m1, m2, flow=40.):
    mt = m1 + m2
    eta = m1 * m2 / mt**2
    mt *= numpy.pi * flow * __dim_mass
    return (__prefac_0 / flow / eta * mt**(-5./3), __prefac_3 / flow / eta * mt**(-2./3))

__prefac_tau = 5. / 32 / numpy.pi
def transform_tau0tau3_m1m2(tau0, tau3, flow=40.):
    mt = __prefac_tau / numpy.pi / flow * tau3 / tau0
    eta = 1.0 / 8 / flow / tau3 * (tau0 / __prefac_tau / tau3)**(2./3)
    m1, m2 = m1m2(mt*eta**(3./5), eta)
    return m1 / __dim_mass, m2 / __dim_mass

def transform_s1zs2z_chi(m1, m2, s1z, s2z):
    return (m1 * s1z + m2 * s2z) / (m1 + m2)

def transform_m1m2_mcq(m1, m2):
    mc = (m1 * m2)**(3./5) / (m1 + m2)**(1./5)
    q = np.min([m1, m2], axis=0) / np.max([m1, m2], axis=0)
    return mc, q

def transform_mcq_m1m2(mc, q):
    m = mc * (1 + q)**(1./5)
    return m / q**(3./5), m * q**(2/5.)

#
# Coordinate transformation boundaries
#

__a0 = __prefac_0 * numpy.pi
__a3 = __prefac_3 * numpy.pi
def check_tau0tau3(tau0, tau3, flow=40):
    A_0 = __a0 / (numpy.pi * flow)**(8./3)
    A_3 = __a3 / (numpy.pi * flow)**(5./3)

    # Lower bound (\eta > 1/4)
    tau3_bnd = 4 * A_3 * (tau0 / 4 / A_0)**(2.0/5)

    # Note, there are two other bounds, but we don't care about them so much
    return tau3 > tau3_bnd

def check_mchirpeta(mchirp, eta):
    return (numpy.array(eta) <= 0.25) & (numpy.array(eta) > 0) & (numpy.array(mchirp) >= 0)

def check_spins(spin):
    return numpy.sqrt(numpy.atleast_2d(spin**2).sum(axis=0)) <= 1

# Make sure the new grid points are physical
def check_grid(grid, intr_prms, distance_coordinates):
    """
    Check the validity of points in various coordinate spaces to ensure they are physical.
    """
    m1_axis, m2_axis = intr_prms.index("mass1"), intr_prms.index("mass2")
    grid_check = numpy.array(grid).T
    if distance_coordinates == "tau0_tau3":
        bounds_mask = check_tau0tau3(grid_check[m1_axis], grid_check[m2_axis])
    elif distance_coordinates == "mchirp_eta":
        bounds_mask = check_mchirpeta(grid_check[m1_axis], grid_check[m2_axis])

    # FIXME: Needs general spin
    if "spin1z" in intr_prms:
        s1_axis = intr_prms.index("spin1z")
        bounds_mask &= check_spins(grid_check[s1_axis])
    if "spin2z" in intr_prms:
        s2_axis = intr_prms.index("spin2z")
        bounds_mask &= check_spins(grid_check[s2_axis])
    if "chi_z" in intr_prms:
        chi_axis = intr_prms.index("chi_z")
        bounds_mask &= check_spins(grid_check[chi_axis])
    return bounds_mask

VALID_TRANSFORMS_MASS = { \
    "mchirp_eta": transform_m1m2_mceta,
    "mchirp_q": transform_m1m2_mcq,
    "tau0_tau3": transform_m1m2_tau0tau3,
    None: None
}

INVERSE_TRANSFORMS_MASS = { \
    transform_m1m2_mceta: transform_mceta_m1m2,
    transform_m1m2_mcq: transform_mcq_m1m2,
    transform_m1m2_tau0tau3: transform_tau0tau3_m1m2,
    None: None
}

def apply_transform(pts, intr_prms, mass_transform=None, spin_transform=None):
    # You know what... recarrays are dumb, and so's your face.
    # FIXME: Why does numpy want me to repack this into tuples!?
    #tpts = numpy.array([tuple(pt) for pt in pts.T], dtype = numpy.dtype([(a ,"float64") for a in intr_prms]))
    m1_idx, m2_idx = intr_prms.index("mass1"), intr_prms.index("mass2")
    if spin_transform:
        if spin_transform == "chi_z":
            s1z_idx, s2z_idx = intr_prms.index("spin1z"), intr_prms.index("spin2z")
            chi_z = transform_s1zs2z_chi(pts[:,m1_idx], pts[:,m2_idx], pts[:,s1z_idx], pts[:,s2z_idx]) 
            pts = numpy.vstack((pts.T, chi_z)).T
            intr_prms.append("chi_z")

    if mass_transform:
       pts[:,m1_idx], pts[:,m2_idx] = VALID_TRANSFORMS_MASS[mass_transform](pts[:,m1_idx], pts[:,m2_idx])

    # Independent transforms go here

    return pts

def apply_inv_transform(pts, intr_prms, mass_transform=None):
    m1_idx, m2_idx = intr_prms.index("mass1"), intr_prms.index("mass2")
    if mass_transform:
       pts[:,m1_idx], pts[:,m2_idx] = INVERSE_TRANSFORMS_MASS[VALID_TRANSFORMS_MASS[mass_transform]](pts[:,m1_idx], pts[:,m2_idx])

    # Independent transforms go here

    return pts
