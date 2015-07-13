import itertools
import copy
import numpy
import lal
from lalinference.rapid_pe import lalsimutils

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

# Take midpoint between target point and neighbor
def midpoint(pt1, pt2):
    diff = pt2 - pt1
    mid = diff / 2
    return pt1 + mid

class Cell(object):
    def __init__(self, boundaries, center=None):
        self._bounds = boundaries
        if center is not None:
            self._center = center
            assert len(self._center) == len(self._bounds)
            assert all([b[0] < c < b[1] for c, b in zip(self._center, self._bounds)])
        else:
            self._center = [(a+b)/2.0 for a, b in self._bounds]

    def area(self):
        return numpy.abs(numpy.diff(self._bounds)).prod()

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

    def refine2(self):
        cells = []
        extent = numpy.diff(self._bounds).flatten() / 2
        # Iterate through all possible offsets for this dimension -- creating a
        # new cell and translating it along this offset
        # Note the factor of two required in the offset is preincluded above.
        for offset in ndim_offsets(len(self._center)):
            cell = Cell(numpy.copy(self._bounds))
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

#
# Utilities
#
def grid_to_indices(pts, region, grid_spacing):
    """
    Convert points in a grid to their 1-D indices according to the grid extent in region, and grid spacing.
    """
    region[:,0] = region[:,0] - grid_spacing
    region[:,1] = region[:,1] + grid_spacing
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
    return eta <= 0.25 and mchirp >= 0

VALID_TRANSFORMS_MASS = { \
    "mchirp_eta": transform_m1m2_mceta,
    "tau0_tau3": transform_m1m2_tau0tau3,
    None: None
}

INVERSE_TRANSFORMS_MASS = { \
    transform_m1m2_mceta: transform_mceta_m1m2,
    transform_m1m2_tau0tau3: transform_tau0tau3_m1m2,
    None: None
}

def apply_transform(pts, intr_prms, mass_transform=None):
    # You know what... recarrays are dumb, and so's your face.
    # FIXME: Why does numpy want me to repack this into tuples!?
    #tpts = numpy.array([tuple(pt) for pt in pts.T], dtype = numpy.dtype([(a ,"float64") for a in intr_prms]))
    m1_idx, m2_idx = intr_prms.index("mass1"), intr_prms.index("mass2")
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
