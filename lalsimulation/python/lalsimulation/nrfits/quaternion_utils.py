import numpy as np

#-----------------------------------------------------------------------------
def multiply_quats(q1, q2):
    """q1, q2 must be [scalar, x, y, z] but those may be arrays or scalars"""
    return np.array([
            q1[0]*q2[0] - q1[1]*q2[1] - q1[2]*q2[2] - q1[3]*q2[3],
            q1[2]*q2[3] - q2[2]*q1[3] + q1[0]*q2[1] + q2[0]*q1[1],
            q1[3]*q2[1] - q2[3]*q1[1] + q1[0]*q2[2] + q2[0]*q1[2],
            q1[1]*q2[2] - q2[1]*q1[2] + q1[0]*q2[3] + q2[0]*q1[3]])

#-----------------------------------------------------------------------------
def quat_inv(q):
    """Returns QBar such that Q*QBar = 1"""
    qConj = -q
    qConj[0] = -qConj[0]
    normSqr = multiply_quats(q, qConj)[0]
    return qConj/normSqr

#-----------------------------------------------------------------------------
def align_vec_quat(vec):
    """Returns a unit quaternion that will align vec with the z-axis"""
    alpha = np.arctan2(vec[1], vec[0])
    beta = np.arccos(vec[2])
    gamma = -alpha*vec[2]
    cb = np.cos(0.5*beta)
    sb = np.sin(0.5*beta)
    return np.array([cb*np.cos(0.5*(alpha + gamma)),
                     sb*np.sin(0.5*(gamma - alpha)),
                     sb*np.cos(0.5*(gamma - alpha)),
                     cb*np.sin(0.5*(alpha + gamma))])

#-----------------------------------------------------------------------------
def transform_time_dependent_vector(quat, vec, inverse=0):
    """Given (for example) a minimal rotation frame quat, transforms
    vec from the minimal rotation frame to the inertial frame.
    With inverse=1, transforms from the inertial frame to the minimal
    rotation frame."""
    qInv = quat_inv(quat)
    if inverse:
        return transform_time_dependent_vector(qInv, vec, inverse=0)

    return multiply_quats(quat, multiply_quats(np.append(np.array([
            np.zeros(len(vec[0]))]), vec, 0), qInv))[1:]


#-----------------------------------------------------------------------------
def rotate_in_plane(chi, phase):
    """ Rotates a given vector, chi, by a clockwise angle, phase, about the
    z-axis. Can be used for transforming spins from the coprecessing frame to
    the coorbital frame"""
    v = chi.T
    sp = np.sin(phase)
    cp = np.cos(phase)
    res = 1.*v
    res[0] = v[0]*cp + v[1]*sp
    res[1] = v[1]*cp - v[0]*sp
    return res.T

#-----------------------------------------------------------------------------
def transform_vector_coorb_to_inertial(vec_coorb, orbPhase, quat_copr):
    """Given a vector (of size 3) in coorbital frame, orbital phase in
    coprecessing frame and a minimal rotation frame quat, transforms
    the vector from the coorbital to the LAL inertial frame.
    """

    # Transform to coprecessing frame
    vec_copr = rotate_in_plane(vec_coorb, -orbPhase)

    # Transform to inertial frame (for the surrogate)
    vec = transform_time_dependent_vector(np.array([quat_copr]).T,
        np.array([vec_copr]).T).T[0]

    return np.array(vec)
