import math

import numpy

from .Quaternions import Quaternion

thr_norm = 1.0e-6
thr_dot = 1.0e-8
thr_mat = 1.0e-10


def get_arotation_matrix(axes, angle, deg=False):
    """Get the rotation matrix around on of the axes
    Param:axes:one of(0, 1, 2) for (x, y, z)
    Param angle:0.0 -> pi counter-clockwise rotation
    Param:deg:angle is in degrees instead of radians
    """
    if deg:
        angle = math.radians(angle)
    rotmat = numpy.zeros((3, 3))
    index0 = (axes) % 3
    index1 = (axes + 1) % 3
    index2 = (axes + 2) % 3
    rotmat[index0, index0] = 1.0
    rotmat[index1, index1] = math.cos(angle)
    rotmat[index2, index2] = math.cos(angle)
    rotmat[index2, index1] = math.sin(angle)
    rotmat[index1, index2] = -math.sin(angle)
    return rotmat


def orient_from(ini_vector, fin_vector, deg=False):
    """Give the angle and axis of rotation to orient an object starting in the 'ini_vector' direction
        toward the "fin_vector" direction
    :param:ini_vector:initial vector direction
    :param:fin_vector:final vector direction
        The axis of rotation is the vector perpenticular to ini_vector and fin_vector
        angle:0.0 -> pi counter-clockwise rotation
        Param:deg:angle is in degrees instead of radians
    """
    if not isinstance(ini_vector, numpy.ndarray):
        ini_vector = numpy.array(ini_vector)
    else:
        ini_vector = ini_vector.copy()
    if not isinstance(fin_vector, numpy.ndarray):
        fin_vector = numpy.array(fin_vector)
    else:
        fin_vector = fin_vector.copy()

    norm_ini = numpy.linalg.norm(ini_vector)
    norm_fin = numpy.linalg.norm(fin_vector)

    if (norm_ini < thr_norm) or (norm_fin < thr_norm):
        return 0.0, numpy.zeros((3,))
    ini_vector /= norm_ini
    fin_vector /= norm_fin

    #    bisector = (ini_vector + fin_vector) / 2.0
    #    norm = numpy.linalg.norm(bisector)
    #    return 180.0, bisector/norm

    s = numpy.dot(ini_vector, fin_vector)
    angle = math.acos(s)
    if (1.0 - abs(s)) < thr_dot:
        dummyvector = numpy.array([1.0, 1.0, 1.0])
        rotationaxis = dummyvector - numpy.dot(ini_vector, dummyvector) * ini_vector
    else:
        rotationaxis = numpy.cross(ini_vector, fin_vector)

    norm = numpy.linalg.norm(rotationaxis)
    if norm < thr_norm:
        return 0.0, numpy.zeros((3,))
    rotationaxis /= norm
    if deg:
        angle = numpy.degrees(angle)
    return angle, rotationaxis


def get_angle_axis(rotmat, deg=False):
    """
    Get the rotation of an angle theta around an axis of rotation
    Angle:0->pi
    Param:deg:angle is in degrees instead of radians
    """
    eig, eiv = numpy.linalg.eig(rotmat)
    eig = eig.real
    index = numpy.where(abs(eig - 1.0) < thr_dot)[0]
    axes = eiv[:, index[-1]].real
    angle = numpy.arccos(eig[(index[-1] + 1) % 3].real)
    # FIXME: sign of the angle
    rot = Rotation(angle, axes[0], axes[1], axes[2])
    rotmat2 = rot.get_rotation_matrix()
    if (abs(rotmat - rotmat2) < thr_mat).all():
        return rot.get_angle(deg=deg), rot.get_axis()
    elif (abs(rotmat - rotmat2.transpose()) < thr_mat).all():
        return -rot.get_angle(deg=deg), rot.get_axis()
    else:
        raise Exception("Get_angle_axis: rotation matrix not found correctly")


class Rotation(object):
    def __init__(self, theta, l, m, n, deg=False):
        """
        Define a rotation of an angle theta around an axis of rotation
        given by l, m, n direction cosines(dot product with the
        orthogonal unit vectors)
        Param:theta(radian):angle of rotation with positive angle meaning a counter-clockwise rotation
        Param:deg:theta is in degrees instead of radians
        """
        self.theta = theta
        if deg:
            self.theta = math.radians(self.theta)
        norm = math.sqrt(l * l + m * m + n * n)
        self.l = l / norm
        self.m = m / norm
        self.n = n / norm

    def __str__(self):
        return "(%.6f,%.6f,%.6f,%.6f)" % (math.degrees(self.theta), self.l, self.m, self.n)

    def get_rotation_matrix(self):
        """
        return the rotation matrix
        """
        rot = numpy.zeros((3, 3))
        c = math.cos(self.theta)
        s = math.sin(self.theta)
        l = self.l
        m = self.m
        n = self.n
        rot[0, 0] = l * l * (1.0 - c) + c
        rot[1, 1] = m * m * (1.0 - c) + c
        rot[2, 2] = n * n * (1.0 - c) + c
        rot[0, 1] = m * l * (1.0 - c) - n * s
        rot[1, 0] = l * m * (1.0 - c) + n * s
        rot[0, 2] = n * l * (1.0 - c) + m * s
        rot[2, 0] = l * n * (1.0 - c) - m * s
        rot[1, 2] = n * m * (1.0 - c) - l * s
        rot[2, 1] = m * n * (1.0 - c) + l * s
        return rot

    def get_quaternion(self):
        """
        return the quaternion
        """
        c = math.cos(self.theta / 2.0)
        s = math.sin(self.theta / 2.0)
        n = numpy.array([self.l, self.m, self.n])
        return Quaternion(c, s * n)

    def get_angle(self, deg=False):
        """
        Param:deg:angle is in degrees instead of radians
        """
        angle = self.theta
        if deg:
            angle = numpy.degrees(angle)
        return angle

    def get_axis(self):
        return (self.l, self.m, self.n)


class Rotation_frame(object):
    """
    Get the rotation angle between two frames
    """

    def __init__(self, lab, mol):
        """
        Define a fix frame(lab frame) and a moving frame(mol frame)
        param:lab:numpy array of dimension 3x3 with each column defining a vector X, Y, Z
        param:mol:numpy array of dimension 3x3 with each column defining  a vector x, y, z
        """
        self.X = lab[:, 0]
        self.Y = lab[:, 1]
        self.Z = lab[:, 2]
        self.x = mol[:, 0]
        self.y = mol[:, 1]
        self.z = mol[:, 2]
        # check that each vector is normalized
        norm = numpy.linalg.norm(self.X)
        if norm < thr_norm:
            self.X = None
        else:
            self.X /= norm
        norm = numpy.linalg.norm(self.Y)
        if norm < thr_norm:
            self.Y = None
        else:
            self.Y /= norm
        norm = numpy.linalg.norm(self.Z)
        if norm < thr_norm:
            self.Z = None
        else:
            self.Z /= norm
        norm = numpy.linalg.norm(self.x)
        if norm < thr_norm:
            self.x = None
        else:
            self.x /= norm
        norm = numpy.linalg.norm(self.y)
        if norm < thr_norm:
            self.y = None
        else:
            self.y /= norm
        norm = numpy.linalg.norm(self.z)
        if norm < thr_norm:
            self.z = None
        else:
            self.z /= norm

        # check if both frame are right-handed frame
        if (self.X is not None) and (self.Y is not None) and (self.Z is not None):
            if (abs(numpy.cross(self.X, self.Y) - self.Z) > thr_dot).all():
                raise Exception("lab frame is not a right-handed system")
        if (self.x is not None) and (self.y is not None) and (self.z is not None):
            if (abs(numpy.cross(self.x, self.y) - self.z) > thr_dot).all():
                raise Exception("mol frame is not a right-handed system")

    def get_angles_zxz(self, deg=False):
        r"""
        Calculated the alpha, beta, gamma angles in the zx'z'' (or ZXZ) convention
        return:alpha, beta, gamma angles
        where beta is the angle between the z-axis and the Z-axis and range from 0 -> \pi
        where gamma is the angle between the x-axis and the line of nodes and range from -\pi -> \pi
        where alpha is the angle between the X-axis and the line of nodes and range from -\pi -> \pi
        The line of nodes is given by(Z cross z)
        Param:deg:alpha, beta, gamma are in degrees instead of radians
        """
        X = self.X
        Y = self.Y
        Z = self.Z
        x = self.x
        y = self.y
        z = self.z

        if z is None and Z is None:
            print("Both z and Z axis should be given")
            return None, None, None

        # Ambiguity if z and Z are colinear
        # beta = 0 or pi
        # alpha is the angle between X and x
        # gamma = 0
        Zz = numpy.dot(Z, z)
        if (1.0 - Zz) < thr_dot:
            # Z == z
            # We can choose N == x
            alpha = numpy.arctan2(numpy.dot(x, Y), numpy.dot(x, X))
            beta = 0.0
            gamma = 0.0
        elif (Zz + 1.0) < thr_dot:
            # Z == -z
            # We can choose N == x
            alpha = numpy.arctan2(numpy.dot(x, Y), numpy.dot(x, X))
            beta = numpy.pi
            gamma = 0.0
        else:
            #        N = numpy.cross(Z, z)
            #        N = N/numpy.linalg.norm(N)
            # beta = acos(z.Z)
            # beta range 0->\pi because N is defined as (Z cross z)
            beta = numpy.arccos(numpy.dot(Z, z))

            # alpha = atan2(z.X, -z.Y) = acos(X.N) * sign(Z.(X cross N)))
            if X is None and Y is None:
                alpha = None
            else:
                alpha = numpy.arctan2(numpy.dot(z, X), -numpy.dot(z, Y))
            # alpha2 = numpy.arccos(numpy.dot(X, N)) * numpy.sign(numpy.dot(Z, numpy.cross(X, N)))

            # gamma = atan2(x.Z, y.Z) = acos(x.N) * sign(z.(N cross x))
            if x is None and y is None:
                gamma = None
            else:
                gamma = numpy.arctan2(numpy.dot(x, Z), numpy.dot(y, Z))
            #  gamma2 = numpy.arccos(numpy.dot(x, N)) * numpy.sign(numpy.dot(z, numpy.cross(N, x)))

        if deg:
            if alpha is not None:
                alpha = numpy.degrees(alpha)
            if beta is not None:
                beta = numpy.degrees(beta)
            if gamma is not None:
                gamma = numpy.degrees(gamma)
        return alpha, beta, gamma

    def get_angles_zyz(self, deg=False):
        r"""
        Calculated the alpha, beta, gamma angles in the zy'z'' (or ZYZ) convention
        return:alpha, beta, gamma angles
        where beta is the angle between the z-axis and the Z-axis and range from 0 -> \pi
        where gamma is the angle between the y-axis and the line of nodes and range from -\pi -> \pi
        where alpha is the angle between the Y-axis and the line of nodes and range from -\pi -> \pi
        The line of nodes is given by(Z cross z)
        Param:deg:alpha, beta, gamma are in degrees instead of radians
        """
        X = self.X
        Y = self.Y
        Z = self.Z
        x = self.x
        y = self.y
        z = self.z

        if z is None and Z is None:
            print("Both z and Z axis should be given")
            return None, None, None

        # Ambiguity if z and Z are colinear
        # beta = 0 or pi
        # alpha is the angle between Y and y
        # gamma = 0
        Zz = numpy.dot(Z, z)
        if (1.0 - Zz) < thr_dot:
            # Z == z
            # We can choose N == y
            alpha = numpy.arctan2(-numpy.dot(y, X), numpy.dot(y, Y))
            beta = 0.0
            gamma = 0.0
        elif (Zz + 1.0) < thr_dot:
            # Z == -z
            # We can choose N == y
            alpha = numpy.arctan2(-numpy.dot(y, X), numpy.dot(y, Y))
            beta = numpy.pi
            gamma = 0.0
        else:
            #        N = numpy.cross(Z, z)
            #        N = N/numpy.linalg.norm(N)
            # beta = acos(z.Z)
            # beta range 0->\pi because N is defined as (z cross Z)
            beta = numpy.arccos(numpy.dot(Z, z))

            if X is None and Y is None:
                alpha = None
            else:
                alpha = numpy.arctan2(numpy.dot(z, Y), numpy.dot(z, X))

            if x is None and y is None:
                gamma = None
            else:
                gamma = numpy.arctan2(numpy.dot(y, Z), -numpy.dot(x, Z))

        if deg:
            if alpha is not None:
                alpha = numpy.degrees(alpha)
            if beta is not None:
                beta = numpy.degrees(beta)
            if gamma is not None:
                gamma = numpy.degrees(gamma)
        return alpha, beta, gamma

    def get_rotation_matrix(self):
        """
        return the rotation matrix from laboratory to molecular frame
        """
        alpha, beta, gamma = self.get_angles_zxz()
        ER = Euler_Rotation((alpha, beta, gamma), axes=(2, 0, 2), convention="intrinsic")
        rotmat = ER.get_rotation_matrix()
        return rotmat

    def get_angle_axis(self, deg=False):
        """
        Get the rotation of an angle theta around an axis of rotation
        to go from laboratory to molecular frame
        Angle:0->pi
        Param:deg:angle is in degrees instead of radians
        """
        return get_angle_axis(self.get_rotation_matrix(), deg=deg)


class Euler_Rotation(object):
    """
    Define a rotation by individual rotation around three axes
    In the "intrinsic" convention, the rotations are around the moving axes
    In the "extrinsic" convention, the rotations are around the fixed axes
    """

    def __init__(self, angles, axes=(2, 0, 2), convention="intrinsic", deg=False):
        """
        Define the convention, the axes and the angles
        Param:angles:tuple with the three angles in radian
        positive angles mean counter-clockwise rotation
        Param:axes:tuple with numbers for the three axes
               0 mean X, 1 -> Y, 2 -> Z
        Param:convention:either "intrinsic" or "extrinsic"
        By default, the rotation done is a zx'z'' or ZXZ
        Param:deg:angles are in degrees instead of radians
        """
        if len(angles) != 3:
            raise Exception("The angles should be a tuple of 3 values")
        self.angles = angles
        if deg:
            self.angles = tuple([numpy.radians(x) for x in self.angles])
        if len(axes) != 3:
            raise Exception("The axes should be a tuple of 3 values")
        self.axes = axes
        if convention == "intrinsic":
            self.intrinsic = True
        elif convention == "extrinsic":
            self.intrinsic = False
        else:
            raise Exception("The convention should either be 'intrinsic' or 'extrinsic'")

    def get_rotation_matrix(self):
        """
        return the rotation matrix
        For intrinsic:rot = axes[0](angles[0]) * axes[1](angles[1]) * axes[2](angles[2])
        For extrinsic:rot = axes[2](angles[2]) * axes[1](angles[1]) * axes[0](angles[0])
        Rotate the lab frame to go the mol frame = U_lab,mol
        """
        rotations = numpy.zeros((3, 3, 3))  # store the three rotation matrices
        for irot in range(3):
            rotations[irot] = get_arotation_matrix(self.axes[irot], self.angles[irot])

        # rotation is the global rotation i.e. the multiplication of the three
        # individual rotations
        if self.intrinsic:
            rotation = numpy.dot(numpy.dot(rotations[0], rotations[1]), rotations[2])
        else:
            rotation = numpy.dot(numpy.dot(rotations[2], rotations[1]), rotations[0])
        return rotation

    def get_angle_axis(self, deg=False):
        """
        Get the rotation of an angle theta around an axis of rotation
        Angle:0->pi
        Param:deg:angle is in degrees instead of radians
        """
        return get_angle_axis(self.get_rotation_matrix(), deg=deg)
