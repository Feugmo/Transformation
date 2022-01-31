import numpy


class Quaternion(object):
    def __init__(self, ps, pv):
        """
        Define a quaternion(Q=ps + vect(p))
        """
        self.s = ps
        if isinstance(pv, list):
            self.v = numpy.array(pv)
        elif isinstance(pv, numpy.ndarray):
            self.v = pv
        else:
            raise Exception("Error in the type of pv. Not a list or numpy.ndarray")

    def __str__(self):
        return "(%.6f,%.6f,%.6f,%.6f)" % (self.s, self.v[0], self.v[1], self.v[2])

    def __neg__(self):
        """
        Unary "-" operation
        """
        return Quaternion(-self.s, -self.v)

    def print_rotQuaternion(self):
        """
        Q = cos(theta/2) + sin(theta/2) * n
        where n is a unit vector along the rotation axis
        We want theta from -pi -> pi
        But theta, n gives the same quaternion as -theta, -n
        Therefore, theta:0 -> pi, theta/2:0 -> pi/2
        cos > 0 while sin > 0
        If cos < 0, Q' = -Q:correspond to have theta + 2pi
        We can allow theta to be -pi -> pi if we restrict n to be in 4 octant instead of the 8 octant.

        """
        Q = self
        if (Q.norm() - 1.0) > 1.0e-8:
            raise Exception("Q is not a unit quaternion")
        if Q.s < 0.0:
            Q = -Q
        cos = Q.s
        sin = numpy.linalg.norm(Q.v)
        if abs(sin) < 1.0e-8:
            print("Theta= %7.2f" % 0.0)
            return
        v = Q.v / sin
        octant = v // 1
        if sum(octant) < -1:
            v = -v
            sin = -sin
        thetahalf = numpy.arctan2(sin, cos)
        print("Theta= %7.2f, Unit vector = (%.6f,%.6f,%.6f)" % (numpy.degrees(thetahalf * 2.0), v[0], v[1], v[2]))

    def conjugate(self):
        return Quaternion(self.s, -self.v)

    def mult(self, q):
        p = self
        return Quaternion(p.s * q.s - numpy.dot(p.v, q.v), p.s * q.v + q.s * p.v + numpy.cross(p.v, q.v))

    def norm(self):
        return numpy.sqrt(self.s * self.s + numpy.dot(self.v, self.v))

    def normalized(self):
        norm = self.norm()
        self.s = self.s / norm
        self.v = self.v / norm

    def reciprocal(self):
        norm = self.norm()
        return Quaternion(self.s / (norm * norm), -self.v / (norm * norm))

    def Rotate(self, r):
        """
        Do the r'= Q r Q* transformation where Q is a unit quaternion and r
        and r' are vectors. Correspond to a counter-clockwise rotation
        """
        Q = self
        if (Q.norm() - 1.0) > 1.0e-8:
            raise Exception("Rotate: Q is not a unit quaternion")
        invQ = self.conjugate()
        r = Quaternion(0.0, r)
        dummy = r.mult(invQ)
        return Q.mult(dummy).v

    def Reflection(self, r):
        """
        Do the r'= Q r Q transformation where Q is a unit imaginary
        quaternion and r and r' are vectors. Correspond to a
        reflection through the plane perpendicular to the vector given
        by the imaginary part of the quaternion
        """
        Q = self
        if (Q.s > 1.0e-8) or ((Q.norm() - 1.0) > 1.0e-8):
            raise Exception("Reflection: Q is not a imaginary unit quaternion")
        r = Quaternion(0.0, r)
        dummy = r.mult(Q)
        return Q.mult(dummy).v

    def get_rotation_matrix(self):
        """
        if a quaternion is a unit quaternion, it can represent a
        counter-clockwise rotation.
        """
        if (self.norm() - 1.0) > 1.0e-8:
            raise Exception("get_rotation_matrix: Q is not a unit quaternion")
        rot = numpy.zeros((3, 3))
        q0 = self.s
        qx = self.v[0]
        qy = self.v[1]
        qz = self.v[2]
        rot[0, 0] = q0 * q0 + qx * qx - qy * qy - qz * qz
        rot[1, 1] = q0 * q0 - qx * qx + qy * qy - qz * qz
        rot[2, 2] = q0 * q0 - qx * qx - qy * qy + qz * qz
        rot[0, 1] = 2.0 * (qx * qy - q0 * qz)
        rot[1, 0] = 2.0 * (qy * qx + q0 * qz)
        rot[0, 2] = 2.0 * (qx * qz + q0 * qy)
        rot[2, 0] = 2.0 * (qz * qx - q0 * qy)
        rot[1, 2] = 2.0 * (qy * qz - q0 * qx)
        rot[2, 1] = 2.0 * (qz * qy + q0 * qx)
        return rot

    def get_reflection_matrix(self):
        """
        if a quaternion is a unit imaginary quaternion, it can represent a
        reflection.
        """
        if (self.norm() - 1.0) > 1.0e-8:
            raise Exception("get_reflection_matrix: Q is not a unit quaternion")
        ref = numpy.zeros((3, 3))
        q0 = self.s
        if q0 > 1.0e-3:
            raise Exception("Reflection matrix: Q is not a imaginary quaternion")
        qx = self.v[0]
        qy = self.v[1]
        qz = self.v[2]
        ref[0, 0] = -qx * qx + qy * qy + qz * qz
        ref[1, 1] = qx * qx - qy * qy + qz * qz
        ref[2, 2] = qx * qx + qy * qy - qz * qz
        ref[0, 1] = -2.0 * qx * qy
        ref[1, 0] = -2.0 * qy * qx
        ref[0, 2] = -2.0 * qx * qz
        ref[2, 0] = -2.0 * qz * qx
        ref[1, 2] = -2.0 * qy * qz
        ref[2, 1] = -2.0 * qz * qy
        return ref
