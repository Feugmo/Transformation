# from __future__ import absolute_import, division, print_function

import copy
import math
import sys

import numpy as np

from .Quaternions import Quaternion


def in_sphere(point, radius, center=None):
    if center is None:
        center = np.array([0, 0, 0])

    # Calculate the difference between the reference and measuring point
    diff = np.subtract(point, center)

    # Calculate square length of vector (distance between ref and point)^2
    dist = np.sum(np.power(diff, 2))
    # If dist is less than radius^2, return True, else return False
    return dist < radius ** 2


def in_box(points, box3d):
    """
    box3d  =  np array of the shape (8,3) with coordinates in the clockwise order. first the bottom plane is considered then the top one.
    points = array of points with shape (N, 3).

    Returns the indices of the points array which are outside the box3d
    adapted from https://stackoverflow.com/questions/21037241/how-to-determine-a-point-is-inside-or-outside-a-cube
    """
    b1, b2, b3, b4, t1, t2, t3, t4 = box3d

    dir1 = t1 - b1  # z
    size1 = np.linalg.norm(dir1)
    # dir1 = dir1 / size1  # the center is already [0,0,0]

    dir2 = b2 - b1  # y
    size2 = np.linalg.norm(dir2)
    # dir2 = dir2 / size2

    dir3 = b4 - b1  # x
    size3 = np.linalg.norm(dir3)
    # dir3 = dir3 / size3

    box3d_center = (b1 + t3) / 2.0

    dir_vec = points - box3d_center

    cell = np.array([size3, size2, size1])
    # cell = np.vstack([dir3, dir2, dir1])

    res1 = np.where((np.absolute(np.dot(dir_vec, dir1)) * 2) < size1)[0]  # x
    res2 = np.where((np.absolute(np.dot(dir_vec, dir2)) * 2) < size2)[0]  # y
    res3 = np.where((np.absolute(np.dot(dir_vec, dir3)) * 2) < size3)[0]  # z

    return intersection(intersection(res1, res2), res3), cell


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def in_box2(points, box):
    """
    box3d  =  np array of the shape (8,3) with coordinates in the clockwise order. first the bottom plane is considered then the top one.
    points = array of points with shape (N, 3).

    Returns the indices of the points array which are outside the box3d
    adapted from https://stackoverflow.com/questions/21037241/how-to-determine-a-point-is-inside-or-outside-a-cube
    """
    Xmax = np.max(box[:, 0])
    Xmin = np.min(box[:, 0])
    Ymax = np.max(box[:, 1])
    Ymin = np.min(box[:, 1])
    Zmax = np.max(box[:, 2])
    Zmin = np.min(box[:, 2])

    a = abs(Xmax) + abs(Xmin)
    b = abs(Ymax) + abs(Ymin)
    c = abs(Zmax) + abs(Zmin)
    cell = np.array([a, b, c])
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    res1 = np.where((x >= Xmin) & (x <= Xmax))[0]  # z
    res2 = np.where((y >= Ymin) & (y <= Ymax))[0]  # y
    res3 = np.where((z >= Zmin) & (z <= Zmax))[0]  # z

    return intersection(intersection(res1, res2), res3),cell


class MoleculeTools(object):
    def __init__(self, atms=None):
        self.atoms = atms
        self.natoms = len(atms)
        self.atmasses = atms.get_masses()
        self.atnums = atms.get_atomic_numbers()
        self.coordinates = atms.positions
        self.atsymbols = atms.get_chemical_symbols()

    def copy(self):
        mol = copy.deepcopy(self)
        return mol

    def get_centroid(self, weight=None, atomlist=None):
        """
        Param weight:list of weight of size natoms
        Param atomlist:list of atoms
        """
        if weight is None:
            weight = np.ones(self.natoms)
        if atomlist is None:
            atomlist = list(range(self.natoms))
        coords = self.coordinates
        if coords is None:
            return None
        #        centroid = np.zeros((3, ))
        #        for i in atomlist:
        #            centroid[:]=centroid[:]+coords[i, :]*weight[i]
        centroid = np.dot(weight[atomlist], coords[atomlist])
        return centroid / weight[atomlist].sum()

    def get_center_of_mass(self, atomlist=None):
        mass = self.atmasses
        if mass is None:
            return None
        return self.get_centroid(weight=mass, atomlist=atomlist)

    def get_center_of_charge(self, atomlist=None):
        charge = self.atnums
        if charge is None:
            return None
        return self.get_centroid(weight=charge, atomlist=atomlist)

    def set_to_cm(self):
        cm = self.get_center_of_mass()
        self.translate(-cm)

    def set_to_center_of_charge(self):
        cc = self.get_center_of_charge()
        self.translate(-cc)

    def get_Rsquare_matrix(self, weight=None, atomlist=None):
        r"""
        Rsquare matrix is a matrix whose elements are given by:
        Rquare_ {\alpha\beta}= \sum_i  R_ {i\alpha} R_ {i\beta} weight_i
        where R_ {i\alpha} = coord[i, \alpha] - centroid
        Param weight:list of weight of size natoms
        Param atomlist:list of atoms
        """
        if weight is None:
            weight = np.ones(self.natoms)
        if atomlist is None:
            atomlist = list(range(self.natoms))
        coords = self.coordinates
        if coords is None:
            return None
        centroid = self.get_centroid(weight, atomlist)
        coords = coords[atomlist].copy() - centroid
        Rsquare = np.dot(coords.transpose() * weight[atomlist], coords)
        return Rsquare

    def get_moment_of_inertia(self, atomlist=None):
        r"""
         moi_ {\alpha\beta} = sum_i m_i[r^2_i \delta_ {\alpha\beta} - r_ {i\alpha} r_ {i\beta}]
          = (\sum_i m_i r^2_i)\delta_ {\alpha\beta} - Rsquare_ {\alpha \beta}
        moi = \sum_i m_i[\vec {r} _i \cdot  \vec {r} _i \identity - \vec {r} _i \otimes \vec {r} _i]
        """
        if atomlist is None:
            atomlist = list(range(self.natoms))
        mass = self.atmasses
        coords = self.coordinates
        if (coords is None) or (mass is None):
            return None
        R2 = self.get_Rsquare_matrix(weight=mass, atomlist=atomlist)
        return np.trace(R2) * np.identity(3) - R2

    def rotate_to_principle_axis(self, atomlist=None):
        moi = self.get_moment_of_inertia(atomlist)
        _, eigenvecs = np.linalg.eigh(moi)
        # ensure that we have a right-handed system of axes
        eigenvecs[:, 2] = np.cross(eigenvecs[:, 0], eigenvecs[:, 1])
        self.rotate(eigenvecs.transpose())

    def get_moment_of_charge(self, atomlist=None):
        """
        similar to moment of inertia but with the charge instead of the mass
        """
        if atomlist is None:
            atomlist = list(range(self.natoms))
        coords = self.coordinates
        charge = self.atnums
        if (coords is None) or (charge is None):
            return None
        R2 = self.get_Rsquare_matrix(weight=charge, atomlist=atomlist)
        return np.trace(R2) * np.identity(3) - R2

    def rotate_to_moment_of_charge(self, atomlist=None):
        moi = self.get_moment_of_charge(atomlist)
        _, eigenvecs = np.linalg.eigh(moi)
        # ensure that we have a right-handed system of axes
        eigenvecs[:, 2] = np.cross(eigenvecs[:, 0], eigenvecs[:, 1])
        self.rotate(eigenvecs.transpose())

    def translate(self, vec, atomlist=None):
        if atomlist is None:
            self.coordinates = self.coordinates + np.array(vec)
        else:
            self.coordinates[atomlist] = self.coordinates[atomlist] + np.array(vec)

    def rotate(self, rotmat):
        self.coordinates = np.dot(self.coordinates, rotmat.transpose())

    def scale(self, s):
        # s is a scaling factor
        self.coordinates = self.coordinates * s

    def isLinear(self):
        """
        Guess if a molecule is linear
        """
        coords = self.coordinates
        natoms = coords.shape[0]
        if natoms > 1:
            # vector along bond 1-2
            v12 = coords[0] - coords[1]
            v12 = v12 / np.linalg.norm(v12)
            for coord in coords[2:]:
                v1i = coords[0] - coord
                v1i = v1i / np.linalg.norm(v1i)
                prod = np.dot(v12, v1i)  # between -1 -> 1
                val = 1.0 - abs(prod)  # between 0 -> 1 where 0 mean linear
                if val > 1.0e-6:
                    return False
            return True
        else:
            return None  # We have an atom

    def RMSD(self, mol2):
        """
        Calculate the Root Mean Square Deviation:
        Square root of the mean of the square of the distances between the atoms of two molecules
        """
        mol1 = self
        if mol1.natoms != mol2.natoms:
            print("Can't compare molecules with a different number of atoms")
            sys.exit(1)
        coord1 = mol1.coordinates
        coord2 = mol2.coordinates
        rmsd = ((coord1 - coord2) ** 2).sum()
        rmsd = math.sqrt(rmsd / (mol1.natoms))
        return rmsd

    def ChkAnm(self, molref):
        """
        Check no total angular momentum is present along the
        linear path connecting molref to self.
        0 = sum_i  atmass_i(molref x self)
        """
        thr = 1.0e-6
        mol = self
        mass = self.atmasses
        coords = mol.coordinates
        coordsref = molref.coordinates
        if mol.natoms != molref.natoms:
            print("Can't check angular momentum if the two molecules have a different number of atoms")
            sys.exit(1)
        L = np.zeros((3,))
        for iatom in range(mol.natoms):
            L += np.cross(coordsref[iatom], coords[iatom]) * mass[iatom]
        # check the three direction of vector L
        ok = True
        for l in L:
            if abs(l) > thr:
                ok = False
        if not ok:
            print("WARNING: Significant angular momentum along the path")
        print("A= %12.5e  %12.5e  %12.5e" % (L[0], L[1], L[2]))
        return ok

    def distance(self, iatom, jatom):
        """Calculate the distance between two atoms of the same molecule
        :param iatom, jatom:atom index(or list of atom indexes) zero-based
                return:one distance or a list of distances
        """
        if isinstance(iatom, int) or isinstance(iatom, np.int64):
            iatom = [iatom]
        if isinstance(jatom, int) or isinstance(jatom, np.int64):
            jatom = [jatom]
        if len(iatom) != len(jatom):
            raise Exception("distance: both lists should have the same length")
        coords = self.coordinates
        dist = []
        for i, j in zip(iatom, jatom):
            dist.append(np.linalg.norm(coords[i] - coords[j]))
        if len(dist) == 1:
            return dist[0]
        else:
            return dist

    @staticmethod
    def angle_vectors(vec1, vec2, deg=False):
        """Calculate the angle between two vectors.
        :param vec1, vec2:two ndarray vectors
        :return:the angle in radian
        """
        for vec in (vec1, vec2):
            if len(vec) != 3:
                raise Exception("angle_vectors:vectors are invalid")
        scalar_prod = np.dot(vec1, vec2)
        if abs(scalar_prod) < 1.0e-8:
            return np.pi

        val = scalar_prod / math.sqrt(np.dot(vec1, vec1) * np.dot(vec2, vec2))
        # sometimes acos behaves strange...
        if val > 1.0:
            val = 1.0
        elif val < -1.0:
            val = -1
        angle = np.arccos(val)
        if deg:
            angle = np.degrees(angle)
        return angle

    def angle(self, iatom, jatom, katom, deg=False):
        """Calculate the angle between three atoms of the same molecule
        :param iatom, jatom, katom:atom index(or list of atom indexes) zero-based
                return:one angle or a list of angles in radians
        """
        if isinstance(iatom, int) or isinstance(iatom, np.int64):
            iatom = [iatom]
        if isinstance(jatom, int) or isinstance(jatom, np.int64):
            jatom = [jatom]
        if isinstance(katom, int) or isinstance(katom, np.int64):
            katom = [katom]
        if (len(iatom) != len(jatom)) or (len(iatom) != len(katom)):
            raise Exception("distance: all lists should have the same length")
        coords = self.coordinates
        angle = []
        for i, j, k in zip(iatom, jatom, katom):
            angle.append(MoleculeTools.angle_vectors(coords[i] - coords[j], coords[k] - coords[j], deg=deg))
        if len(angle) == 1:
            return angle[0]
        else:
            return angle

    def dihedral(self, iatom, jatom, katom, latom, deg=False):
        """Calculate the dihedral angle between four atoms of the same molecule.
        :param iatom, jatom, katom, latom:atom index(or list of atom indexes) zero-based
                return:one torsion angle or a list of torsion angles in radians
        """
        if isinstance(iatom, int) or isinstance(iatom, np.int64):
            iatom = [iatom]
        if isinstance(jatom, int) or isinstance(jatom, np.int64):
            jatom = [jatom]
        if isinstance(katom, int) or isinstance(katom, np.int64):
            katom = [katom]
        if isinstance(latom, int) or isinstance(latom, np.int64):
            latom = [latom]
        if (len(iatom) != len(jatom)) or (len(iatom) != len(katom)) or (len(iatom) != len(latom)):
            raise Exception("distance: all lists should have the same length")
        coords = self.coordinates
        tangle = []
        for i, j, k, l in zip(iatom, jatom, katom, latom):
            v_ji = coords[i] - coords[j]
            v_jk = coords[k] - coords[j]
            v_kj = coords[j] - coords[k]
            v_kl = coords[l] - coords[k]
            # normal to the plane (i, j, k)
            norm1 = np.cross(v_ji, v_jk)
            # normal to the plane (j, k, l)
            norm2 = np.cross(v_kj, v_kl)
            # scalar triple product which defines the sign of the dihedral
            # angle
            if np.dot(v_jk, np.cross(norm1, norm2)) < 0.0:
                sign = -1.0
            else:
                sign = +1.0
            tangle.append(sign * MoleculeTools.angle_vectors(norm1, norm2, deg=deg))
        if len(tangle) == 1:
            return tangle[0]
        else:
            return tangle

    def distance_vect(self, i, j):
        """
        Return the distance vector between atom i and atom j
        """
        coords = self.coordinates

        return coords[j] - coords[i]

    def distances(self, mol2):
        """
        Return a list of vectors for the distances between each corresponding atoms
        """
        mol1 = self
        if mol1.natoms != mol2.natoms:
            print("Can't work on molecules with a different number of atoms")
            sys.exit(1)
        alist = []
        coord1 = mol1.coordinates
        coord2 = mol2.coordinates
        for iatom in range(mol1.natoms):
            alist.append(coord1[iatom] - coord2[iatom])
        return alist

    def match(self, mol2, improper_rotation=True, scale=None):
        """
        Match the coordinates between two molecules
        The output  rotation, scale, translation give the following transformation:
        mol1' = sR*mol1 + trans = mol2
        param improper_rotation:if yes, allow for improper rotation if they give a better match
        param scale:use this scaling if given, otherwise, get the best scaling
        (1) HORN, B. J Opt Soc Am A 1987, 4, 629-642.
        (2) Coutsias, E.; Seok, C.; Dill, K. J. Comput. Chem. 2004, 25, 1849-1857.
        """
        mol1 = self
        cm1 = mol1.get_center_of_mass()
        cm2 = mol2.get_center_of_mass()
        coord1 = mol1.coordinates
        coord2 = mol2.coordinates
        mass = mol1.atmasses
        if mol1.natoms != mol2.natoms:
            print("Can't match molecules with a different number of atoms")
            sys.exit(1)
        natoms = mol1.natoms
        for iatom in range(natoms):
            coord1[iatom] = coord1[iatom] - cm1
            coord2[iatom] = coord2[iatom] - cm2
        M = np.zeros((3, 3))
        norm1 = 0.0
        norm2 = 0.0
        for iatom in range(natoms):
            M[0, 0] += (coord1[iatom, 0] * coord2[iatom, 0]) * mass[iatom]
            M[1, 1] += (coord1[iatom, 1] * coord2[iatom, 1]) * mass[iatom]
            M[2, 2] += (coord1[iatom, 2] * coord2[iatom, 2]) * mass[iatom]
            M[0, 1] += (coord1[iatom, 0] * coord2[iatom, 1]) * mass[iatom]
            M[1, 0] += (coord1[iatom, 1] * coord2[iatom, 0]) * mass[iatom]
            M[0, 2] += (coord1[iatom, 0] * coord2[iatom, 2]) * mass[iatom]
            M[2, 0] += (coord1[iatom, 2] * coord2[iatom, 0]) * mass[iatom]
            M[1, 2] += (coord1[iatom, 1] * coord2[iatom, 2]) * mass[iatom]
            M[2, 1] += (coord1[iatom, 2] * coord2[iatom, 1]) * mass[iatom]
            norm1 += mass[iatom] * (np.linalg.norm(coord1[iatom])) ** 2
            norm2 += mass[iatom] * (np.linalg.norm(coord2[iatom])) ** 2
        N = np.zeros((4, 4))
        N[0, 0] = M[0, 0] + M[1, 1] + M[2, 2]
        N[1, 1] = M[0, 0] - M[1, 1] - M[2, 2]
        N[2, 2] = -M[0, 0] + M[1, 1] - M[2, 2]
        N[3, 3] = -M[0, 0] - M[1, 1] + M[2, 2]
        N[0, 1] = M[1, 2] - M[2, 1]
        N[1, 0] = M[1, 2] - M[2, 1]
        N[0, 2] = M[2, 0] - M[0, 2]
        N[2, 0] = M[2, 0] - M[0, 2]
        N[0, 3] = M[0, 1] - M[1, 0]
        N[3, 0] = M[0, 1] - M[1, 0]
        N[1, 2] = M[0, 1] + M[1, 0]
        N[2, 1] = M[0, 1] + M[1, 0]
        N[1, 3] = M[2, 0] + M[0, 2]
        N[3, 1] = M[2, 0] + M[0, 2]
        N[2, 3] = M[1, 2] + M[2, 1]
        N[3, 2] = M[1, 2] + M[2, 1]
        eigenvalue, eigenvector = np.linalg.eigh(N)
        #        print "eval", eigenvalue
        #        print "evec", eigenvector
        if scale is None:
            scale = math.sqrt(norm2 / norm1)
        if improper_rotation and (eigenvalue[3] < abs(eigenvalue[0])):
            # Rotoinversion gives the best fit
            rotation = Quaternion(eigenvector[0, 0], eigenvector[1:4, 0])
            scale = -abs(scale)
        else:
            # Rotation gives the best fit
            rotation = Quaternion(eigenvector[0, 3], eigenvector[1:4, 3])
            scale = abs(scale)
        mat = rotation.get_rotation_matrix()
        translation = cm2 - scale * np.dot(mat, cm1)
        return rotation, scale, translation

    def match_kneller(self, molref, scale=None):
        """
                Match the coordinates between two molecules
        (1) Kneller, G. Mol. Simulation 1991, 7, 113-119.
                The output  rotation, scale, translation give the following transformation:
                mol = sR*molref + trans
        """
        mol1 = self
        cm1 = mol1.get_center_of_mass()
        cm2 = molref.get_center_of_mass()
        coord1 = mol1.coordinates
        coord2 = molref.coordinates
        mass = molref.atmasses
        if mol1.natoms != molref.natoms:
            print("Can't match molecules with a different number of atoms")
            sys.exit(1)
        natoms = mol1.natoms
        for iatom in range(natoms):
            coord1[iatom] = coord1[iatom] - cm1
            coord2[iatom] = coord2[iatom] - cm2
        M = np.zeros((3, 3))
        norm1 = 0.0
        norm2 = 0.0
        for iatom in range(natoms):
            M[0, 0] += (coord1[iatom, 0] * coord2[iatom, 0]) * mass[iatom]
            M[1, 1] += (coord1[iatom, 1] * coord2[iatom, 1]) * mass[iatom]
            M[2, 2] += (coord1[iatom, 2] * coord2[iatom, 2]) * mass[iatom]
            M[0, 1] += (coord1[iatom, 0] * coord2[iatom, 1]) * mass[iatom]
            M[1, 0] += (coord1[iatom, 1] * coord2[iatom, 0]) * mass[iatom]
            M[0, 2] += (coord1[iatom, 0] * coord2[iatom, 2]) * mass[iatom]
            M[2, 0] += (coord1[iatom, 2] * coord2[iatom, 0]) * mass[iatom]
            M[1, 2] += (coord1[iatom, 1] * coord2[iatom, 2]) * mass[iatom]
            M[2, 1] += (coord1[iatom, 2] * coord2[iatom, 1]) * mass[iatom]
            norm1 += mass[iatom] * (np.linalg.norm(coord1[iatom])) ** 2
            norm2 += mass[iatom] * (np.linalg.norm(coord2[iatom])) ** 2
        N = np.zeros((4, 4))
        K = norm1 + norm2
        N[0, 0] = K - 2.0 * M[0, 0] - 2.0 * M[1, 1] - 2.0 * M[2, 2]
        N[1, 1] = K - 2.0 * M[0, 0] + 2.0 * M[1, 1] + 2.0 * M[2, 2]
        N[2, 2] = K + 2.0 * M[0, 0] - 2.0 * M[1, 1] + 2.0 * M[2, 2]
        N[3, 3] = K + 2.0 * M[0, 0] + 2.0 * M[1, 1] - 2.0 * M[2, 2]
        N[0, 1] = 2.0 * (M[1, 2] - M[2, 1])
        N[1, 0] = 2.0 * (M[1, 2] - M[2, 1])
        N[0, 2] = 2.0 * (M[2, 0] - M[0, 2])
        N[2, 0] = 2.0 * (M[2, 0] - M[0, 2])
        N[0, 3] = 2.0 * (M[0, 1] - M[1, 0])
        N[3, 0] = 2.0 * (M[0, 1] - M[1, 0])
        N[1, 2] = -2.0 * (M[0, 1] + M[1, 0])
        N[2, 1] = -2.0 * (M[0, 1] + M[1, 0])
        N[1, 3] = -2.0 * (M[2, 0] + M[0, 2])
        N[3, 1] = -2.0 * (M[2, 0] + M[0, 2])
        N[2, 3] = -2.0 * (M[1, 2] + M[2, 1])
        N[3, 2] = -2.0 * (M[1, 2] + M[2, 1])
        #        print "N", N
        _, eigenvector = np.linalg.eigh(N)
        #        print "eval", eigenvalue
        #        print "evec", eigenvector
        # the rotation is given from molref->mol
        rotation = Quaternion(eigenvector[0, 0], eigenvector[1:4, 0])
        if scale is None:
            scale = math.sqrt(norm1 / norm2)
        translation = cm1 - scale * np.dot(rotation.get_rotation_matrix(), cm2)
        return rotation, scale, translation

    def match_reflexion(self, mol2):
        """
        Match the coordinates between two molecules using a reflexion
        instead of rotation
        """
        mol1 = self
        cm1 = mol1.get_center_of_mass()
        cm2 = mol2.get_center_of_mass()
        coord1 = mol1.coordinates
        coord2 = mol2.coordinates
        if mol1.natoms != mol2.natoms:
            print("Can't match molecules with a different number of atoms")
            sys.exit(1)
        natoms = mol1.natoms
        for iatom in range(natoms):
            coord1[iatom] = coord1[iatom] - cm1
            coord2[iatom] = coord2[iatom] - cm2
        M = np.zeros((3, 3))
        norm1 = 0.0
        norm2 = 0.0
        for iatom in range(natoms):
            M[0, 0] += coord1[iatom, 0] * coord2[iatom, 0]
            M[1, 1] += coord1[iatom, 1] * coord2[iatom, 1]
            M[2, 2] += coord1[iatom, 2] * coord2[iatom, 2]
            M[0, 1] += coord1[iatom, 0] * coord2[iatom, 1]
            M[1, 0] += coord1[iatom, 1] * coord2[iatom, 0]
            M[0, 2] += coord1[iatom, 0] * coord2[iatom, 2]
            M[2, 0] += coord1[iatom, 2] * coord2[iatom, 0]
            M[1, 2] += coord1[iatom, 1] * coord2[iatom, 2]
            M[2, 1] += coord1[iatom, 2] * coord2[iatom, 1]
            norm1 += (np.linalg.norm(coord1[iatom])) ** 2
            norm2 += (np.linalg.norm(coord2[iatom])) ** 2
        N = np.zeros((4, 4))
        N[0, 0] = M[0, 0] + M[1, 1] + M[2, 2]
        N[1, 1] = -M[0, 0] + M[1, 1] + M[2, 2]
        N[2, 2] = M[0, 0] - M[1, 1] + M[2, 2]
        N[3, 3] = M[0, 0] + M[1, 1] - M[2, 2]
        N[0, 1] = -M[1, 2] + M[2, 1]
        N[1, 0] = M[1, 2] - M[2, 1]
        N[0, 2] = -M[2, 0] + M[0, 2]
        N[2, 0] = M[2, 0] - M[0, 2]
        N[0, 3] = -M[0, 1] + M[1, 0]
        N[3, 0] = M[0, 1] - M[1, 0]
        N[1, 2] = -(M[0, 1] + M[1, 0])
        N[2, 1] = -(M[0, 1] + M[1, 0])
        N[1, 3] = -(M[2, 0] + M[0, 2])
        N[3, 1] = -(M[2, 0] + M[0, 2])
        N[2, 3] = -(M[1, 2] + M[2, 1])
        N[3, 2] = -(M[1, 2] + M[2, 1])
        _, eigenvector = np.linalg.eigh(N)
        #        print "eval", eigenvalue
        #        print "evec", eigenvector
        reflexion = Quaternion(eigenvector[0, 3], eigenvector[1:4, 3])
        scale = math.sqrt(norm2 / norm1)
        translation = cm2 - scale * np.dot(reflexion.get_reflection_matrix(), cm1)
        return reflexion, scale, translation

    def nbatomtype(self, atnum):
        """
        Get the number of atoms of atomic number atnum in the molecule
        """
        atnums = self.atnums
        return np.where(atnums == atnum)[0].shape[0]

    def getatomtype(self, atnum):
        """
        Get the atoms of atomic number atnum in the molecule
        """
        atnums = self.atnums
        alist = np.where(atnums == atnum)[0]
        return self.get_fragment(alist)

    def FitCylinder(self, atomlist, guessVector=None):
        """"""
        from scipy.optimize import basinhopping

        def f(x):
            """
            x[0] = theta angle
            x[1] = phi angle
            """
            theta = x[0]
            phi = x[1]
            currentVector = np.array([math.cos(theta) * math.sin(phi), math.sin(theta) * math.sin(phi), math.cos(phi)])
            error, _, _ = self.getGErrorForCylinder(currentVector, atomlist=atomlist)
            return error

        if atomlist is None:
            atomlist = list(range(self.natoms))
        if self.coordinates is None:
            raise Exception("FitCylinder: no coordinate")
        #            return None
        centroid = self.get_centroid(atomlist=atomlist)

        guess = [0, 0]
        if guessVector is not None:
            guess = [math.acos(guessVector[2]), math.atan2(guessVector[1], guessVector[0])]
        #        print "guess", guess
        #        res = minimize(f, guess, method="Powell")
        res = basinhopping(f, guess, minimizer_kwargs={"method": "Powell"}, niter=100)
        #        if not res.success:
        #            raise Exception("Error in the minimization procedure of the cylinder")
        theta = res.x[0]
        phi = res.x[1]
        #        print "Result of optimization"
        #        print res
        vector = np.array([math.cos(theta) * math.sin(phi), math.sin(theta) * math.sin(phi), math.cos(phi)])
        minError, C, r2 = self.getGErrorForCylinder(vector, atomlist=atomlist)
        # transform minError into RMSD of r2
        minError = math.sqrt(minError / len(atomlist))
        return minError, vector, C + centroid, r2

    def FitHelicalAxis(self, atomlist, rotmatrix, center, radius, rightHelix):
        """
        Param rotmatrix:Rotation matrix(new, old) which rotate the system in order that Z is the helical Axis
        Param center:Center of the Helical Axis
        Param radius:radius of the Helical Axis
        Param rightHelix:right(P) Helix:True or left(M) Helix:False
        Return:return the error, the pitch and the phase of the Helical Axis
        """
        # coord wrt center of the cylinder
        coords = self.coordinates
        if coords is None:
            raise Exception("FitHelicalAxis: no coordinate")
        #            return None
        coords = coords[atomlist].copy() - center
        coords = np.dot(coords, rotmatrix.transpose())
        #        print coords
        sign = 1.0
        if not rightHelix:
            sign = -1.0
        from scipy.optimize import basinhopping

        def f(x):
            """
            x[0] = 1/lambda
            x[1] = 0 ->1
            """
            errorx = coords[:, 0] - radius * np.cos(2.0 * math.pi * (x[0] * coords[:, 2] + x[1]))
            errory = coords[:, 1] - sign * radius * np.sin(2.0 * math.pi * (x[0] * coords[:, 2] + x[1]))
            error = (errorx * errorx).sum() + (errory * errory).sum()
            #            print x, error
            return error

        class MyBounds(object):
            def __init__(self, xmax=[1.0, 1.0], xmin=[0.01, 0.0]):
                self.xmax = np.array(xmax)
                self.xmin = np.array(xmin)

            def __call__(self, **kwargs):
                x = kwargs["x_new"]
                tmax = bool(np.all(x <= self.xmax))
                tmin = bool(np.all(x >= self.xmin))
                return tmax and tmin

        mybounds = MyBounds()

        #        res = minimize(f, [2.0, 0.0], method="Powell")
        res = basinhopping(
            f,
            [1.0, 0.0],
            minimizer_kwargs={"method": "SLSQP", "bounds": ((0.01, 1.0), (0.0, 1.0))},
            niter=100,
            accept_test=mybounds,
        )
        #        if not res.success:
        #            raise Exception("Error in the minimization procedure of the Helical Axis")
        #        print "Result of optimization"
        #        print res
        # transform minError into RMSD
        minError = math.sqrt(res.fun / len(atomlist))
        pitch = 1.0 / res.x[0]
        phase = (res.x[1] % 1) * 2.0 * math.pi  # 0 -> 2 pi
        return minError, pitch, phase

    def getGErrorForCylinder(self, vector, atomlist=None):
        """
        Calculate the G function
        return the error, the Center and the radius square
        """
        if atomlist is None:
            atomlist = list(range(self.natoms))
        coords = self.coordinates
        if coords is None:
            raise Exception("getGerrorForCylinder: no coordinate")
        #            return None
        coords = coords[atomlist].copy() - self.get_centroid(atomlist=atomlist)
        natoms = len(atomlist)
        # P matrix = I - W W^T
        P = np.identity(3) - np.outer(vector, vector)
        # Define S, the skew symmetry matrix
        S = np.array([[0, -vector[2], vector[1]], [vector[2], 0, -vector[0]], [-vector[1], vector[0], 0]])
        # Define A = \sum_i Y_i Y_i^T where Y_i = P X_i the projection of X_i on the plane perpendicular to vector
        # Define B = (Y_i. Y_i) Y_i
        A = np.zeros((3, 3))
        B = np.zeros((3))
        Y = np.dot(coords, P)
        sqrLength = np.diag(np.dot(Y, Y.transpose()))
        meansqrLength = sqrLength.sum() / natoms
        for i in range(natoms):
            A += np.outer(Y[i], Y[i])
            B += sqrLength[i] * Y[i]
        A /= natoms
        B /= natoms
        # Define Ahat = -S A S
        Ahat = -np.dot(np.dot(S, A), S)
        # Calculate the Center on the plane PC
        PC = np.dot(Ahat, B) / np.trace(np.dot(Ahat, A))
        # calculate the radius square
        # calculate the error
        error = 0.0
        r2 = 0.0
        for i in range(natoms):
            term = sqrLength[i] - meansqrLength - 2.0 * np.dot(Y[i], PC)
            error += term * term
            diff = PC - Y[i]
            r2 += np.dot(diff, diff)
        error /= natoms
        r2 /= natoms
        return error, PC, r2
