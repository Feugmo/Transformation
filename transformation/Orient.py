#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from optparse import OptionParser

import os.path
from VibTools import Parser
from VibTools import Rotations
import numpy
import math

usage = "Usage: %prog [options] \nRead a QM file, and perform a rotation and a scaling around a given origin followed by a translation.\nThe transformation is: m' = sR(m-ori) + ori + T"
parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputfile", dest="inputfile",
                  help="Input xyz file, or any QM file",
                  action="store", type="string", default="none")

parser.add_option("-o", "--outputfile", dest="outputfile",
                  help="Outputfile xyz file. Default is the standard output",
                  action="store", type="string", default="none")

parser.add_option("-v", "--verbose", dest="verbose",
                  help="Print the position of the center of mass, center of charge, ...",
                  action="store_true", default=False)

parser.add_option("-s", "--scale", dest="scale",
                  help="Define a scaling. Default: 1.0",
                  action="store", type="string", default="1.0")

parser.add_option("-r", "--rotation", dest="rotation",
                  help="Define a counter-clockwise Rotation by an angle and an axis. Ex 30 deg along 'x' axis: '30.0:1.0:0.0:0.0'",
                  action="store", type="string", default="none")

parser.add_option("", "--origin", dest="origin",
                  help="Define the origin from which the rotation and scaling are performed. Can be any vector or one of the keyword ('cm','cc'). Default: '0.0:0.0:0.0'",
                  action="store", type="string", default="0.0:0.0:0.0")

parser.add_option("-t", "--translation", dest="translation",
                  help="Define a Translation. Can be any vector or one of the keyword ('-cm','-cc', '-ori'). Default: '0.0:0.0:0.0'",
                  action="store", type="string", default="0.0:0.0:0.0")

parser.add_option("", "--principal", dest="principal",
                  help="Rotate the molecule around the center of mass to align it along the principal axis of inertia. Overwrite origin and rotation options.",
                  action="store_true", default=False)

parser.add_option("", "--moment_of_charge", dest="moc",
                  help="Rotate the molecule aounnd the center of charge to align it along the principal axis of moment of charge. Overwrite origin and rotation options",
                  action="store_true", default=False)

parser.add_option("-f", "--fragmentfile", dest="fragmentfile",
                  help="File containing atoms to orient while the other atoms are kept fixed",
                  action="store", type="string", default="none")

parser.add_option('-a', "--axis", dest="axis",
                  help="Rotate the molecule around an axis defined by two atoms. The origin is set to the middle position between the two atoms. Ex 30 degrees along the line between atom 1 and 5: '30.0:1:5'",
                  action="store", type="string", default="none")

parser.add_option("", "--bestvector", dest="bestvector",
                  help="Rotate the molecule to align its z axis with the vector that minimize the (parallel or perpendicular) error between the square distances of coordinates of a list of atoms. Give the list of atoms (separated by comma) than'par' or 'per' (parallel or perpendicular), than an atom number that define the direction of the vector.  Example: 1,2,3:par:10",
                  action="store", type="string", default="none")

parser.add_option("", "--bestcylinder", dest="bestcylinder",
                  help="Rotate the molecule to align its z axis with the vector that fit the best the axial direction of a cylinder. Give the list of atoms (separated by comma), than an atom number that define the direction of the vector. Example: 1,2,3:10",
                  action="store", type="string", default="none")

parser.add_option("", "--helical", dest="helical",
                  help="Try to fit an helical axis on top of the bestcylinder approach. Give P or M to fit a right- or left-handed helicoidal axis. Default: None",
                  action="store", type="string", default="none")

(options, args) = parser.parse_args()

# Treat inputfile option
if options.inputfile == "none":
    parser.error("No inputfile given")

# open the file if it exists and read the geom
if not os.path.isfile(options.inputfile):
    parser.error("The file ['%s'] does not exist" % options.inputfile)

file = Parser.parse_file(options.inputfile, XYZfile=True)
file.read()

# initialize molecule
molecule = file.mol

atomlist = None
# define the fragments
if options.fragmentfile != "none":
    if not os.path.isfile(options.fragmentfile):
        parser.error("The file ['%s'] does not exist" % options.fragmentfile)
    fragfile = open(options.fragmentfile, 'r')
    atomlist = []
    for line in fragfile:
        if len(line.strip()) != 0 and ("#" not in line):
            data = [int(d) - 1 for d in line.strip().split()]  # zero-based numbers
            atomlist.extend(data)
    fragfile.close()
    otheratomlist = list(set(range(file.mol.natoms)) - set(atomlist))
    molecule = file.mol.get_fragment(atomlist)

# Treat other options
opt = {}
try:
    opt['scale'] = float(options.scale)
except ValueError:
    parser.error("The scale factor cannot be converted into float")

if options.rotation == "none":
    opt['rotation'] = numpy.identity(3)
else:
    try:
        q = numpy.array([float(x) for x in options.rotation.split(":")])
    except ValueError:
        parser.error("The rotation values cannot be converted into floats")
    angle = math.radians(q[0])
    vector = q[1:4] / numpy.linalg.norm(q[1:4])
    R = Rotations.Rotation(angle, vector[0], vector[1], vector[2])
    opt['rotation'] = R

if options.origin == 'cm':
    opt['origin'] = molecule.get_center_of_mass()
elif options.origin == 'cc':
    opt['origin'] = molecule.get_center_of_charge()
else:
    try:
        opt['origin'] = numpy.array([float(x) for x in options.origin.split(":")])
    except ValueError:
        parser.error("The origin values cannot be converted into floats")

# rotate around an axis
if options.axis != "none":
    try:
        vals = options.axis.split(":")
        angle = math.radians(float(vals[0]))
        at1, at2 = [int(x) - 1 for x in vals[1:3]]  # zero based number
    except ValueError:
        parser.error("The axis option is not valid")
    if (at1 < 0) or (at1 >= file.mol.natoms):
        parser.error("The atom1 in axis option is out of range")
    if (at2 < 0) or (at2 >= file.mol.natoms):
        parser.error("The atom2 in axis option is out of range")
    if at1 == at2:
        parser.error("atom1 and atom2 are the same in axis option")
    c1 = file.mol.coordinates[at1]
    c2 = file.mol.coordinates[at2]
    dir = c2 - c1
    opt['origin'] = (c1 + c2) / 2.0
    vector = dir / numpy.linalg.norm(dir)
    R = Rotations.Rotation(angle, vector[0], vector[1], vector[2])
    opt['rotation'] = R

# rotate to align the molecule along best vector that minimize some error
if options.bestvector != "none":
    try:
        listatom, error, atom = options.bestvector.split(":")
        atom = int(atom) - 1
        listatom = [int(x) - 1 for x in listatom.split(",")]  # zero-based
        if error == "par":  # Normal
            eigcolumns = (0, 2)
        elif error == "per":  # Along
            eigcolumns = (2, 0)
        else:
            raise ValueError
    except ValueError:
        parser.error("The bestvector option is not correct.")
    if (len(listatom) > 1 and eigcolumns[0] == 2) or (len(listatom) > 2 and eigcolumns[0] == 0):
        R2 = file.mol.get_Rsquare_matrix(atomlist=listatom)
        v, vv = numpy.linalg.eigh(R2)
        opt['origin'] = file.mol.get_centroid(atomlist=listatom)
        rotmatrix = numpy.zeros((3, 3))
        # Z = eigenvector with smallest (largest) eigenvalue
        veca = file.mol.coordinates[atom] - opt['origin']
        rotmatrix[2] = numpy.sign(numpy.dot(veca, vv[:, eigcolumns[0]])) * vv[:, eigcolumns[0]]
        rotmatrix[0] = vv[:, eigcolumns[1]]
        # Y = Z cross X
        rotmatrix[1] = numpy.cross(rotmatrix[2], rotmatrix[0])
        # fix the sign of Y (and therefore the sign of X) by using veca
        # which is the vector from centroid of points to a given atom specify by the user
        # reverse both X and Y if needed
        thesign = numpy.sign(numpy.dot(veca, rotmatrix[1]))
        rotmatrix[1] *= thesign
        rotmatrix[0] *= thesign
        opt['rotation'] = rotmatrix  # new, old
    else:
        parser.error(
            "The listatom for bestvector shoud have at least two (three) elements for perpendicular (parallel) error")

# rotate to align the molecule along best cylinder that minimize some error
if options.bestcylinder != "none":
    try:
        listatom, atom = options.bestcylinder.split(":")
        atom = int(atom) - 1
        listatom = [int(x) - 1 for x in listatom.split(",")]  # zero-based
    except ValueError:
        parser.error("The bestcylinder option is not correct.")
    if len(listatom) > 3:
        v, vv = numpy.linalg.eigh(file.mol.get_Rsquare_matrix(atomlist=listatom))
        error, vector, opt["origin"], r2 = file.mol.FitCylinder(atomlist=listatom, guessVector=vv[:, 2])
        rotmatrix = numpy.zeros((3, 3))
        # Z = vector
        veca = file.mol.coordinates[atom] - opt['origin']
        rotmatrix[2] = numpy.sign(numpy.dot(veca, vector)) * vector
        print("Radius of the cylinder: %.3f" % (math.sqrt(r2)))
        print("Error on the cylinder: %14.6e" % (error))
        rotmatrix[0] = vv[:, 0]  # X is the axis that minimize the parallel error of bestvector
        # make sure X is perpendicular to Z
        rotmatrix[0] = numpy.dot(numpy.identity(3) - numpy.outer(rotmatrix[2], rotmatrix[2]), rotmatrix[0])
        # renormalize X
        rotmatrix[0] = rotmatrix[0] / numpy.linalg.norm(rotmatrix[0])
        # Y = Z cross X
        rotmatrix[1] = numpy.cross(rotmatrix[2], rotmatrix[0])
        # fix the sign of Y (and therefore the sign of X) by using veca
        # which is the vector from centroid of points to a given atom specify by the user
        # reverse both X and Y if needed
        thesign = numpy.sign(numpy.dot(veca, rotmatrix[1]))
        rotmatrix[1] *= thesign
        rotmatrix[0] *= thesign
        opt['rotation'] = rotmatrix  # new, old
        # helical option
        if options.helical != "none":
            if options.helical.upper() == "P":
                rightHelix = True
            elif options.helical.upper() == "M":
                rightHelix = False
            else:
                parser.error("The helical option require either 'P' or 'M' value.")
            error, pitch, phase = file.mol.FitHelicalAxis(atomlist=listatom, rotmatrix=rotmatrix, center=opt["origin"],
                                                          radius=math.sqrt(r2), rightHelix=rightHelix)
            print("Error on the Helical Axis: %14.6e" % (error))
            print("Helical Pitch of the Helical Axis (Ang): %.3f" % (pitch))
            print("Phase of the Helical Axis (deg): %.3f" % (phase * 180.0 / math.pi))
    else:
        parser.error("The listatom for bestcylinder shoud have at least four elements")

# Rotate to principal axis of inertia
if options.principal:
    opt['origin'] = molecule.get_center_of_mass()
    moi = molecule.get_moment_of_inertia()
    eigenvals, eigenvecs = numpy.linalg.eigh(moi)
    # ensure that we have a right-handed system of axes
    eigenvecs[:, 2] = numpy.cross(eigenvecs[:, 0], eigenvecs[:, 1])
    opt['rotation'] = eigenvecs.transpose()
elif options.moc:
    opt['origin'] = molecule.get_center_of_charge()
    moi = molecule.get_moment_of_charge()
    eigenvals, eigenvecs = numpy.linalg.eigh(moi)
    # ensure that we have a right-handed system of axes
    eigenvecs[:, 2] = numpy.cross(eigenvecs[:, 0], eigenvecs[:, 1])
    opt['rotation'] = eigenvecs.transpose()

# translation option
if options.translation == '-cm':
    opt['translation'] = -molecule.get_center_of_mass()
elif options.translation == '-cc':
    opt['translation'] = -molecule.get_center_of_charge()
elif options.translation == '-ori':
    opt['translation'] = -opt["origin"]
else:
    try:
        opt['translation'] = numpy.array([float(x) for x in options.translation.split(":")])
    except ValueError:
        parser.error("The translation values cannot be converted into floats")

## APPLY THE TRANSFORMATION ###

# Rotation and scaling
print("The scaling factor is: %.6f" % (opt['scale']))
print("The rotation operation is:")
if isinstance(opt['rotation'], Rotations.Rotation):
    opt['rotation'].get_quaternion().print_rotQuaternion()
    opt['rotation'] = opt['rotation'].get_rotation_matrix()
print(opt['rotation'])
print("The origin is:", opt['origin'])
molecule.translate(-opt['origin'])
molecule.rotate(opt['rotation'] * opt['scale'])
molecule.translate(opt['origin'])

# translation
print("The translation operation is:", opt['translation'])
molecule.translate(opt['translation'])

# Print info for the final geom
if options.verbose:
    print("Final geometry mol' = sR(mol-ori) + ori + T ")
    cm = molecule.get_center_of_mass()
    cc = molecule.get_center_of_charge()
    inertia = molecule.get_moment_of_inertia()
    print("Center of mass (Ang): %10.6f  %10.6f  %10.6f" % (cm[0], cm[1], cm[2]))
    print("Center of charge (Ang): %10.6f  %10.6f  %10.6f" % (cc[0], cc[1], cc[2]))
    print("Inertia tensor (Ang**2 amu):")
    print("%12.4f  %12.4f  %12.4f" % (inertia[0, 0], inertia[0, 1], inertia[0, 2]))
    print("%12.4f  %12.4f  %12.4f" % (inertia[1, 0], inertia[1, 1], inertia[1, 2]))
    print("%12.4f  %12.4f  %12.4f" % (inertia[2, 0], inertia[2, 1], inertia[2, 2]))
    moc = molecule.get_moment_of_charge()
    print("moment of charge tensor (Ang**2 q_e):")
    print("%12.4f  %12.4f  %12.4f" % (moc[0, 0], moc[0, 1], moc[0, 2]))
    print("%12.4f  %12.4f  %12.4f" % (moc[1, 0], moc[1, 1], moc[1, 2]))
    print("%12.4f  %12.4f  %12.4f" % (moc[2, 0], moc[2, 1], moc[2, 2]))

# Add the fixed coordinates to the new coordinates
if atomlist is not None:
    molecule = molecule + file.mol.get_fragment(otheratomlist)

atvalues = molecule.atvalues

# print or write the new coordinates
if options.outputfile == "none":
    print("Atomic number      Coordinates (Ang)")
    iatom = 0
    for z, c in zip(molecule.atnums, molecule.coordinates):
        print("%3i   %14.8f %14.8f %14.8f" % (z, c[0], c[1], c[2]), end=' ')
        if atvalues[iatom] is not None:
            print(" %14.8f" % (atvalues[iatom]), end=' ')
        print()
        iatom = iatom + 1
else:
    # write the molecule to outputfile
    molecule.write_to_xyz(options.outputfile)
