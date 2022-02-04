#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : Jan 05 5:45 p.m. 2022
@Author  : Conrard TETSASSI
@Email   : giresse.feugmo@gmail.com
@File    : align.py.py
@Project : Transformation
@Software: PyCharm
"""
import math
import os
from pathlib import Path
import numpy as np
from ase.io import read, write
from transformation.tools import Molecule, Rotations
from transformation.utils.arg_helper import get_params, parse_arguments

aseformat = {"POSCAR": "vasp", "xyz": "xyz"}


def transform(param_file=None, replicate=None, writeoutput=True):
    args = parse_arguments()

    if param_file is not None:
        params = get_params(param_file)
        path = Path(param_file).parent
    else:
        raise Exception("Please Provide a parameter file")
    # else:
    #     params = get_params(args.param_file)
    #     path = Path(param_file).parent
    for key, value in vars(args).items():
        if key not in params.keys():
            params[key] = value
    # Treat inputfile paramsion
    if params["inputfile"] is None:
        raise Exception("No inputfile given")

    # open the file if it exists and read the geom
    if os.path.isfile(params["inputfile"]):
        atm = read(params["inputfile"])
    elif os.path.isfile(os.path.join(path, params["inputfile"])):
        params["inputfile"] = os.path.join(path, params["inputfile"])
        atm = read(params["inputfile"])
    else:
        raise Exception(f" {params['inputfile']} does not exist")

    if replicate:
        try:
            atm = atm * replicate
        except Exception as er:
            print(f"{er}")
    mol = Molecule.MoleculeTools(atm)
    along = None
    # # write('atm.xyz', atm, format='vasp')

    # write('atm_input.xyz', atm)

    # Treat other args
    try:
        params["scale"] = float(params["scale"])
    except ValueError:
        raise Exception("The scale factor cannot be converted into float")

    if params["rotation"] is None:

        params["rotation"] = np.identity(3)
    else:
        try:
            q = np.array([float(x) for x in args.rotation.split(":")])
        except ValueError:
            raise Exception("The rotation values cannot be converted into floats")
        angle = math.radians(q[0])
        vector = q[1:4] / np.linalg.norm(q[1:4])
        R = Rotations.Rotation(angle, vector[0], vector[1], vector[2])
        params["rotation"] = R

    if args.origin == "cm":
        params["origin"] = atm.get_center_of_mass()
    else:
        try:
            params["origin"] = np.array([float(x) for x in args.origin.split(":")])
        except ValueError:
            raise Exception("The origin values cannot be converted into floats")

    # rotate around an axis
    if params["axis"] is not None:
        try:
            vals = params["axis"].split(":")
            angle = math.radians(float(vals[0]))
            at1, at2 = [int(x) - 1 for x in vals[1:3]]  # zero based number
        except ValueError:
            raise Exception("The axis option is not valid")
        if (at1 < 0) or (at1 >= len(atm)):
            raise Exception("The atom1 in axis option is out of range")
        if (at2 < 0) or (at2 >= len(atm)):
            raise Exception("The atom2 in axis option is out of range")
        if at1 == at2:
            raise Exception("atom1 and atom2 are the same in axis option")
        c1 = atm.positions[at1]
        c2 = atm.positions[at2]
        dir = c2 - c1
        params["origin"] = (c1 + c2) / 2.0
        vector = dir / np.linalg.norm(dir)
        R = Rotations.Rotation(angle, vector[0], vector[1], vector[2])
        params["rotation"] = R

    # rotate to align the molecule along best vector that minimize some error
    if params["bestvector"] is not None:

        try:
            listatom, error, atom, along = params["bestvector"].split(":")
            print(f"rotate to align the molecule along best vector [{along}]")
            atom = int(atom) - 1
            listatom = [int(x) - 1 for x in listatom.split(",")]  # zero-based
            if error == "par":  # Normal
                eigcolumns = (0, 2)
            elif error == "per":  # Along
                eigcolumns = (2, 0)  # z= (2, 0) , y = (0, 2) x = (2,1)
            else:
                raise ValueError
        except ValueError:
            raise Exception("The bestvector option is not correct.")
        if (len(listatom) > 1 and eigcolumns[0] == 2) or (len(listatom) > 2 and eigcolumns[0] == 0):
            R2 = mol.get_Rsquare_matrix(atomlist=listatom)
            v, vv = np.linalg.eigh(R2)
            params["origin"] = mol.get_centroid(atomlist=listatom)
            rotmatrix = np.zeros((3, 3))
            # Z = eigenvector with smallest (largest) eigenvalue
            veca = mol.coordinates[atom] - params["origin"]
            # veca = atm.positions[atom] - params['origin']

            rotmatrix[2] = np.sign(np.dot(veca, vv[:, eigcolumns[0]])) * vv[:, eigcolumns[0]]
            rotmatrix[0] = vv[:, eigcolumns[1]]
            # Y = Z cross X
            rotmatrix[1] = np.cross(rotmatrix[2], rotmatrix[0])
            # fix the sign of Y (and therefore the sign of X) by using veca
            # which is the vector from centroid of points to a given atom specify by the user
            # reverse both X and Y if needed
            thesign = np.sign(np.dot(veca, rotmatrix[1]))
            rotmatrix[1] *= thesign
            rotmatrix[0] *= thesign
            params["rotation"] = rotmatrix  # new, old
        else:
            raise Exception(
                "The listatom for bestvector shoud have at least two (three) elements for perpendicular (parallel) error"
            )

    # rotate to align the molecule along best cylinder that minimize some error
    if params["bestcylinder"] is not None:
        try:
            listatom, atom, along = params["bestcylinder"].split(":")
            print(f"rotate to align the molecule along best cylinder [{along}]")
            atom = int(atom) - 1
            listatom = [int(x) - 1 for x in listatom.split(",")]  # zero-based
        except ValueError:
            raise Exception("The bestcylinder option is not correct.")
        if len(listatom) > 3:
            v, vv = np.linalg.eigh(mol.get_Rsquare_matrix(atomlist=listatom))
            error, vector, params["origin"], r2 = mol.FitCylinder(atomlist=listatom, guessVector=vv[:, 2])
            rotmatrix = np.zeros((3, 3))
            # Z = vector
            veca = mol.coordinates[atom] - params["origin"]
            # veca = atm.positions[atom] - params['origin']
            rotmatrix[2] = np.sign(np.dot(veca, vector)) * vector
            print("Radius of the cylinder: {:.3f}".format(math.sqrt(r2)))
            print("Error on the cylinder: {:14.6e}".format(error))
            # X is the axis that minimize the parallel error of bestvector
            rotmatrix[0] = vv[:, 0]
            # make sure X is perpendicular to Z
            rotmatrix[0] = np.dot(np.identity(3) - np.outer(rotmatrix[2], rotmatrix[2]), rotmatrix[0])
            # renormalize X
            rotmatrix[0] = rotmatrix[0] / np.linalg.norm(rotmatrix[0])
            # Y = Z cross X
            rotmatrix[1] = np.cross(rotmatrix[2], rotmatrix[0])
            # fix the sign of Y (and therefore the sign of X) by using veca
            # which is the vector from centroid of points to a given atom specify by the user
            # reverse both X and Y if needed
            thesign = np.sign(np.dot(veca, rotmatrix[1]))
            rotmatrix[1] *= thesign
            rotmatrix[0] *= thesign
            params["rotation"] = rotmatrix  # new, old
            # helical option
            if params["helical"] is not None:
                if params["helical"].upper() == "P":
                    rightHelix = True
                elif params["helical"].upper() == "M":
                    rightHelix = False
                else:
                    raise Exception("The helical option require either 'P' or 'M' value.")
                error, pitch, phase = mol.FitHelicalAxis(
                    atomlist=listatom,
                    rotmatrix=rotmatrix,
                    center=params["origin"],
                    radius=math.sqrt(r2),
                    rightHelix=rightHelix,
                )
                print("Error on the Helical Axis: {:14.6e}".format(error))
                print("Helical Pitch of the Helical Axis (Ang): {:.3f}".format(pitch))
                print("Phase of the Helical Axis (deg): {:.3f}".format(phase * 180.0 / math.pi))
        else:
            raise Exception("The listatom for bestcylinder shoud have at least four elements")

    # Rotate to principal axis of inertia
    if params["principal"]:
        params["origin"] = atm.get_center_of_mass()
        moi = mol.get_moment_of_inertia()
        eigenvals, eigenvecs = np.linalg.eigh(moi)
        # ensure that we have a right-handed system of axes
        eigenvecs[:, 2] = np.cross(eigenvecs[:, 0], eigenvecs[:, 1])
        params["rotation"] = eigenvecs.transpose()

    # translation option
    if params["translation"] == "-cm":
        params["translation"] = -atm.get_center_of_mass()
    elif params["translation"] == "-ori":
        params["translation"] = -params["origin"]
    else:
        try:
            params["translation"] = np.array([float(x) for x in params["translation"].split(":")])
        except ValueError:
            raise Exception("The translation values cannot be converted into floats")

    ## APPLY THE TRANSFORMATION ###

    # Rotation and scaling
    print("The scaling factor is: {:.6f}".format(params["scale"]))
    print("The rotation operation is:")
    if isinstance(params["rotation"], Rotations.Rotation):
        params["rotation"].get_quaternion().print_rotQuaternion()
        params["rotation"] = params["rotation"].get_rotation_matrix()
    print(params["rotation"])
    print("The origin is:", params["origin"])

    mol.translate(-params["origin"])

    mol.rotate(params["rotation"] * params["scale"])
    # mol.translate(params['origin'])

    # translation
    print("The translation operation is:", params["translation"])
    mol.translate(params["translation"])
    # print(atm.positions)
    atm.set_positions(mol.coordinates)
    # print(atm.positions)

    if along == "x":

        atm.rotate(90, "y")
    elif along == "y":
        atm.rotate(90, "x")

    p = Path(params["inputfile"])
    parent = str(p.parent.absolute())
    stem = p.stem
    # suffix = p.suffix
    outfile = os.path.join(parent, stem + "_Transformed." + params["outputformat"])

    if writeoutput:
        write(outfile, atm, format=aseformat[params["outputformat"]])

    return atm

# if __name__ == "__main__":
