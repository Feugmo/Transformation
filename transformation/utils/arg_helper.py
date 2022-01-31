#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : Jan 05 6:33 p.m. 2022
@Author  : Conrard TETSASSI
@Email   : giresse.feugmo@gmail.com
@File    : arg_helper.py.py
@Project : Transformation
@Software: PyCharm
"""


import argparse
from transformation import __version__
import yaml


def parse_arguments():
    parser = argparse.ArgumentParser(description="Read a file, and perform a rotation and a scaling around a given "
                                                 "origin followed by a translation."
                                                 "The transformation is: "
                                                 "m = sR(m-ori) + ori + T")

    parser.add_argument('--version', action='version', version=__version__)

    parser.add_argument('--param_file', type=str, default="../params.yaml",
                        required=False, help="Path of config file")

    parser.add_argument("-i", "--inputfile", dest="inputfile",
                        help="Input  file, read  with ase",
                        action="store", type=str, default=None)

    parser.add_argument("--outputformat", dest="outpuformat",
                        help="Output format type . Default is xyz",
                        action="store", type=str, default="xyz")

    parser.add_argument("-o", "--outputfile", dest="outputfile",
                        help="Outputfile file. Default is the standard output",
                        action="store", type=str, default=None)

    parser.add_argument("-v", "--verbose", dest="verbose",
                        help="Print the position of the center of mass, center of charge, ...",
                        action="store_true", default=False)

    parser.add_argument("-s", "--scale", dest="scale",
                        help="Define a scaling. Default: 1.0",
                        action="store", type=str, default="1.0")

    parser.add_argument("-r", "--rotation", dest="rotation",
                        help="Define a counter-clockwise Rotation by an angle and an axis. Ex 30 deg along 'x' axis: "
                             "'30.0:1.0:0.0:0.0'",
                        action="store", type=str, default=None)

    parser.add_argument("--origin", dest="origin",
                        help="Define the origin from which the rotation and scaling are performed. Can be any vector "
                             "or one of the keyword ('cm'). Default: '0.0:0.0:0.0'",
                        action="store", type=str, default="0.0:0.0:0.0")

    parser.add_argument("-t", "--translation", dest="translation",
                        help="Define a Translation. Can be any vector or one of the keyword ('-cm', '-ori'). "
                             "Default: '0.0:0.0:0.0'",
                        action="store", type=str, default="0.0:0.0:0.0")

    parser.add_argument("--principal", dest="principal",
                        help="Rotate the molecule around the center of mass to align it along the principal axis of "
                             "inertia. Overwrite origin and rotation options.",
                        action="store_true", default=False)

    parser.add_argument('-a', "--axis", dest="axis",
                        help="Rotate the molecule around an axis defined by two atoms. The origin is set to the "
                             "middle position between the two atoms. Ex 30 degrees along the line between atom 1 and "
                             "5: '30.0:1:5'",
                        action="store", type=str, default=None)

    parser.add_argument("--bestvector", dest="bestvector",
                        help="Rotate the molecule to align its x, y, z, z axis with the vector that minimize the ("
                             "parallel or "
                             "perpendicular) error between the square distances of coordinates of a list of atoms. "
                             "Give the list of atoms (separated by comma) than'par' or 'per' (parallel or "
                             "perpendicular), than an atom number that define the direction of the vector.  Example: "
                             "1,2,3:par:10:z",
                        action="store", type=str, default=None)

    parser.add_argument("--bestcylinder", dest="bestcylinder",
                        help="Rotate the molecule to align its  x, y, z  axis with the vector that fit the best "
                             "the axial "
                             "direction of a cylinder. Give the list of atoms (separated by comma), than an atom "
                             "number that define the direction of the vector. Example: 1,2,3:10:z",
                        action="store", type=str, default=None)

    parser.add_argument("--helical", dest="helical",
                        help="Try to fit an helical axis on top of the bestcylinder approach. Give P or M to fit a "
                             "right- or left-handed helicoidal axis. Default: None",
                        action="store", type=str, default=None)

    args = parser.parse_args()

    return args


def get_params(param_file):
    param = yaml.safe_load(open(param_file, 'r'))

    return param
# parser = ar
