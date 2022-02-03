Transformation
==============

Packages to re-orient molecule

Support and Documentation
-------------------------
see docs for documentation, reporting bugs, and getting support.


Usage
-------------------------
+ Read the file  transformation/utils/arg_helper.py to find all the options



+ edit the file params.yaml to change the parameter
    + **inputfile:**  Input  file, read  with ase
    + **outpuformat :** Output format type . Default is xyz",
    + **outputfile:** Outputfile file. Default is the standard output",
    + **scale:** Define a scaling. Default: 1.0",
    + **rotation :**   Define a counter-clockwise Rotation by an angle and an axis. Ex 30 deg along 'x' axis:  '30.0:1.0:0.0:0.0'
    + **origin:**  Define the origin from which the rotation and scaling are performed. Can be any vector
                             or one of the keyword ('cm'). Default: '0.0:0.0:0.0',

    + **translation:**  Define a Translation. Can be any vector or one of the keyword ('-cm', '-ori'). "
                             "Default: '0.0:0.0:0.0'",
    + **principal:** Rotate the molecule around the center of mass to align it along the principal axis of "
                             "inertia. Overwrite origin and rotation options."

    + **axis:** "Rotate the molecule around an axis defined by two atoms. The origin is set to the "
                             "middle position between the two atoms. Ex 30 degrees along the line between atom 1 and "
                             "5: '30.0:1:5'",

    + **bestvector:** "Rotate the molecule to align its x, y, z, z axis with the vector that minimize the ("
                             "parallel or "
                             "perpendicular) error between the square distances of coordinates of a list of atoms. "
                             "Give the list of atoms (separated by comma) than'par' or 'per' (parallel or "
                             "perpendicular), than an atom number that define the direction of the vector.  Example: "
                             "1,2,3:par:10:z",
    + **bestcylinder:** "Rotate the molecule to align its  x, y, z  axis with the vector that fit the best "
                             "the axial "
                             "direction of a cylinder. Give the list of atoms (separated by comma), than an atom "
                             "number that define the direction of the vector. Example: 1,2,3:10:z",




+ run the file transformation/align.py



Installation
-------------------------
- Clone the repo

.. sourcecode:: bash

    $ git clone https://github.com/Feugmo/Transformation.git

    $ cd Transformation




