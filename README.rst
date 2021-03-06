Transformation
==============

Package to re-orient molecules or crystal structures



Usage
-------------------------

+ python transformation/align.py  --param_file  params.yaml

+
.. code-block:: python

    from transformation.align import transform
    from transformation.tools import Molecule

    # re-orient the molecule
    atms = transform(param_file=<<path to params file>>, writeoutput=True)

    # find the atoms inside a box : box is an array with 8 components
    inside, cell = Molecule.in_box2(atms.positions, box)


+ Edit the file params.yaml to change the parameter
    + **inputfile:**  Input  file (readable by ASE)
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








Installation
-------------------------
- Clone the repo

.. sourcecode:: bash

    $ git clone https://github.com/Feugmo/Transformation.git

    $ cd Transformation

    $Transformation poetry install

    $Transformation poetry check

    $Transformation poetry run pytest




Support and Documentation
-------------------------
see docs for documentation, reporting bugs, and getting support.


