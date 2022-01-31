Transformation
==============

Packages to re-orient molecule

Support and Documentation
-------------------------
see docs for documentation, reporting bugs, and getting support.


Example
-------------------------
Read the file  transformation/utils/arg_helper.py to find all the options

run the file align.py as example

edit the file params.yaml to change the parameter



Installation
-------------------------
- Clone the repo

.. sourcecode:: bash

    $ git clone https://github.com/Feugmo/Transformation.git

    $ cd Transformation


- Install poetry (https://github.com/python-poetry/poetry)

    + osx / linux / bashonwindows install instructions
        .. sourcecode:: bash

            recommended
            $ curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py | python -

            other possibility
            $ pip install poetry

    + windows powershell install instructions
        .. sourcecode:: bash

            recommended
            $ (Invoke-WebRequest -Uri https://raw.githubusercontent.com/python-poetry/poetry/master/install-poetry.py -UseBasicParsing).Content | python -

    + other possibility
            $ pip install poetry

- Once Poetry is installed you can execute the following:

.. sourcecode:: bash

    $ poetry --version

    $ poetry self update



- Create a Virtual Environment

    By default, Poetry create virtual environment in $HOME/.poetry/env or  $HOME/.cache/pypoetry/virtualenvs for cahcing/sharing purpose
        - use poetry env use python_version to specify the Python version to use for the project.
        .. sourcecode:: bash

            $Transformation poetry config virtualenvs.in-project true

        -   To change or otherwise add a new configuration setting,
        .. sourcecode:: bash

            $Transformation  poetry config virtualenvs.path /path/to/cache/directory/virtualenvs


- install the packages
.. sourcecode:: bash

    $Transformation poetry install

    $Transformation  poetry check

    $Transformation  poetry build


+ Listing the current configuration

    .. sourcecode:: bash

        $Transformation  poetry config --list


    which will give you something similar to this

    .. sourcecode:: bash

        cache-dir = "/path/to/cache/directory"
        virtualenvs.create = true
        virtualenvs.in-project = null
        virtualenvs.path = "{cache-dir}/virtualenvs"  # /path/to/cache/directory/virtualenvs

+ Show Information of the Vitual Environment

    .. sourcecode:: bash

        poetry env info


        Virtualenv
        Python:         x.x.x
        Implementation: CPython
        Path:           "/path/to/cache/poetry virtual environment"
        Valid:          True

        System
        Platform: linux
        OS:       posix
        Python:   /path/to/python

+ Activate Virtual Environment

.. sourcecode:: bash

    $Transformation  poetry shell


