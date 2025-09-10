Installation
============

.. contents::
   :local:
   :depth: 2
   :backlinks: entry

.. important::

   GUANACO requires **Python 3.10 or higher**.  
   Please install a compatible version before proceeding.

Install from Source
-------------------

Step 1. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/Systems-Immunometabolism-Lab/guanaco-viz
   cd guanaco-viz

Step 2. Set up a Python environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can use **venv**, **Conda**, or **Pixi** to create an isolated environment.

**Using venv (built-in):**

.. code-block:: bash

   # Linux / macOS
   python3 -m venv myenv
   source myenv/bin/activate

   # Windows
   python -m venv myenv
   myenv\Scripts\activate

**Using Conda:**

.. code-block:: bash

   conda create -n guanaco-env
   conda activate guanaco-env

**Using Pixi:**

.. code-block:: bash

   pixi shell
   pixi add pip
   pip install .

Step 3. Install GUANACO
^^^^^^^^^^^^^^^^^^^^^^^

Inside the activated environment:

.. code-block:: bash

   pip install .

For development (editable mode):

.. code-block:: bash

   pip install -e .

Step 4. Run GUANACO
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   guanaco -c your_config.json

Additional options:

* ``guanaco -h`` → show help  
* ``guanaco -d DATA_DIR`` → specify data directory  
* ``guanaco -p PORT`` → change port (default: 4399)  
* ``guanaco --host HOST`` → change host (default: 0.0.0.0)  
* ``guanaco --debug`` → run in debug mode  
* ``guanaco --max-cells MAX_CELLS`` → adjust max cells (default: 8000)  
* ``guanaco --seed SEED`` → set random seed for subsampling  
* ``guanaco --backed-mode`` → recommended for large datasets, avoids memory issues  

.. note::

   It is advised to keep the data file and the config file in the same directory.

The application will be running at: http://0.0.0.0:4399/


Install with Docker
-------------------

Docker eliminates the need for Python setup. You only need your Json file and AnnData/MuData.

Step 1. Install Docker
^^^^^^^^^^^^^^^^^^^^^^

If not already installed, download Docker from https://www.docker.com/products/docker-desktop/.  

For Linux:

.. code-block:: bash

   sudo apt install -y docker.io

Verify installation:

.. code-block:: bash

   docker --version

Step 2. Run Docker Image
^^^^^^^^^^^^^^^^^^^^^^^^

GUANACO can also be run directly via docker container here: https://hub.docker.com/r/xiyuanzchloe/guanaco.
It is necessary to bind the directory containing the data and configuration file to a path inside the container.  

.. code-block:: bash

   docker run -p 8080:8080 \
   -v {directory of your data}:/app/data \
   -v {your json file}.json:/app/config.json \
   xiyuanzchloe/guanaco

The application will be available at: http://0.0.0.0:8080/