Running GUANACO on a Server
===========================

GUANACO can be deployed on cloud servers (e.g., AWS EC2) or on-premise servers, either using the system’s Python environment or Docker.

Before starting, connect to your server:

.. code-block:: bash

   ssh -i ~/filepath/my-aws-key.pem instance-name@<EC2_PUBLIC_IP>

.. note::
   - The instance name is usually the OS user (e.g. ``ubuntu`` or ``ec2-user``).  
   - ``ubuntu`` is preferable because it typically comes with updated Python.  
   - NGINX is optional but recommended for URL configuration and extra security.  


Deployment with System Python
-----------------------------

Step 1. Clone the repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   git clone https://github.com/Systems-Immunometabolism-Lab/guanaco-viz
   cd guanaco-viz

Upload your ``.json`` config and AnnData/MuData files into this folder, e.g. with:

.. code-block:: bash

   scp -i /path/to/key.pem /path/to/file ubuntu@<EC2_PUBLIC_IP>:/home/ubuntu/

Step 2. Create a virtual environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

From source:

.. code-block:: bash

   pip install .
   pip install -e .   # development mode


Step 4. Run GUANACO
^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   nohup guanaco -c your_config.json --backed-mode &

The app will be running at: http://<instance-public-ip>:4399/

Deployment with Docker
----------------------

Docker simplifies setup by bundling dependencies. Python installation is not required.

Step 1. Install Docker
^^^^^^^^^^^^^^^^^^^^^^

If not already installed:

.. code-block:: bash

   sudo apt update
   sudo apt install -y docker.io
   sudo systemctl enable docker
   sudo systemctl start docker

(Optional) Add Docker to user group:

.. code-block:: bash

   sudo usermod -aG docker $USER

Verify installation:

.. code-block:: bash

   docker run hello-world

Step 2. Run Docker image
^^^^^^^^^^^^^^^^^^^^^^^^

Upload your files to the server:

.. code-block:: bash

   scp -i /path/to/key.pem /path/to/guanaco.json ubuntu@<EC2_PUBLIC_IP>:/home/ubuntu/

Run the container (replace paths with your own):

.. code-block:: bash

   docker run -p 8080:8080 \
       -v {directory of your data}:/app/data \
       -v {your json file}.json:/app/config.json \
       systemsimmunometabolismlab/guanaco-viz

If Docker requires ``sudo``:

.. code-block:: bash

   sudo docker run -p 8080:8080 \
       -v {directory of your data}:/app/data \
       -v {your json file}.json:/app/config.json \
       systemsimmunometabolismlab/guanaco-viz

The app will be running at: http://<instance-public-ip>:8080/
