# Getting started on NERSC

Below is a tutorial to set up a Python environment for the hackathon and use it with NERSC's Jupyter server.

## Build a Python environment

1. Log in to NERSC. Load the NERSC Python module:

`module load python`

2. Download the hackathon python environment package list from [this Github link](https://github.com/digital-earths-global-hackathon/tools/blob/main/python_envs/environment.yaml). Or run this command from your terminal:

`wget https://raw.githubusercontent.com/digital-earths-global-hackathon/tools/refs/heads/main/python_envs/environment.yaml`

You should see a file called `environment.yaml` in your current directory.

3. Set up a Python environment:

`mamba env create -f environment.yaml`

4. Activate the Python environment:

`conda activate easy`

Optional: run this command to test your environment:

`python -c "import easygems"`

If you don't see any error, you have all the packages needed for the hackathon.

## Set up Jupyter kernel

1. Install Jupyter kernel for your Python environment

`python -m ipykernel install --user --name=easy --display-name hackathon`

You should see a message like this: 

"Installed kernelspec easy in /global/../username/.local/share/jupyter/kernels/easy"

2. From your browser, connect to [NERSC Jupyter Hub](https://jupyter.nersc.gov/hub/home). Click The banner "National Energy Research Scientific Computing Center", then log in with your NERSC credential. 

3. You should see several buttons with different Node options. For typical use, click the `start` button under **Login Node**. The Jupyter server will start, and you will enter the Jupyter Hub interface.

## Run sample notebook

1. Download the sample notebook from [this Github link](https://github.com/digital-earths-global-hackathon/hk25-teams/blob/main/hk25-tutorials/simple_plot.ipynb). Or run this command from your terminal:

`wget https://raw.githubusercontent.com/digital-earths-global-hackathon/hk25-teams/refs/heads/main/hk25-tutorials/simple_plot.ipynb`

You should see a file called `simple_plot.ipynb` in your Jupyter Hub (left column under a list of files in your home directory). Double click it to open the Notebook.

![](https://raw.githubusercontent.com/digital-earths-global-hackathon/hk25-teams/main/hk25-USWest/images/getstart_img1.gif)

2. Click the upper right corner with the text "NERSC Python ◯", a pop-up menu will open, select the kernel "hackathon".

![](https://raw.githubusercontent.com/digital-earths-global-hackathon/hk25-teams/main/hk25-USWest/images/getstart_img2.gif)

3. Run the sample notebook.

More detailed documentation of Jupyter on NERSC is available [here](https://docs.nersc.gov/services/jupyter/).