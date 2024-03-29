# dtu-digital-twin-mm

This exercise evaluates the learning objectives for the Digital Twin Micro Module. After completing the module, you should be able to
- interact with a digital twin in order to predict the response of a system and suggest improvements. 
- define the key components of a digital twin and identify areas where digital twins can be applied,
- explain the difference between a digital model, a digital shadow, and a digital twin.​

The exercise consists of two sets of multiple-choice questions:

**Part 1** considers general definitions and concepts of digital twins that reinforce the following learning objectives from the micromodule.

**Part 2** contains a simulated interaction with a digital twin. You will need to execute and interact with Python scripts, which are available at the github repository (https://github.com/olgasolodovnikova/dtu-digital-twin-mm.git). You can watch an introductory video covering the installation as well as how you can interact with the code. The files can be evaluated using the included `digital_twin_exercise.ipynb` Jupyter notebook, or executed in your preferred Python IDE. 

## Installation
This guide assumes you have a working installation of Python, the [conda](https://docs.conda.io/en/latest/) package manager and git.

To clone the repository, copy the https url link of the repository (see green code button). 

1. Using the terminal (for Windows users, try the Anaconda prompt), navigate to the directory you want to clone the repository to.

`$ cd Desktop`

2. Clone the git repository to your Desktop.

`$ git clone https://github.com/olgasolodovnikova/dtu-digital-twin-mm.git`

If you are struggling with authentication, I recommend getting [Git Credential Manager](https://github.com/git-ecosystem/git-credential-manager/blob/main/README.md).

3. Make a designated conda environment for this project.

`$ conda create --name MyEnv`

`$ conda activate MyEnv`

Change to `dtu-digital-twin-mm` directory

`$ cd dtu-digital-twin-mm`

4. Install modules from the `requirements.txt` file

`$ conda install --file requirements.txt`

Run `test.py` or `task_c.py` to see that everything is working. 

`$ python3 test.py`

`$ python3 task_c.py`

5. If you wish to use Jupyter notebook, ensure that you install the kernel for your environment by running the following python command.

`$ python3 -m ipykernel install --user --name=MyEnv`

6. You can view and edit the scripts using [vim](https://vim.rtorr.com),

`$ vim task_a.py`

or with Jupyter Lab,

`$ jupyter lab`








