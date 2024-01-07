# dtu-digital-twin-mm

# What is a Digital Twin?

After completing this exercise, you should be able to
- interact with a digital twin in order to predict the response of a system and suggest improvements. 

The exercise consists of two sets of multiple-choice questions:

**Part 1** considers general definitions and concepts of digital twins that reinforce the following learning objectives from the micromodule;
- define the key components of a digital twin and identify areas where digital twins can be applied,
- explain the difference between a digital model, a digital shadow, and a digital twin.​

**Part 2** contains a simulated interaction with a digital twin. You will need to execute and interact with Python scripts from this repository. You can watch an introductory video covering the installation as well as how you can interact with the code. The files can be evaluated using the included `digital_twin_exercise.ipynb` Jupyter notebook, or executed in your preferred Python IDE.

## Installation
This guide assumes you have a working installation of Python, the [conda](https://docs.conda.io/en/latest/) package manager and git.

To clone the repository, copy the https url link of the repository (see green <code> button). 

Using the terminal, navigate to the directory you want to clone the repository to, .e.g open your terminal and run the following commands for installing the repository on your Desktop.
>`cd Desktop`

>`git clone https://github.com/olgasolodovnikova/dtu-digital-twin-mm.git`

If you are struggling with authentication, I recommend getting [Git Credential Manager](https://github.com/git-ecosystem/git-credential-manager/blob/main/README.md).

Make a designated conda environment for this project.
>`conda create --name MyEnv`
>`conda activate MyEnv`
>`conda install modules_req.txt`

Change to `dtu-digital-twin-mm` directory
>`cd dtu-digital-twin-mm`

Run `test.py` or `task_c.py` to see that everything is working. 

>`python3 test.py`
>`python3 task_c.py`

If you wish to use Jupyter notebook, ensure that you activate the kernel for your environment. First install ipykernel (it's bad form to install something after installing from .txt file, so consider adding `ipykernel` to `modules_req.txt`. 
>conda install ipykernel
>python3 -m ipykernel install --user --name=MyEnv





