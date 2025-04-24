# CRISPR_guides_design
Automated pipeline for customizable CRISPRa and CRISPRi guides design

# environment set up and installation
All the installations are provided in the docker published on our repository as hedgelab/crispr:image1.
To use this pipeline is highly recommended to use it through the docker image to ensure reproducibility, alernatively it is possible to manually repeat the installation on your computer following the steps in the Dockerfile (optimised for a linux bash shell).

If you do not have Docker installed on your computer yet, you can follow the instructions in the following link: https://docs.docker.com/get-started/get-docker/.
Ensure that the RAM and CPU allocated to Docker are sufficient to run the analysis. If you wish to run the analysis on an already indexed genome, 10 GB of RAM will be enough. Else, if you will need to run the indexing step, at least 30 GB are needed. To set the proper RAM limits for your analysis, some versions of Docker Desktop allow you to directly edit these parameters in the settings, while others require the use of an external settings file. In the latter case, create a .wslconfig file in the User/user_name folder of your computer. You can follow the example template provided. After making these changes, restart the Docker engine to apply them effectively.

Once Docker is installed and properly configured, activate the Docker engine directly opening the Docker Desktop application or by running on your terminal: start docker
