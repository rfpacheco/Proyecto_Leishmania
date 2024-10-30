# Use a base image with Miniconda
FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the environment.yml file into the container
COPY environment.yml .

# Create the Conda environment from the environment.yml file
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "SIDER_RepetitiveSearcher", "/bin/bash", "-c"]

# Copy the application code into the container
COPY . .

# Ensure the environment is activated whenever a command is executed
# This is necessary because CMD will not have access to the conda environment variables by default
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "SIDER_RepetitiveSearcher"]

# Define the default command to run your application
CMD ["python", "SIDER_RepetitiveSearcher.py"]