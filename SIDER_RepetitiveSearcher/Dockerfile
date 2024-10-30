# Use a base image with Miniconda or Anaconda
FROM continuumio/miniconda3

# Set up the working directory
WORKDIR /usr/src/app

# Copy the environment.yml to the Docker image
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "SIDER_RepetitiveSearcher", "/bin/bash", "-c"]

# Install pip packages
RUN conda run -n SIDER_RepetitiveSearcher pip install numpy pandas pytz six wheel

# Copy your application code to the Docker image
COPY . .

# Define the default command to run your application
CMD ["conda", "run", "--no-capture-output", "-n", "SIDER_RepetitiveSearcher", "python", "your_script.py"]FROM ubuntu:latest
LABEL authors="viskuit"

ENTRYPOINT ["top", "-b"]