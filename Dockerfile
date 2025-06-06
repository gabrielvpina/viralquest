# Use miniconda3 as the base image
FROM continuumio/miniconda3:latest

# Set environment variables
ENV CONDA_ENV_NAME=viralquest
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1

# Set working directory
WORKDIR /app

# Update conda and install mamba for faster package installation
RUN conda update -n base -c defaults conda && \
    conda install -n base -c conda-forge mamba

# conda environment with Python 3.12
RUN mamba create -n ${CONDA_ENV_NAME} python=3.12 -y

# RUN commands use the new environment
SHELL ["conda", "run", "-n", "viralquest", "/bin/bash", "-c"]

# bioconda packages
RUN mamba install -n ${CONDA_ENV_NAME} -c bioconda cap3 blast diamond=2.1.11 -y

# git if not already available
RUN apt-get update && apt-get install -y git && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/gabrielvpina/viralquest.git .

# requirements.txt
RUN conda run -n ${CONDA_ENV_NAME} pip install -r requirements.txt

# environment activated by default
RUN echo "conda activate ${CONDA_ENV_NAME}" >> ~/.bashrc

# Test
RUN conda run -n ${CONDA_ENV_NAME} python viralquest.py --version

# default command to activate conda environment and start bash
CMD ["conda", "run", "-n", "viralquest", "/bin/bash"]

# Alternative: run viralquest directly, uncomment the following line:
# ENTRYPOINT ["conda", "run", "-n", "viralquest", "python", "viralquest.py"]
