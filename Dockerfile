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

# Create conda environment with Python 3.12
RUN conda create -n ${CONDA_ENV_NAME} python=3.12 -y

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "viralquest", "/bin/bash", "-c"]

# Install bioconda packages (without diamond)
RUN conda install -n ${CONDA_ENV_NAME} -c bioconda cap3 blast -y

# Download and install Diamond executable
RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.12/diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && \
    chmod +x diamond && \
    mv diamond /opt/conda/envs/${CONDA_ENV_NAME}/bin/ && \
    rm diamond-linux64.tar.gz

# Install git, wget, and build tools
RUN apt-get update && apt-get install -y \
    git \
    wget \
    build-essential \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Clone the repository
RUN git clone https://github.com/gabrielvpina/viralquest.git .

# Install pip packages from requirements.txt
RUN conda run -n ${CONDA_ENV_NAME} pip install -r requirements.txt

# Make sure the conda environment is activated by default
RUN echo "conda activate ${CONDA_ENV_NAME}" >> ~/.bashrc

# Test the installation
RUN conda run -n ${CONDA_ENV_NAME} python viralquest.py --version

# Set the default command to activate conda environment and start bash
CMD ["conda", "run", "-n", "viralquest", "/bin/bash"]

# Alternative: If you want to run viralquest directly, uncomment the following line:
# ENTRYPOINT ["conda", "run", "-n", "viralquest", "python", "viralquest.py"]
