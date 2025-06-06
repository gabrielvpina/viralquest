# Use miniconda3 as the base image
FROM continuumio/miniconda3:latest

# Set environment variables
ENV CONDA_ENV_NAME=viralquest
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=utf-8
ENV TERM=xterm-256color

# Set working directory
WORKDIR /app

# Update conda and install mamba for faster package installation
RUN conda update -n base -c defaults conda && \
    conda install -n base -c conda-forge mamba

# Create conda environment with Python 3.12
RUN mamba create -n ${CONDA_ENV_NAME} python=3.12 -y

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "viralquest", "/bin/bash", "-c"]

# Install bioconda packages (without diamond)
RUN mamba install -n ${CONDA_ENV_NAME} -c bioconda cap3 blast -y

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

# Create a wrapper script for easier execution with proper user
RUN echo '#!/bin/bash\ncd /workspace\nconda run -n viralquest python /app/viralquest.py "$@"' > /usr/local/bin/viralquest-wrapper && \
    chmod +x /usr/local/bin/viralquest-wrapper

# Create a non-root user to avoid permission issues
RUN useradd -m -u 1000 -s /bin/bash viraluser && \
    chown -R viraluser:viraluser /app

# Set default working directory to workspace
WORKDIR /workspace

# Switch to non-root user
USER viraluser

# Set the default command to activate conda environment and start bash
CMD ["conda", "run", "-n", "viralquest", "/bin/bash"]

# Alternative: Set viralquest.py as the default entrypoint
# ENTRYPOINT ["conda", "run", "-n", "viralquest", "python", "viralquest.py"]
