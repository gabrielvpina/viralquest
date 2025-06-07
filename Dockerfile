FROM continuumio/miniconda3:latest

ENV CONDA_ENV_NAME=viralquest
ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV PYTHONIOENCODING=utf-8
ENV TERM=xterm-256color

WORKDIR /app

RUN conda update -n base -c defaults conda && \
    conda install -n base -c conda-forge mamba

RUN mamba create -n ${CONDA_ENV_NAME} python=3.12 -y

SHELL ["conda", "run", "-n", "viralquest", "/bin/bash", "-c"]

RUN mamba install -n ${CONDA_ENV_NAME} -c bioconda cap3 blast -y

RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.1.12/diamond-linux64.tar.gz && \
    tar -xzf diamond-linux64.tar.gz && \
    chmod +x diamond && \
    mv diamond /opt/conda/envs/${CONDA_ENV_NAME}/bin/ && \
    rm diamond-linux64.tar.gz

RUN apt-get update && apt-get install -y \
    git \
    wget \
    build-essential \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/gabrielvpina/viralquest.git .

RUN conda run -n ${CONDA_ENV_NAME} pip install -r requirements.txt

RUN echo "conda activate ${CONDA_ENV_NAME}" >> ~/.bashrc

RUN conda run -n ${CONDA_ENV_NAME} python viralquest.py --version

RUN echo '#!/bin/bash\ncd /workspace\nconda run -n viralquest python /app/viralquest.py "$@"' > /usr/local/bin/viralquest-wrapper && \
    chmod +x /usr/local/bin/viralquest-wrapper

RUN useradd -m -u 1000 -s /bin/bash viraluser && \
    chown -R viraluser:viraluser /app

WORKDIR /workspace

USER viraluser

CMD ["conda", "run", "-n", "viralquest", "/bin/bash"]
