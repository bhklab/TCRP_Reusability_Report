# hash:sha256:ebed32e72fa6f36d267b175ab20565476b118fefa9487ef56f90b42f4c214533
FROM registry.codeocean.com/codeocean/miniconda3:4.9.2-cuda11.7.0-cudnn8-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.8ubuntu1.1 \
        curl=7.68.0-1ubuntu2.15 \
        libgit2-dev=0.28.4+dfsg.1-2 \
        libjpeg-dev=8c-2ubuntu8 \
        libssl-dev=1.1.1f-1ubuntu2.16 \
        pandoc=2.5-3build2 \
        pkg-config=0.29.1-0ubuntu4 \
        wget=1.20.3-1ubuntu2 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        jupyterlab==3.5.2 \
    && conda clean -ya

RUN pip install -U --no-cache-dir \
    matplotlib==3.6.2 \
    mygene==3.2.2 \
    networkx==3.0 \
    pandas==1.5.2 \
    pillow==9.4.0 \
    scikit-learn==1.2.0 \
    scipy==1.10.0 \
    seaborn==0.12.2 \
    torch==1.13.1

ADD "https://github.com/coder/code-server/releases/download/v4.9.0/code-server-4.9.0-linux-amd64.tar.gz" /.code-server/code-server.tar.gz
	
RUN cd /.code-server \
	&& tar -xvf code-server.tar.gz \
	&& rm code-server.tar.gz \
	&& ln -s /.code-server/code-server-4.9.0-linux-amd64/bin/code-server  /usr/bin/code-server
