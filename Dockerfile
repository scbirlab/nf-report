FROM mambaorg/micromamba:1.5.6

ENV HOME=/home/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR $HOME

# Add a dummy user entry to avoid `passwd` errors in Docker
USER root
RUN echo "user:x:1001:1001::/home/user:/bin/bash" >> /etc/passwd && \
    mkdir -p /home/user && chown -R 1001:1001 /home/user

RUN mkdir -p $HOME/.conda && chown -R 1000:1000 $HOME
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    gcc \
    g++ \
    perl liburi-perl \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*
USER 1000

COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes
RUN mkdir -p /tmp/{xdg,fontconfig,mpl,numba}

# force OS binaries to take precednt over mamba (curl issues)
ENV PATH=$PATH:$MAMBA_ROOT_PREFIX/envs/env/bin
ENV MPLBACKEND=Agg
ENV XDG_CACHE_HOME=/tmp/xdg
ENV MPLCONFIGDIR=/tmp/mpl
ENV NUMBA_CACHE_DIR=/tmp/numba
