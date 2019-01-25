FROM nfcore/base

LABEL maintainer="Hadrien Gourlé <hadrien.gourle@slu.se>"
LABEL description="Docker image containing all requirements for nf-core/mag pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mag-1.0.0/bin:$PATH

RUN mkdir -p /opt/checkm_data && \
    cd /opt/checkm_data && \
    curl -L -O https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz && \
    tar xzf checkm_data_2015_01_16.tar.gz && \
    rm checkm_data_2015_01_16.tar.gz && \
    cd .. && \
    printf "/opt/checkm_data\n/opt/checkm_data\n" | checkm data setRoot
