FROM nfcore/base
MAINTAINER Hadrien Gourlé <hadrien.gourle@slu.se>
LABEL authors="hadrien.gourle@slu.se" \
    description="Docker image containing all requirements for nf-core/mag pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/mag-0.1.0/bin:$PATH
