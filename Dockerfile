FROM nfcore/base:1.7

LABEL authors="Hadrien Gourlé <hadrien.gourle@slu.se>, Daniel Straub <d4straub@gmail.com>" \
    description="Docker image containing all requirements for nf-core/mag pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-mag-1.0.0/bin:$PATH
