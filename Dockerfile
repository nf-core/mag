FROM nfcore/base:1.9

LABEL authors="Hadrien Gourlé <hadrien.gourle@slu.se>, Daniel Straub <d4straub@gmail.com>" \
    description="Docker image containing all requirements for nf-core/mag pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-mag-1.0.0 > nf-core-mag-1.0.0.yml
ENV PATH /opt/conda/envs/nf-core-mag-1.0.0/bin:$PATH
