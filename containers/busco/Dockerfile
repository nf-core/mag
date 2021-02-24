FROM nfcore/base:1.12.1
LABEL authors="Hadrien Gourlé <hadrien.gourle@slu.se>, Daniel Straub <d4straub@gmail.com>, Sabrina Krakau <sabrinakrakau@gmail.com>" \
      description="Docker image containing BUSCO requirements for the nf-core/mag pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create --quiet -f /environment.yml && conda clean -a

# For BUSCO to avoid segfault error with R 4.0.2
RUN apt-get update
RUN apt-get install -y libxt6

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-mag-busco-1.3.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-mag-busco-1.3.0dev > nf-core-mag-busco-1.3.0dev.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron
