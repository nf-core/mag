From:nfcore/base
Bootstrap:docker

%labels
=======
    MAINTAINER Hadrien Gourlé <hadrien.gourle@slu.se>
    DESCRIPTION Singularity image containing all requirements for nf-core/mag pipeline
    VERSION 1.0.0

%environment
    PATH=/opt/conda/envs/nf-core-mag-1.0.0/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
