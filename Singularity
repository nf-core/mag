From:nfcore/base
Bootstrap:docker

%labels
=======
    MAINTAINER Hadrien Gourlé <hadrien.gourle@slu.se>
    DESCRIPTION Singularity image containing all requirements for nf-core/mag pipeline
    VERSION 0.1.0dev

%environment
    PATH=/opt/conda/envs/nf-core-mag-0.1.0dev/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a
