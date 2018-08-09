From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER Hadrien Gourlé <hadrien.gourle@slu.se>
    DESCRIPTION Container image containing all requirements for the nf-core/mag pipeline
    VERSION 0.1.0

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
    /opt/conda/bin/conda clean -a
