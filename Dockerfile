FROM ubuntu

RUN apt-get update && \
    apt-get -y install software-properties-common && \
    apt-get -y install wget && \
    apt-get -y install build-essential

RUN add-apt-repository -y ppa:staticfloat/juliareleases && \
    add-apt-repository -y ppa:staticfloat/julia-deps && \
    apt-get update && \
    apt-get install -y julia && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get -y install libhdf5-7

RUN julia -e 'Pkg.init()' && \
    julia -e 'Pkg.add("Blosc")' && \
    julia -e 'Pkg.add("HDF5")' && \
    julia -e 'Pkg.add("JLD")' && \
    julia -e 'Pkg.add("ArgParse")' && \
    julia -e 'Pkg.add("Options")' && \
    julia -e 'Pkg.add("Distributions")' && \
    julia -e 'Pkg.add("JSON")'

WORKDIR /root/.julia/v0.4/
RUN git clone https://github.com/leviboyles/CPABS.git

WORKDIR /opt
RUN mkdir cpabs
RUN cp /root/.julia/v0.4/CPABS/cpabs_cmd.jl cpabs/
RUN cp /root/.julia/v0.4/CPABS/*.xml cpabs/
