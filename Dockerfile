FROM ubuntu

RUN sudo add-apt-repository ppa:staticfloat/juliareleases && apt-get update && apt-get install julia
RUN sudo apt-get install libhdf5-7 hdf5-tools

RUN julia -e 'Pkg.init()'
RUN julia -e 'Pkg.add("HDF5")'
RUN julia -e 'Pkg.add("JLD")'
RUN julia -e 'Pkg.add("ArgParse")'
RUN julia -e 'Pkg.add("Options")'
RUN julia -e 'Pkg.add("Distributions")'
RUN julia -e 'Pkg.add("JSON")'

WORKDIR ~/.julia/v0.4/
RUN git clone https://github.com/leviboyles/CPABS.git

WORKDIR /opt
RUN mkdir cpabs
RUN cp ~/.julia/v0.4/CPABS/cpabs_cmd.jl cpabs/
RUN cp ~/.julia/v0.4/CPABS/*.xml cpabs/
