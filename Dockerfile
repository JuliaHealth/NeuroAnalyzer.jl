# syntax=docker/dockerfile:1

FROM julia:latest

LABEL maintainer="Adam Wysokiński <adam.wysokinski@umed.lodz.pl>"

COPY misc/install.jl /

# install dependencies
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get install -y --no-install-recommends libgl1-mesa-dev libglu1-mesa-dev freeglut3-dev libxt6 libxrender1 libxext6 libgl1-mesa-glx libqt5widgets5
RUN apt-get install -y --no-install-recommends xvfb xauth
RUN apt-get install -y --no-install-recommends git hdf5-tools
RUN rm -rf /var/lib/apt/lists/*

RUN mkdir ~/NeuroAnalyzer
RUN mkdir -p ~/.julia/config

ENV DISPLAY=:0
ENV install='using Pkg; Pkg.add("Revise"); Pkg.add(url="https://codeberg.org/AdamWysokinski/FIRLSFilterDesign.jl"); Pkg.add("NeuroAnalyzer"); Pkg.instantiate(); Pkg.resolve(); Pkg.update()'
RUN xvfb-run -s '-screen 0 1024x768x24' julia -q --color=yes -e "$install"

COPY misc/startup.jl /root/.julia/config

CMD ["bash", "-c", "julia -q --color=yes"]
