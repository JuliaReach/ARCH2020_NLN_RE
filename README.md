# ARCH2020_NLN_RE

This is the repeatability evaluation package for the
ARCH-COMP20 Category Report: Continuous and Hybrid Systems with Nonlinear Dynamics
of the 4th International
Competition on Verifying Continuous and Hybrid Systems Friendly Competition
(ARCH-COMP '20).

## Installation

To build the container, you need the program `docker`.
For installation instructions on different platforms, consult
[the Docker documentation](https://docs.docker.com/install/).
For general information about `Docker`, see
[this guide](https://docs.docker.com/get-started/).

Once you have installed Docker, clone this repository and build the container:

```shell
$ git clone https://github.com/JuliaReach/ARCH2020_NLN_RE.git
$ cd ARCH2020_NLN_RE
$ docker build . -t juliareach
```

To run the container with the benchmarks, type:

```shell
$ docker run -i juliareach
```

Or to run it interactively, type:

```shell
$ docker run -it juliareach bash

$ julia

julia> include("startup.jl")
```

## Models

**Important Note.** The `master` branch in this repository provides a `Manifest.toml` file that is meant to be a
"screenshot" of the package ecosystem for this RE. If you intend to use
`JuliaReach` for other purposes than RE, we strongly recommend that you instead
use the current long-term-support Julia version and follow the installation
instructions in
[ReachabilityAnalysis.jl](https://github.com/JuliaReach/ReachabilityAnalysis.jl).

## Plots

The plots are stored in the main folder as `*.png` files.
To obtain these files, run the image interactively, specify an output volume,
and in the end copy the files to the output volume.
In the example below we call the output volume `result`.

```shell
$ docker run -it -v result:/result juliareach bash

$ julia -e 'include("startup.jl")'

$ cp *.png /result
```

Then one can access the plots via the following command from outside:

```shell
$ docker cp FANCY_NAME:/result/ .
```

Here `FANCY_NAME` is the name of the Docker container, which can be found via:

```shell
$ docker container ls --all
```

Here is an example of the output of that command:

```shell
CONTAINER ID        IMAGE               COMMAND             CREATED             STATUS                     PORTS               NAMES
1f661e604079        juliareach          "bash"              15 minutes ago      Exited (0) 2 minutes ago                       upbeat_ride
```

## About this RE

This repeatability evaluation package relies on Julia's package manager, `Pkg`, to create an environment that can be used to exactly use the same package versions in different machines. For further instructions on creating this RE, see [this wiki](https://github.com/JuliaReach/ARCH2020_NLN_RE/wiki/Instructions-for-creating-this-RE).
