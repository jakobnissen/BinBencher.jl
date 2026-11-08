# BinBencher
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jakobnissen.github.io/BinBencherBackend.jl/dev)

BinBencher is a command-line interface (CLI) tool to evaluate metagenomic binnings given a ground-truth reference.
It is the CLI version of the Julia library [BinBencherBackend.jl](https://github.com/jakobnissen/BinBencherBackend.jl).
For more customised, flexible analysis of binnings, it may be preferred to use the backend directly in e.g. an interactive notebook.

Read [the documentation on the BinBencherBackend.jl doc page.](https://jakobnissen.github.io/BinBencherBackend.jl/dev)

## Installation of the CLI program
1. Install the [Julia programming language](https://julialang.org/), version 1.14 or newer.
2. Make sure `julia` is on your `$PATH`, e.g. test it by running `julia -v`
3. Install BinBencher as a Julia app with the following command:

```shell
julia --startup=no -e 'using Pkg; Pkg.Apps.add(url="https://github.com/jakobnissen/BinBencher.jl")'
```

You will now have an executable called `binbench` in the `bin` subdirectory your Julia home directory.
For example, on my computer the executable is at `~/.julia/bin/binbench`.
Add the `bin` folder to your `$PATH` environmental variable to be able run to BinBencher as `binbench` from the shell.
Run `binbench --help` to see the available subcommands.

## Running on example files
1. Go to [the releases page](https://github.com/jakobnissen/BinBencher.jl/releases/tag/filesv1) and download the file `inputs.tar.gz`.
2. Extract the tar archive with `tar -xzvf inputs.tar.gz`
3. Run BinBencher with the command `binbench bench out reference.json bins.tsv`

## Existing datasets
* CAMI2 toy human microbiome short read gold standard assembly: https://zenodo.org/records/15083711
* NCBI file: https://github.com/jakobnissen/BinBencher.jl/releases/download/filesv1/ncbi.tsv.gz
