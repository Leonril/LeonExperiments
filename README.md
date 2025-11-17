# LeonExperiments.jl

A personal Julia playground for experiments with geometry, meshes, Finite Element Analysis (FEA), visualization, and the Comodo library.

## Purpose

This repository is meant for:
- Trying out ideas using Julia & Comodo
- Running numerical or geometric experiments
- Testing visualization workflows with GLMakie
- Storing reusable utilities and prototypes
- Creating optimisation methods for my designs and meshes to be additively manufactured

It is **not** meant to be a polished package (yet!), but it may evolve into one.

## Project Setup

To work in this environment:

```julia
using Pkg
Pkg.activate(".")     # activate this project
Pkg.instantiate()     # install dependencies
