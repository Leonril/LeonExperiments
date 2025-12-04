# LeonExperiments.jl

A personal Julia playground for experiments with geometry, meshes, Finite Element Analysis (FEA), visualization, and the Comodo library.

## Purpose

This repository is meant for:
- Trying out ideas using Julia & Comodo
- Running Finite Element Analysis simulations using open source FEA software FEBio
- Testing visualization workflows with GLMakie
- Storing reusable utilities and prototypes
- Creating optimisation methods for my designs and meshes to be additively manufactured

It is **not** meant to be a polished package, but it may evolve into one.

## Currently Working On
- [ ] Creating 3D meshes using levelsets and isosurfaces
- [ ] Implementing topology optimisation
- [ ] A reverse finite element demo targeting both shape and/or material properties
- [ ] Creating a 3D mesh from image data (CT scan etc.)

## Project Setup

To work in this environment:

```julia
using Pkg
Pkg.activate(".")     # activate this project
Pkg.instantiate()     # install dependencies
