# demo_cube.jl
# Description: A demo script to create and display a cube using Comodo.

#using LeonExperiments
using Comodo
using GLMakie

F,V = Comodo.cube(7)
fig = Figure()
ax=Axis3(fig[1,1])
Comodo.meshplot!(ax, F, V)
display(fig)