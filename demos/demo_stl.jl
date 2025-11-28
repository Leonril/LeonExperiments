# demo_stl.jl
# Description: A demo script to create and save stl of test shapes.


# Future Improvements to Add: 
# - Add more shape variants for demo


#using LeonExperiments
using Comodo
using GLMakie
using GeometryBasics
using MeshIO
using FileIO

testCase = 1 # Change this value to switch between test cases
#testCase = 1: Create and save an stl of a cube
#testCase = 2: Create and save an stl of an extruded batman logo from Comodo



if testCase == 1
    #generate cube as tet mesh
    E,V,Fb,Cb = Comodo.tetbox((5,5,5),2.5) 
    
    #plotting
    fig = Figure()
    ax=Axis3(fig[1,1])
    Comodo.meshplot!(ax, Fb, V)
    display(fig)

    #convert to triangular mesh
    surfaceMesh=GeometryBasics.Mesh(Point.(V),Fb)

    #Save as .stl
    save("C:/Users/riley/Docs PC/GitHub/LeonExperiments.jl/exports/demo_stl_cube.stl", surfaceMesh)
    

    #generate F,V of cube 
    #plot cube
    #convert to triangular mesh if needed
    #plot again as triangular mesh
    #convert to stl
    #save stl

elseif testCase == 2
    #generate F,V of batman 
    #plot batman
    #convert to triangular if needed
    #plot again as triangular mesh
    #convert to stl
    #save stl
end