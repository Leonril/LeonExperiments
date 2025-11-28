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
using FilePathsBase


testCase = 2 # Change this value to switch between test cases
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
    export_dir = joinpath(@__DIR__, "..", "exports") #Adds exports folder in parent directory
    isdir(export_dir) || mkdir(export_dir) #Create exports folder if it doesn't exist
    File_name = joinpath(export_dir,"demo_stl_cube.stl") #Define file path and name
    save(File_name, surfaceMesh); #Save stl

elseif testCase == 2
    Vb=batman(50)


    fig = Figure(size=(1000,500))

    ax1 = AxisGeom(fig[1, 1]; title="stepwise=true, approximately n points, anti-clockwise", azimuth=-pi/2, elevation=pi/2)
    hp1 = lines!(ax1, Vb,linewidth=3,color=:blue)
    hp2 = scatter!(ax1, Vb,markersize=8,color=:red)
    hp2 = scatter!(ax1, Vb[1],markersize=15,color=:yellow)
    hp2 = scatter!(ax1, Vb[2],markersize=15,color=:orange)
    fig
    #generate F,V of batman 
    #plot batman
    #convert to triangular if needed
    #plot again as triangular mesh
    #convert to stl
    #save stl
end