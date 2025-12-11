# demo_stl.jl
# Description: A demo script to create and save stl of test shapes.


# Future Improvements to Add: 
# - Finish batman logo example
#      - Batman curve working 
#      - Need to make 3D and add tetgen 
# - Finish holly Christmas ornament example
#      - Just beginning

# - Add more shape variants for demo


#using LeonExperiments
using Comodo
using GLMakie
using GeometryBasics
using MeshIO
using FileIO
using FilePathsBase

GLMakie.closeall() #Close any open Makie windows

testCase = 2# Change this value to switch between test cases
#testCase = 1: Create and save an stl of a cube
#testCase = 2: Create and save an stl of an extruded batman logo from Comodo
#testCase = 3: Create and save an stl of a Christmas ornament using an overengineered levelset method, for practice


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
    # Generating Batman curve
    Vb=batman(50) #List of nodes of the Batman logo

    #scaling curve
    width = 40.0 # desired width of batman logo in mm
    Vb = Vb .* (width / 2) # Scaling, default width is 2 mm, so dividing by 2

    # Visualling Vb Batman curve
    fig = Figure(size=(1000,500))
    ax1 = AxisGeom(fig[1, 1]; title="Boundary Vertices", azimuth=-pi/2, elevation=pi/2)
    hp1 = lines!(ax1, Vb,linewidth=3,color=:blue)
    hp2 = scatter!(ax1, Vb,markersize=8,color=:red)
    screen = display(GLMakie.Screen(), fig)

    # Setting up regiontrimesh inputs
    VT = (copy(Vb),) # Tuple of the vertex array, regiontrimesh seems to mutate Vb so use a copy
    R = ([1],) # Tuple of the region to be meshed
    P = pointspacingmean(copy(Vb); close_loop=true) # Average spacing of the points on the curve
    
    # Creating triangular mesh of batman logo
    F,Vbottom,C = regiontrimesh(VT,R,P)
   
    # Plotting triangular mesh
    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig2[1, 1], title="Triangular mesh", azimuth=-pi/2, elevation=pi/2)
    hp1 = meshplot!(ax1, F, Vbottom, strokewidth=1.0, color=:lightblue)    
    screen = display(GLMakie.Screen(), fig)

    # Copying mesh for top of extrude
    height = 5.0
    Vtop = [Point(v[1], v[2], v[3] + height) for v in Vbottom]
    Ctop = C .+ 1 # Shift color indices for top faces

    #Plotting top and bottom meshes
    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Top and bottom faces", azimuth=-pi/2, elevation=pi/2)
    hp1 = meshplot!(ax1, F, Vbottom, strokewidth=1.0, color=:red)
    hp1 = meshplot!(ax1, F, Vtop, strokewidth=1.0, color=:lightblue)
    screen = display(GLMakie.Screen(), fig)

    # shift Ftop
    offset = length(Vbottom) # number of vertices in bottom face
    Ftop = [ # shifting indices of top face triangles
        TriangleFace(f[1] + offset, f[2] + offset, f[3] + offset)
    for f in F]

    # Finding the boundary edges for extrude sides
    indexmap = Dict(Vbottom[i] => i for i in eachindex(Vbottom))

    Vb_bottomInd = [ get(indexmap, p, nothing) for p in Vb ] # Find indexes of boundary points Vb in Vbottom
    Vb_bottomInd = filter(!isnothing, Vb_bottomInd) # Remove any missing points

    Vb_topInd = [i + offset for i in Vb_bottomInd] #Indexes of boundary points for top of the curve 
    

    # Creating side faces of extrude
    n = length(Vb) # number of vertices in bottom face
    Fs1 = [TriangleFace(Vb_bottomInd[i], Vb_bottomInd[i+1], Vb_topInd[i]) for i in 1:(n-1)] # one half of quad split into two tris
    Fs2 = [TriangleFace(Vb_topInd[i], Vb_topInd[i+1], Vb_bottomInd[i+1]) for i in 1:(n-1)] # other half of quad split into two tris
    
    push!(Fs1, TriangleFace(Vb_bottomInd[n], Vb_bottomInd[1], Vb_topInd[n])) # closing the loop of faces
    push!(Fs2, TriangleFace(Vb_topInd[n], Vb_topInd[1], Vb_bottomInd[1]))
    
    C_sides = fill(maximum(Ctop)+1, length(Fs1) + length(Fs2)) # New color for side faces

    # combine geometry
    V = [Vbottom; Vtop]
    Fb = [F; Ftop; Fs1; Fs2]
    Cb = [C; Ctop; C_sides]

    # Plotting full extruded mesh
    Fp,Vp = separate_vertices(Fb,V) # Give each face its own point set 
    Cp = simplex2vertexdata(Fp,Cb) # Convert face color data to vertex color data 

    fig4 = Figure(size=(1200,1000))
    ax4 = AxisGeom(fig4[1, 1], title="Batman mesh 3D", azimuth=-pi/2, elevation=pi/2)
    hp4 = meshplot!(ax4, Fp, Vp, strokewidth=1.0, color=Cp, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))    
    Colorbar(fig4[1, 1][1, 2], hp4)
    screen = display(GLMakie.Screen(), fig4)

    #convert to triangular mesh
    surfaceMesh=GeometryBasics.Mesh(Point.(V),Fb)

    #Save as .stl
    export_dir = joinpath(@__DIR__, "..", "exports") #Adds exports folder in parent directory
    isdir(export_dir) || mkdir(export_dir) #Create exports folder if it doesn't exist
    File_name = joinpath(export_dir,"demo_stl_batman.stl") #Define file path and name
    save(File_name, surfaceMesh); #Save stl

elseif testCase == 3
    # Creating holly berries 
    FBerryBase, VBerryBase = Comodo.geosphere(3,10) # Creates a sphere mesh

    fig = Figure(size=(1000,500))
    ax=Axis3(fig[1,1])
    Comodo.meshplot!(ax, FBerryBase, VBerryBase)
    display(fig)

    # Creating the leave levelsets

    
    
    #generate holly ornament levelset
    #convert to triangular mesh
    #plot holly ornament
    #convert to stl
    #save stl

end