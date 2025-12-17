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
using LinearAlgebra

GLMakie.closeall() #Close any open Makie windows

testCase = 3 # Change this value to switch between test cases
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
    offset = length(Vbottom) # number of vertices in bottom tri mesh faces
    Ftop = [ # shifting indices of top face triangles
        TriangleFace(f[1] + offset, f[2] + offset, f[3] + offset)
    for f in F]

    # Finding the boundary edges for extrude sides, regiontrimesh does not return boundary vertex order
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
    rBerry = 5.0
    pointSpacingBerry = 1.0
    z_threshold = -(rBerry*0.85) # Height at which berry meets base
    FBerry, VBerry = geosphere(2,rBerry) # Creates a sphere mesh

    # Plotting berry tri sphere
    fig = Figure(size=(1000,500))
    ax=Axis3(fig[1, 1], title="Berry (tri sphere)", azimuth=-pi/2, elevation=pi/2)
    Comodo.meshplot!(ax, FBerry, VBerry)
    dscreen = display(GLMakie.Screen(), fig)

    # Remove lower part of the berry to intergrate with base
    idx = findall(getindex.(VBerry, 3) .< z_threshold) # Find indices of verticecs below z_threshold
    setRemove = Set(idx) # Converting to a set  
    F_keep = [f for f in FBerry if all(v -> !(v in setRemove), f)] # Using set to remove faces with vertices below threshold

    # Plotting berry with bottom removed
    fig = Figure(size=(1000,500))
    ax=Axis3(fig[1, 1], title="Berry jagged cut ", azimuth=-pi/2, elevation=pi/2)
    Comodo.meshplot!(ax, F_keep, VBerry)
    dscreen = display(GLMakie.Screen(), fig)

    # Even out z height of bottom of berry edges
    F_keep_cleaned, VBerryTemplate  = remove_unused_vertices(F_keep, VBerry) # Remove unused vertices after face removal
    FaceConnectivity = con_vertex_face(F_keep_cleaned, VBerryTemplate)# List of how many faces touch each vertex

    # Suppose connectivityV[i] is a vector of connected vertices for vertex i
    edge_indices = [i for i in 1:length(VBerryTemplate) if length(FaceConnectivity[i]) < 5] #Edges have less than 5 connected faces

    for i in edge_indices # Update z value of edge vertices from cut faces
    v = VBerryTemplate[i]
    VBerryTemplate[i] = Point(v[1], v[2], z_threshold)  # create a new Point with updated z
    end

    VBerry_Edge = VBerryTemplate[edge_indices] # Store edge vertices 

    # Plotting berry with cleaned bottom edges
    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Berry  with cleaned cut ", azimuth=-pi/2, elevation=pi/2)
    hp1 = Comodo.meshplot!(ax, F_keep_cleaned, VBerryTemplate, color=:red)
    hp2 = scatter!(ax, VBerry_Edge,markersize=10,color=:blue)
    dscreen = display(GLMakie.Screen(), fig)

    # Creating the leaf boundary
    rLeaf1 = 8.0 # radius of the curves of the leaves
    rLeaf2 = 5.0 # radius of the smaller curves of the leaves
    nLeafPoints = 40 # refinement control of the mesh, use multiple of 4 for symmetry
    cut = rLeaf2/5  # cut to make smaller curve fit regularly 
    errror = 0.01  # small error tolerance for point inclusion checks

    V_leaf_circle1 = circlepoints(rLeaf1,nLeafPoints)
    V_leaf_circle2 = circlepoints(rLeaf2,nLeafPoints)
    
    Varc1_unsorted = [v for v in V_leaf_circle1 if v[1] > 0-error && v[2] < 0+error]
    Varc1 = [Varc1_unsorted[2:end]; Varc1_unsorted[1:1]] # Sorting points to go counter clockwise
    Varc2 = [v for v in V_leaf_circle2 if v[2] > 0+cut] # Cutting the points for a curve from a circle

    # Plotting leaf arcs
    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Leaf arcs ", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Varc1 ,markersize=10,color=:blue)
    hp2 = scatter!(ax, Varc2 ,markersize=10,color=:red)
    dscreen = display(GLMakie.Screen(), fig)

    # Using the arcs to create the leaf shape
    VLeafCurve2 = reverse([Point(v[1] + (rLeaf1+rLeaf2), -v[2], v[3]) for v in Varc2])
    VLeafCurve3 = reverse([Point(-v[1] + 2*(rLeaf1+rLeaf2), v[2], v[3]) for v in Varc1])
    
    VLeafHalf_unshifted = [Varc1; VLeafCurve2; VLeafCurve3] # Combine arcs to make half leaf
    VLeafHalf1 = [Point(v[1], v[2]+rLeaf1, v[3]) for v in VLeafHalf_unshifted] # Shift leaf up to center around origin
    VLeafHalf2 = reverse([Point(v[1], -v[2], v[3]) for v in VLeafHalf1])
    VLeaf = [VLeafHalf1; VLeafHalf2[2:end-1]] # Complete leaf by mirroring half and trimming start and stop

    # Plotting leaf arcs
    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Leaf arcs ", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Varc1 ,markersize=10,color=:blue)
    hp2 = scatter!(ax, VLeafCurve2 ,markersize=10,color=:red)
    hp3 = scatter!(ax, VLeafCurve3 ,markersize=10,color=:purple)
    dscreen = display(GLMakie.Screen(), fig)

    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Leaf V", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, VLeaf ,markersize=10,color=:darkgreen)
    dscreen = display(GLMakie.Screen(), fig)

    # Creating the base boundary
    spacingBerry = 2.0 # spacing between berries
    dist_berryCentres = 2*rBerry + spacingBerry # distance between berry centers
    
    # Radius of the circular base to fit three berries
    radiusBerryCentre= ceil((sqrt(3)/3) * dist_berryCentres)

    # first vertex at negative x-axis
    centreBerry1 = [-radiusBerryCentre, 0.0]

    # rotation matrix function
        function rotate(v, θ)
        [cos(θ) -sin(θ);
         sin(θ)  cos(θ)] * v
    end

    θ = 2π/3
    centreBerry2 = rotate(centreBerry1, θ)
    centreBerry3 = rotate(centreBerry1, 2θ)

    radiusBase = radiusBerryCentre + rBerry + spacingBerry 

    # Trimming leaf and base to fit together
    VLeafShifted = [Point(v[1]+0.8*radiusBase, v[2], v[3]) for v in VLeaf] # Shift leaf down to original position
    VLeafTrimmed = [v for v in VLeafShifted if norm(v) > radiusBase]

    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Leaf V, positioned and trimmed for base", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, VLeafTrimmed ,markersize=10,color=:darkgreen)
    dscreen = display(GLMakie.Screen(), fig)

    #Creating base arc points
    VLeadEnd = VLeafTrimmed[end]
    angleLeaf = atan.(VLeadEnd[2], VLeadEnd[1])-0.02 # small offset to ensure no overlap
    VBaseArc = reverse([Point(cos(θ)*radiusBase, sin(θ)*radiusBase, 0.0) for θ in ((2*(2π/3))-angleLeaf):pi/30:2π+angleLeaf])

    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Leaf and Base before rotation", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, VLeafTrimmed ,markersize=10,color=:darkgreen)
    hp2 = scatter!(ax, VBaseArc ,markersize=10,color=:red)
    dscreen = display(GLMakie.Screen(), fig)

    # Using rotation to close the base shape
    VBase1 = [VLeafTrimmed;VBaseArc]
    
    R2 = [cos(θ) -sin(θ); # Rotation matrix 
     sin(θ)  cos(θ)]

    # Rotate all points
    VBase2 = [Point(R * [p[1]; p[2]]..., p[3]) for p in VBase1]
    VBase3 = [Point(R * [p[1]; p[2]]..., p[3]) for p in VBase2]

    # Combine all base points
    Vb_Base = [VBase1; VBase2; VBase3] # Complete base points

    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Base Boundary Points", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Vb_Base ,markersize=10,color=:darkgreen)
    dscreen = display(GLMakie.Screen(), fig)

    # Positioning berries on base
    VBerry1 = [Point(v[1]+centreBerry1[1], v[2]+centreBerry1[2], v[3]+abs(z_threshold)) for v in VBerryTemplate] # Shift Berrys to new positions
    VBerry2 = [Point(v[1]+centreBerry2[1], v[2]+centreBerry2[2], v[3]+abs(z_threshold)) for v in VBerryTemplate] # Shift Berrys to new positions
    VBerry3 = [Point(v[1]+centreBerry3[1], v[2]+centreBerry3[2], v[3]+abs(z_threshold)) for v in VBerryTemplate] # Shift Berrys to new positions
    
    VBerry1_Edge =  VBerry1[edge_indices] # Store edge vertices for regiontrimesh
    VBerry2_Edge =  VBerry2[edge_indices]
    VBerry3_Edge =  VBerry3[edge_indices]

    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Base Boundary Points", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Vb_Base ,markersize=10,color=:darkgreen)
    hp2 = scatter!(ax, [VBerry1_Edge;VBerry2_Edge;VBerry3_Edge] ,markersize=10,color=:blue)
    hp3 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry1, color=:red)
    hp4 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry2, color=:red)
    hp5 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry3, color=:red)
    dscreen = display(GLMakie.Screen(), fig)

    # Put hole in base for string
    HoleSpacing = 2.0 # distance betweeen two HoleSpacing
    HoleRadius = 2.0 # radius of hole
    V_HoleTemplate = circlepoints(HoleRadius,20)

    VHole = [Point(v[1] + radiusBase - (2*HoleRadius), v[2], v[3]) for v in V_HoleTemplate]
    
    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Base Boundary Points", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Vb_Base ,markersize=10,color=:darkgreen)
    hp2 = scatter!(ax, [VBerry1_Edge;VBerry2_Edge;VBerry3_Edge] ,markersize=10,color=:blue)
    hp3 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry1, color=:red)
    hp4 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry2, color=:red)
    hp5 = Comodo.meshplot!(ax, F_keep_cleaned, VBerry3, color=:red)
    hp6 = scatter!(ax, VHole ,markersize=10,color=:purple)
    dscreen = display(GLMakie.Screen(), fig)

    # Making underside of ornament
    HollyThickness = 4.0
    testCircle = circlepoints(radiusBase,40)
    VtestCircle_Underside = [Point(v[1], v[2], -HollyThickness) for v in testCircle]
    Vb_Base_Underside = [Point(v[1], v[2], v[3]-HollyThickness) for v in Vb_Base]
    VHole_Underside = [Point(v[1], v[2], v[3]-HollyThickness) for v in VHole]

    # Setting up regiontrimesh inputs
    #VT_top = [copy(Vb_Base),copy(VBerry1_Edge),copy(VBerry2_Edge),copy(VBerry3_Edge),copy(VHole),] # Tuple of the vertex array, regiontrimesh seems to mutate Vb so use a copy
    VT_bottom = [copy(Vb_Base_Underside),]
    R_top = ([1,2,3,4,5],)
    R_bottom = ([1],)
    PointSpacing = (1.0)

    # regiontrimesh
    #F_top,V_top,C_top = regiontrimesh(VT_top,R_top,PointSpacing)
    F_bottom,V_bottom,C_bottom = regiontrimesh(VT_bottom,R_bottom,PointSpacing)
    
    # Plotting 
    fig = Figure(size=(1000,500))
    ax=AxisGeom(fig[1, 1], title="Base Boundary Points", azimuth=-pi/2, elevation=pi/2)
    hp1 = scatter!(ax, Vb_Base ,markersize=10,color=:darkgreen)
    hp2 = scatter!(ax, [VBerry1_Edge;VBerry2_Edge;VBerry3_Edge] ,markersize=10,color=:blue)
    #hp3 = Comodo.meshplot!(ax, F_top, V_top, color=:red)
    hp4 = Comodo.meshplot!(ax, F_bottom, V_bottom, color=:green)
    hp6 = scatter!(ax, VHole ,markersize=10,color=:purple)
    dscreen = display(GLMakie.Screen(), fig)

   # To do from here:
   # - position and track berries
   # - Put hole in base for a string and track
   # - copy base to negative thickness 
   # - region tri mesh top and bottom
   # - refind boundary edges 
   # - create side faces
   # - combine all
   # - convert to stl
 
end