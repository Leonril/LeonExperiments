# demo_stl.jl
#= Description: This code recreates the work of my Masters project "Design
of a Modelling Framwork for the Perfomance Evaluation of Soft Robotic
Bending Actuators". It allows the parametric creation of elastomeric 
soft pneunet actuators with various design and loading parameters. 
Contact is modelled and adjustable, which is essential for pnuenet simulation.
The code also allows for the creation of simple rigid bodies for the pneunet
to interact with. =#

# Based on Revision I of my Masters Matlab code available at:https://github.com/Leonril/SoRoPneu/

## Keywords
#
# * febio_spec version 3.0
# * febio, FEBio
# * pressure loading
# * hexahedral elements, hex8
# * quad elements, quad4
# * hexahedral elements, hex8
# * pneunet actuator
# * soft robotic
# * static, solid
# * hyperelastic, Ogden
# * neo-Hookean
# * multiple steps
# * displacement logfile
# * stress logfile
##

using Comodo
using GLMakie
using GeometryBasics
using MeshIO
using FileIO
using FilePathsBase
using LinearAlgebra
using Statistics


GLMakie.closeall() #Close any open Makie windows

## Inputs:
## Plot settings
fontSize=20;
faceAlpha1=0.8;
markerSize=40;
markerSize2=20;
lineWidth=3;



## Control parameters

#Mesh Tool
MeshCreationTool=2;#
# == 1: Uses 'must points' at key variable locations
# == 2: Used designated 'mesh size' value to add more homogentity to the
# mesh - this should reduce some convergence issues

#SLL method
SLL_control=3;
# ==1: Adds a second material layer on the bottom of the pneunet of
# specified thickness
# ==2: Adds a thin shell element the base of the structure
# ==3: Adds thin shell SLL elements at in between the top and bottom layers
# created in method one. Bottom layer is now the sane material as the top
# and material 2 is the shell.

#Object type
object_control=1;
# ==1: Pneunet fixed in Z axis at time=1. Least computationally expensive
# ==2: Rigid body. At a specificed Z coordinate or plane, a ridid face is added to restrict bending and simulate a rigid object.
# ==3: Soft object: A 3D soft body is modelled. Most computationally expensive.
softObject_shape=1;
# ==1: Half of a sphere is modelled at a specified height.
# ==2: Half of a cube is modelled, at specified height. 

#Chamber Wall Contact
contact_control=2;
# ==1: Contact modelling inactive
# ==2: Contact modelling active

#Plots
plot_control=1;
# ==1: Plots on
# ==2: Plots off
plot_vector=[];# list of plots to be active regardless of plot control - see each plot's code for identifiers

#Fabrication file outputs
stl_control=0;
# ==1: Body .stls (for full print of body and/or SLL) - dependant on chosen
# SLL method
# ==2: Moulds (for moulding pneunet, as per SLL method 3)


## File saving locations
# Path names

#= UPDATE
defaultFolder = fileparts(fileparts(mfilename('fullpath')));
savePath=fullfile(defaultFolder,'data','temp');
=#


# Defining file names    UPDATE

#=
febioFebFileNamePart='tempModel';
febioFebFileName=fullfile(savePath,[febioFebFileNamePart,'.feb']); #FEB file name
febioLogFileName=fullfile(savePath,[febioFebFileNamePart,'.txt']); #FEBio log file name
febioLogFileName_disp=[febioFebFileNamePart,'_disp_out.txt']; #Log file name for exporting displacement
febioLogFileName_stress=[febioFebFileNamePart,'_stress_out.txt']; #Log file name for exporting stress
febioLogFileName_force=[febioFebFileNamePart,'_force_out.txt']; #Log file name for exporting force

STL_Pnuenet_Body=fullfile(savePath,[febioFebFileNamePart,'_PneunetBody.stl']); #Log file name for exporting .stl
STL_Pnuenet_SLL=fullfile(savePath,[febioFebFileNamePart,'_PneunetSLL.stl']); #Log file name for exporting .stl

STL_TopAMould=fullfile(savePath,[febioFebFileNamePart,'_TopAMould.stl']); #Log file name for exporting .stl
STL_TopBMould=fullfile(savePath,[febioFebFileNamePart,'_TopBMould.stl']); #Log file name for exporting .stl
STL_BottomMould=fullfile(savePath,[febioFebFileNamePart,'_BottomMould.stl']); #Log file name for exporting .stl

=#

## FEA control settings
numTimeSteps=50; #Number of time steps desired
opt_iter=7; #Optimum number of iterations
max_refs=opt_iter*4; #Max reforms
max_ups=0; #Set to zero to use full-Newton iterations
max_retries=10; #Maximum number of retires
dtmin=(1/numTimeSteps)/100; #Minimum time step size
dtmax=(1/numTimeSteps)*2;#*4; #Maximum time step size


#= UPDATE
runMode='external';#'internal';
=#


## Contact Parameters 
contactPenalty=5;#5
laugon=0;
minaug=1;
maxaug=15;
fric_coeff=0;


## Load Inputs
#Load
designPressureAngle=0.008; #(MPa)-target pressure for required bend
designPressureForce=0.014; #0.1


## Pneunet Material Properties
#Material_Bank - takes from bank of previous materials
#materialBank=Material_Bank; UPDATE

#Material parameter set

# #Elastosil
c1=75*(10^(-3)); #Shear-modulus-like parameter
m1=2.749; #Material parameter setting degree of non-linearity
# 
# c1=1;
# m1=2;
k_factor=100; #Bulk modulus factor 
k1=c1*k_factor; #Bulk modulus

#Paper
E_material2=6.5*(10^3);
Poissons_material2=0.2;

#Backup sample material 2
c2=c1*50; #Shear-modulus-like parameter
m2=2; #Material parameter setting degree of non-linearity
k2=c2*k_factor; #Bulk modulus



## Pneunet Geometry Inputs:

n=5; #no. of chambers 
pointSpacing=2;#Only active if MeshCreationTool == 2, designated aproximate size of desired cube element length and width if possible (mm) 

# X direction 

Chamber_length=8; # internal chamber length (mm)
Chamber_wt=1.2; # wall thickness expanding (mm)
Gap_length=2.5; #(mm)

# y direction 

Side_wall_thickness=2;
Chamber_width=15; # internal chamber width (mm)
Channel_width=2; # width of the intenal channel (mm)

# Z direction 

SLL_thickness=1; #Strain Limiting Layer (mm)# for Method 2, this could be made zero and Mat1 increased
SLLThicknessShell=0.1;

Mat1_base_thickness=1; #mat 1 basse layer thickness (mm)
Channel_height=3; #height of internal channel between chambers (mm)
Channel_roof_thickness=1;# (mm)
Chamber_roof_thickness=4; # (mm)
Chamber_height=25; #internal height of the chamber (mm)

## Object controls

# Soft object controls
ObjectCentreDist=[0, 0, -30] # distance from Pneuent tip to object centre (mm)

# Rigid plane distance from Pneunet bottom (Z) (for object_control==1)
rigidLength=30
rigidWidth=20

# Soft object material properties
#Object_properties=materialBank.SoftMaterialSample1; UPDATE


# Hemisphere shape controls
sphereRadius=15; # (mm)


# Half cuboid shape controls
SoftObject_PointSpacing=1; #mesh refinement parameter
SoftObjectLength=15;
SoftObjectWidth=30;
SoftObjectHeight=15; #from midplane - so half of true body height  

## Desired Mesh Inputs
# If MeshCreationTool == 1 input desired mesh geometry manually below

if MeshCreationTool == 1
# X direction

Chamber_length_elements=4; # no of mesh elements on this length
Chamber_wt_elements=2; # no of mesh elements on this length
Gap_length_elements=3; # no of mesh elements on this length

# y direction 

Side_wall_thickness_elements=1; # no of mesh elements on this length
Chamber_width_elements=6; # no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=2; # no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

# Z direction 

SLL_thickness_elements=1; #Elements used in SLL layer if method 1 applied

Mat1_base_thickness_elements=1; #mat 1 basse layer thickness (mm)
Channel_height_elements=1; #height of internal channel between chambers (mm)
Channel_roof_thickness_elements=2;# (mm)
Chamber_roof_thickness_elements=2; # (mm)
Chamber_height_elements=8; # no of mesh elements on this length (Inclusive of channel parameters)

# Generalised mesh spacing calculations
elseif MeshCreationTool == 2
    
Chamber_length_elements=ceil(Int64, Chamber_length/pointSpacing); # no of mesh elements on this length
Chamber_wt_elements=ceil(Int64,Chamber_wt/pointSpacing); # no of mesh elements on this length
Gap_length_elements=ceil(Int64,Gap_length/pointSpacing); # no of mesh elements on this length

# y direction 

Side_wall_thickness_elements=ceil(Int64,Side_wall_thickness/pointSpacing); # no of mesh elements on this length
Chamber_width_elements=ceil(Int64,Chamber_width/pointSpacing); # no of mesh elements on this length (Inclusive of channel width)
Channel_width_elements=ceil(Int64,Channel_width/pointSpacing); # no of mesh elements on this length must even if two Chamber width is even, off if chamber width is odd

if mod(Chamber_width_elements,2) != mod(Channel_width_elements,2)
    Chamber_width_elements=Chamber_width_elements+1;
end

# Z direction 

SLL_thickness_elements=ceil(Int64,SLL_thickness/pointSpacing); #Elements used in SLL layer if method 1 applied

Mat1_base_thickness_elements=ceil(Int64,Mat1_base_thickness/pointSpacing); #mat 1 basse layer thickness (mm)
Channel_height_elements=ceil(Int64,Channel_height/pointSpacing); #height of internal channel between chambers (mm)
Channel_roof_thickness_elements=ceil(Int64,Channel_roof_thickness/pointSpacing);# (mm)
Chamber_roof_thickness_elements=ceil(Int64,Chamber_roof_thickness/pointSpacing); # (mm)
Chamber_height_elements=ceil(Int64,Chamber_height/pointSpacing); # no of mesh elements on this length (Inclusive of channel parameters)
end 

if SLL_control==2
    SLL_thickness=0
    SLL_thickness_elements=0
end

## Simplfied geometry creation
# Simple geometry is defined using mesh elements only

# X direction elements needed
Length=(n*((2*Chamber_wt_elements)+Chamber_length_elements))+((n-1)*Gap_length_elements)

# Y direction elements needed
Width=(2*Side_wall_thickness_elements)+Chamber_width_elements

# Z direction elements needed
Height=SLL_thickness_elements+Mat1_base_thickness_elements+Chamber_height_elements+Chamber_roof_thickness_elements

# Creating hexahedral box 
boxDim=[Float64(Length); Float64(Width); Float64(Height)]
boxEl=[Length; Width; Height]

E_bar,V_bar,F_bar,Fb_bar,Cb_bar = hexbox(boxDim,boxEl)

# Adjusting alignment
xmin = minimum(p -> p[1], V_bar)
zmin = minimum(p -> p[3], V_bar)

V_bar = [Point(p[1] - xmin, p[2], p[3]-zmin) for p in V_bar] # undoes centre alignment on x and z axes

# Finding element centres
VE_bar = [sum(V_bar[i] for i in elem) / length(elem) for elem in E_bar]# creats a list of element centres
VE_bar = [Point(abs(p[1]), abs(p[2]), abs(p[3])) for p in VE_bar] # Abs values for easier control

#X axis logic - Inner & Outer 
LowerLimChannel=Chamber_wt_elements;
UpperLimChannel=((Length)-Chamber_wt_elements);

CX = [p[1] for p in VE_bar]
CXGap_length = zeros(length(VE_bar))
CXChamber = zeros(length(VE_bar))
CXChannel = zeros(length(VE_bar))


for i=1:1:length(VE_bar)
    
    for j=1:1:n
        #Setting Changing Limits
        LowerLimChamber=Chamber_wt_elements + ((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements))
        UpperLimChamber=Chamber_wt_elements+Chamber_length_elements  + ((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements))
        
        LowerLimGap_length=((2*Chamber_wt_elements)+Chamber_length_elements)+((j-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements))
        UpperLimGap_length=((j)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements))
        
        #Outer geometry -  if needs removal -> mark as 1
        if CX[i] > LowerLimGap_length && CX[i] < UpperLimGap_length
            CXGap_length[i]=1
        end
        
        if CX[i] > LowerLimChamber && CX[i] < UpperLimChamber
           CXChamber[i]=1
        end
        
        if CX[i] > LowerLimChannel && CX[i] < UpperLimChannel
            CXChannel[i]=1
        end
        
    end
    
end

logicDeleteOuterX = CXGap_length .!= 0
logicDeleteInnerX1 = CXChamber .!= 0
logicDeleteInnerX2 = CXChannel .!= 0


# Y and Z logic -  Inner and Outer _ simpler as not effected by no. of chambers

logicDeleteOuterZ = [p[3] > (SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements+Channel_roof_thickness_elements) for p in VE_bar]
logicKeepOuter = .!(logicDeleteOuterX .& logicDeleteOuterZ)

logicDeleteInnerY1 = [!(p[2] > Chamber_width_elements / 2) for p in VE_bar]
logicDeleteInnerY2 = [!(p[2] > Channel_width_elements / 2) for p in VE_bar]

logicDeleteInnerZ1 = [
    (p[3] > (SLL_thickness_elements+Mat1_base_thickness_elements)) & (p[3] < (Height-Chamber_roof_thickness_elements))
    for p in VE_bar
]

logicDeleteInnerZ2 = [
    (p[3] > (SLL_thickness_elements+Mat1_base_thickness_elements)) & (p[3] < (SLL_thickness_elements+Mat1_base_thickness_elements+Channel_height_elements))
    for p in VE_bar
]

# uses logic vectors to choose which elements to keep
logicKeepChamber = .!(logicDeleteInnerX1 .& logicDeleteInnerY1 .& logicDeleteInnerZ1)

logicKeepChannel = .!(logicDeleteInnerX2 .& logicDeleteInnerY2 .& logicDeleteInnerZ2)

logicKeepInner = logicKeepChamber .& logicKeepChannel
logicKeepInner = logicKeepInner[logicKeepOuter]


E1 = E_bar[logicKeepOuter] # removes elements between chambers
F1 = element2faces(E1) # gets faces of remaining elements


function boundary_face_indices(F) # Function to find boundary faces of a set of faces
    counts = Dict{Tuple{Vararg{Int}} , Int}()

    # Count occurrences of each face (order-independent)
    for f in F
        key = Tuple(sort(collect(f)))   # canonical representation
        counts[key] = get(counts, key, 0) + 1
    end

    # Collect indices of faces that occur exactly once
    idx = Int[]
    for (i, f) in enumerate(F)
        if counts[Tuple(sort(collect(f)))] == 1
            push!(idx, i)
        end
    end

    return idx
end

# Using boundary_face_indices to find outer faces 
indBoundary1 =  boundary_face_indices(F1) # Finds the boundary faces of the elements to find the outer surfaces

# Using boundary_face_indices to find outer faces and cavity faces
F2 = element2faces(E1[logicKeepInner]) # removes internal elements to create inner shape
indBoundary2 = boundary_face_indices(F2) # Finds the boundary faces of the elements to find the outer surfaces and inner surfaces
Fb = F2[indBoundary2] # All boundary faces of the final shape

# Assigning face colors based on original boundary faces
Cb = fill(7, length(Fb))  # Sample C vector for face labelling

sameface(f1, f2) = sort(collect(f1)) == sort(collect(f2)) # Function to check if two faces are the same

for q in 1:6 # for each color in the original Cb_bar
    Fq = Fb_bar[Cb_bar .== q]   # faces with color q
    for (i, f) in enumerate(Fb)
        for fq in Fq
            if sameface(f, fq)
                Cb[i] = q
                break # saves time by exiting after first match
            end
        end
    end
end


# Inline canonical function
canonical(f) = Tuple(sort(collect(f)))

# Create the set of faces to keep
keep_set = Set(canonical.(F1[indBoundary1]))

# Use set to assign inner faces to color 0
for i in eachindex(Fb)
    if !(canonical(Fb[i]) in keep_set)
        Cb[i] = 0
    end
end

# Remove unused nodes and clean up index matrices
E,V,indFix2=remove_unused_vertices(E1[logicKeepInner],V_bar);

# Function to remap faces after vertex removal
function remap_faces(F::Vector{QuadFace{Int}}, indFix::Vector{Int})
    return [QuadFace(indFix[f[1]], indFix[f[2]], indFix[f[3]], indFix[f[4]]) for f in F]
end

# Using remap_faces to update Fb and F for new V
Fb = remap_faces(Fb, indFix2)
F = remap_faces(F2, indFix2)

# Plotting vertex based colour data
Fp,Vp = separate_vertices(Fb,V) # Give each face its own point set 
Cp = simplex2vertexdata(Fp,Cb) # Convert face color data to vertex color data 

# Plotting 
plot_number=1; # Plot of the unscaled mesh
if plot_control == 1 || plot_number in plot_vector

    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Unscaled simplified geometry", azimuth=-pi/4, elevation=pi/4)
    hp1 = meshplot!(ax1, Fp, Vp, strokewidth=1.0, color=Cp, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))    
    hp2 = scatter!(ax1, V ,markersize=markerSize/4,color=:black)
    Colorbar(fig[1, 1][1, 2], hp1)
    screen = display(GLMakie.Screen(), fig)

end

# Creating an unscaled V for easy reference later
ymin = minimum(p -> p[2], V_bar)
V_element = [Point(p[1], p[2]-ymin, p[3]) for p in V]#undoing centre alignment in Z axis


## Defining the boundary conditions
# The visualization of the model boundary shows colours. These labels can be used to define boundary conditions. 

#Define supported node sets
bcSupportList= unique(vcat([collect(f) for f in Fb[Cb.==6]]...)) #Node set part of selected face (face =1 in MATLAB)

#Get pressure faces
F_pressure=Fb[Cb.==0]

#Get end face
F_end=Fb[Cb.==5]# face = 2 in MATLAB

bcTipList=unique(vcat([collect(f) for f in F_end]...)) #Unique vertexes on end face

#Get end face corner nodes
Spine_list=unique(vcat([collect(f) for f in Fb[Cb.==1]]...))#List of nodes along bottom surface of Pneunet

indForceNodes = intersect(bcTipList,Spine_list);#index list of Nodes where reaction force will be investigated


#=

## Defining Contact surfaces - UPDATE LATER  - Continued beneath
logicContactSurf=Cb==7;#chamber wall faces
Normals=patchNormal(Fb,V);#gets normal vector of all facets

logicContactPosX=Normals(:,1)==1;#where norm v is pos x 
logicContactNegX=Normals(:,1)==-1;#where norm v is neg x 

logicContactAllPrimarySets=logical(logicContactSurf.*logicContactPosX);
logicContactAllSecondarySets=logical(logicContactSurf.*logicContactNegX);

# F_contactPrimary=Fb(logicContactPrimary,:);
# F_contactSecondary=Fb(logicContactSecondary,:);


Vm=patchCentre(Fb,V);# centre location of faces, where used here should not affect x axis as contact faces are parrallel, this is used to find the faces on the x ccordinate of interest
ChamberAngleLeft={1 n-1};#empty for analysis of bending angle
ChamberAngleRight={1 n-2};



if n>1
for q=1:1:n-1
#where the primary and secondary for each contact pair should be found
desiredPrimary=(2*Chamber_wt_elements)+Chamber_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));
desiredSecondary=(2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements+((q-1)*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements));

logicContactCoordinatePrimary=Vm(:,1)==desiredPrimary;
logicContactCoordinateSecondary=Vm(:,1)==desiredSecondary;

logicContactPrimary=logical(logicContactCoordinatePrimary.*logicContactAllPrimarySets);
logicContactSecondary=logical(logicContactCoordinateSecondary.*logicContactAllSecondarySets);


if q<n
ChamberAngleLeft{q}=F(logicContactSecondary,:);
end
if q>1
ChamberAngleRight{q-1}=F(logicContactPrimary,:);
end

ContactPair.Primary{q}=Fb(logicContactPrimary,:);
ContactPair.Secondary{q}=Fb(logicContactSecondary,:);

end

#display contact pairs
plot_number=2;
if plot_control ==1 || ismember(plot_number,plot_vector)==1
cFigure;
hold on;
title('Contact Pair Surfaces')
gpatch(Fb,V,'bw','k',0.25);
for q=1:1:n-1
    gpatch(ContactPair.Primary{q},V,'rw','k');
    gpatch(ContactPair.Secondary{q},V,'y','k');
end
axisGeom;
end


end

=#

## SLL Layer
if SLL_control == 1 # Thick SLL material

    VE = [sum(V[i] for i in elem) / length(elem) for elem in E]# creats a list of element centres
    logicSLL= [p[3] < SLL_thickness_elements for p in VE]
    E1=E[.!logicSLL]#main body
    E2=E[logicSLL]#SLL layer
    E=[E1;E2]


elseif SLL_control == 2 # Shell SLL material

    FSLL=Fb[Cb.==1]# Spine faces 

    E1=E
    E2=FSLL

    E=(E1, E2) # E1 is the hex main body, E2 is the SLL quads

    # Plotting
    plot_number=3
    if plot_control == 1 || plot_number in plot_vector

        fig = Figure(size=(1200,1000))
        ax1 = AxisGeom(fig[1, 1], title="SLL shell elements", azimuth=-pi/4, elevation=pi/4)
        hp1 = meshplot!(ax1, Fp, Vp, strokewidth=1.0, color=Cp, alpha = 0.2, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))    
        hp2 = meshplot!(ax1, E[2], V, strokewidth=1.0, color=:green)
        Colorbar(fig[1, 1][1, 2], hp1)
        screen = display(GLMakie.Screen(), fig)

    end



elseif SLL_control == 3 
    
    FE =  [sum(V[i] for i in elem) / length(elem) for elem in F]
    normals_internal =  facenormal(F,V)
    keepSLL = [normals_internal[i][3] == 1 && FE[i][3] == SLL_thickness_elements for i in eachindex(F)]
    FSLL = F[keepSLL]

    E1 = E
    E2 = FSLL
    E = (E1, E2)

    # Plotting
    plot_number=4
    if plot_control == 1 || plot_number in plot_vector

        fig = Figure(size=(1200,1000))
        ax1 = AxisGeom(fig[1, 1], title="SLL shell elements", azimuth=-pi/4, elevation=pi/4)
        hp2 = meshplot!(ax1, E[2], V, strokewidth=1.0, color=:green)
        hp1 = meshplot!(ax1, Fp, Vp, strokewidth=1.0, color=Cp, alpha = 0.2, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))    
        Colorbar(fig[1, 1][1, 2], hp1)
        screen = display(GLMakie.Screen(), fig)

    end
end


## Visualizing boundary conditions. Markers plotted on the semi-transparent
# model denote the nodes in the various boundary condition lists. 

plot_number=5
if plot_control == 1 || plot_number in plot_vector

    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Boundary conditions", azimuth=-pi/4, elevation=pi/4)
    hp2 = scatter!(ax1, V[bcSupportList] ,markersize=markerSize/4,color=:black)    
    hp3 = meshplot!(ax1, F_pressure, V, strokewidth=1.0, color=:red)
    hp4 = scatter!(ax1, V[indForceNodes] ,markersize=markerSize/4,color=:red)    
    hp1 = meshplot!(ax1, Fp, Vp, strokewidth=1.0, color=:green, alpha = 0.2, )  
    Legend(fig[1, 2], [hp2, hp3, hp4] , ["BC full support", "Pressure surfaces", "End points" ])
    screen = display(GLMakie.Screen(), fig)

end


## Scaling to match desired geometry
#Creating vectors of needed Coordinates
V_height=zeros(1,Height)
V_length=zeros(1,Length)
V_width=zeros(1,Width)

#Height scaling vector creation
current_height=0
for i in 1:Height
    
    if i <= SLL_thickness_elements
        V_height[i]=current_height+(SLL_thickness*(1/SLL_thickness_elements))
        
    elseif i <= Mat1_base_thickness_elements+SLL_thickness_elements
        V_height[i]=current_height+(Mat1_base_thickness*(1/Mat1_base_thickness_elements))
        
    elseif i <= Channel_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height[i]=current_height+(Channel_height*(1/Channel_height_elements))
        
    elseif i <= Channel_roof_thickness_elements+Channel_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements 
        V_height[i]=current_height+(Channel_roof_thickness*(1/Channel_roof_thickness_elements))
        
    elseif i <= Chamber_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height[i]=current_height+((Chamber_height-Channel_height-Channel_roof_thickness)*(1/(Chamber_height_elements-Channel_height_elements-Channel_roof_thickness_elements)))
        
    elseif  i <= Chamber_roof_thickness_elements+Chamber_height_elements+Mat1_base_thickness_elements+SLL_thickness_elements
        V_height[i]=current_height+(Chamber_roof_thickness*(1/Chamber_roof_thickness_elements))
        
    end
    global current_height=V_height[i]
    
end
V_height=[0 V_height]

# Width scaling vector creation
V = [Point(p[1], p[2]-ymin, p[3]) for p in V] # undoes centre alignment on y axis

current_width=0;
for i in 1:Width
    
    if i <= Side_wall_thickness_elements
        V_width[i]=current_width+(Side_wall_thickness*(1/Side_wall_thickness_elements));
    
    elseif i <= Side_wall_thickness_elements+((Chamber_width_elements/2)-(Channel_width_elements/2))
        V_width[i]=current_width+((Chamber_width-Channel_width)*(1/(Chamber_width_elements-Channel_width_elements)));
        
    elseif i <= Side_wall_thickness_elements+((Chamber_width_elements/2)+(Channel_width_elements/2))
        V_width[i]=current_width+((Channel_width)*(1/(Channel_width_elements)));
        
    elseif i <= Side_wall_thickness_elements+Chamber_width_elements
        V_width[i]=current_width+((Chamber_width-Channel_width)*(1/(Chamber_width_elements-Channel_width_elements)));
    
    elseif i<= (2*Side_wall_thickness_elements)+Chamber_width_elements
        V_width[i]=current_width+(Side_wall_thickness*(1/Side_wall_thickness_elements));
    end
    global current_width=V_width[i];
    
end
V_width=[0 V_width]

#Length scaling vector creation
current_length=0;
for i in 1:Length
    itemp=i;
    total_chamber_length_elements=((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements)
    chamber_count=floor(itemp/total_chamber_length_elements)
    itemp=itemp-(chamber_count*((2*Chamber_wt_elements)+Chamber_length_elements+Gap_length_elements))
    
    
        if itemp == 0
            V_length[i]=current_length+(Gap_length*(1/Gap_length_elements));
            
        elseif itemp <= Chamber_wt_elements
            V_length[i]=current_length+(Chamber_wt*(1/Chamber_wt_elements));

        elseif itemp <= Chamber_wt_elements+Chamber_length_elements
            V_length[i]=current_length+(Chamber_length*(1/Chamber_length_elements));

        elseif itemp <= (2*Chamber_wt_elements)+Chamber_length_elements
            V_length[i]=current_length+(Chamber_wt*(1/Chamber_wt_elements));

        elseif itemp <= Gap_length_elements+(2*Chamber_wt_elements)+Chamber_length_elements
            V_length[i]=current_length+(Gap_length*(1/Gap_length_elements));

        end
      global current_length=V_length[i];
end
V_length=[0 V_length]

# Scaling V
V = [Point(V_length[Int(p[1]) + 1], V_width[Int(p[2]) + 1], V_height[Int(p[3]) + 1],) for p in V]

# Plotting 
Fp,Vp = separate_vertices(Fb,V) # Give each face its own point set 
Cp = simplex2vertexdata(Fp,Cb) # Convert face color data to vertex color data 

plot_number=6 # Plot of the unscaled mesh
if plot_control == 1 || plot_number in plot_vector

    fig = Figure(size=(1200,1000))
    ax1 = AxisGeom(fig[1, 1], title="Scaled geometry", azimuth=-pi/4, elevation=pi/4)
    hp1 = meshplot!(ax1, Fp, Vp, strokewidth=1.0, color=Cp, colormap=Makie.Categorical(Makie.Reverse(:Spectral)))    
    hp2 = scatter!(ax1, V ,markersize=markerSize/4,color=:black)
    Colorbar(fig[1, 1][1, 2], hp1)
    screen = display(GLMakie.Screen(), fig)

end
#=


## Object
if object_control==2
    
    # Drawing the rigid body
    V_rigid=[-rigidLength/2 rigidWidth/2; -rigidLength/2 -rigidWidth/2; rigidLength/2 -rigidWidth/2; rigidLength/2 rigidWidth/2]; #unshifted shape of rigid face
    
    regionCellRigid={V_rigid}; #A region between V1 and V2 (V2 forms a hole inside V1)
    plotOnRigid=0; #This turns on/off plotting
    pointSpacingRigid=1; #Desired point spacing
    resampleCurveOpt=1; #Option to turn on/off resampling of input boundary curves

    [F_rigid,V_rigid]=regionTriMesh2D(regionCellRigid,pointSpacingRigid,resampleCurveOpt,plotOnRigid); #creates tri mesh of rigid body
    
    #shifting position of the rigid body
    V_rigidZ=zeros(size(V_rigid,1),1); #Vector of Z heights of rigid plate
    V_rigid=[V_rigid(:,1) V_rigid(:,2) V_rigidZ]; # Adding the Z coords to the plate
    
    V_rigid(:,1)=V_rigid(:,1)+max(V(:,1))++ObjectCentreDist(1);# X direction shift
    V_rigid(:,2)=V_rigid(:,2)+(max(V(:,2))/2)+ObjectCentreDist(2);# Y direction shift - centre of Y
    V_rigid(:,3)=V_rigid(:,3)+min(V(:,3))+ObjectCentreDist(3);# Z direction shift
    
    # Merging F,V,
    F_rigid=F_rigid+size(V,1);
    Cb_rigid=ones(size(F_rigid,1),1)*(max(Cb)+1); #adding a boundary colour for plotting
    
    V=[V; V_rigid]; #updated V to include rigid body
    center_of_mass_rigid=mean(V_rigid);
    
    #contact surfaces
    contactObjectPrimary=Fb(Cb==2 | Cb==5,:); # object contact faces of the Pneunet
    contactObjectSecondary=F_rigid;
    
    
    cFigure;
    gpatch(Fb,V,Cb);
    hold on;
    gpatch(F_rigid, V, Cb_rigid)
    axisGeom; icolorbar;
    
    
elseif object_control == 3
    
    #Sphere
    if softObject_shape == 1 
        
        SoftObjectStruct.sphereRadius=sphereRadius;
        SoftObjectStruct.coreRadius=SoftObjectStruct.sphereRadius/2;
        SoftObjectStruct.numElementsMantel=3;
        SoftObjectStruct.numElementsCore=SoftObjectStruct.numElementsMantel*2; 
        
        [SoftObjectStruct]=hexMeshHemiSphere(SoftObjectStruct);
        
        Eso=SoftObjectStruct.E;
        Vso=SoftObjectStruct.V;
        Fso=SoftObjectStruct.F;
        Fbso=SoftObjectStruct.Fb;
        Cbso=SoftObjectStruct.faceBoundaryMarker;
        
    end
    
    # Cuboid
    if softObject_shape == 2 
        
        boxDim=[SoftObjectLength SoftObjectWidth SoftObjectHeight];
        boxEl=[ceil(SoftObjectLength/SoftObject_PointSpacing) ceil(SoftObjectWidth/SoftObject_PointSpacing) ceil(SoftObjectHeight/SoftObject_PointSpacing)];

        [SoftObjectStruct]=hexMeshBox(boxDim,boxEl);
        
        Eso=SoftObjectStruct.E;
        Vso=SoftObjectStruct.V;
        Fso=SoftObjectStruct.F;
        Fbso=SoftObjectStruct.Fb;
        Cbso=SoftObjectStruct.faceBoundaryMarker;
        
    end
   
    # Shifting the location of the Soft Object
    Vso(:,1)=Vso(:,1)+max(V(:,1))+ObjectCentreDist(1);# X direction shift
    Vso(:,2)=Vso(:,2)+(max(V(:,2))/2)+ObjectCentreDist(2);# Y direction shift - centre of Y
    Vso(:,3)=Vso(:,3)+min(V(:,3))+ObjectCentreDist(3);# Z direction shift
    
    cFigure;
    gpatch(Fb,V,'rw');
    hold on;
    gpatch(Fbso,Vso,Cbso)
    axisGeom;
    icolorbar;
    
    # Taking indices of object boundary conditions
    if softObject_shape ==1 #sphere
        
        bcMidPlane=unique(Fbso(Cbso==2,:))+size(V,1);
        contactObjectSecondary=Fbso(Cbso==1,:)+size(V,1); #outer contact face of hemisphere
    
    elseif softObject_shape == 2 #cuboid
        
        bcMidPlane=unique(Fbso(Cbso==5))+size(V,1);
        contactObjectSecondary=Fbso(Cbso==6 | Cbso==1 | Cbso==2,:)+size(V,1); #outer contact faces of cuboid
        
    end
    
    contactObjectPrimary=Fb(Cb==2 | Cb==5,:); # object contact faces of the Pneunet
    
    # Merging V, E, F matrices
    
    Fbso2=Fbso+size(V,1);
    E3=Eso+size(V,1);
    V=[V;Vso];
    
    cFigure;
    gpatch(Fbso2,V,'bw');
    hold on;
    gpatch(Fb,V,'rw');
    axisGeom;
end

## .stl file creation
# Pneunet body .stl files
if stl_control==1# if 3D printing .stls are needed
    [F_PneunetBody,~,~]=element2patch(E1,[],'hex8');
    Fb_PneunetBody=F_PneunetBody(tesBoundary(F_PneunetBody),:);#getting boundary of Pneunet(no SLL)
    Fb_PneunetBodyTri=quad2tri(Fb_PneunetBody,V);#changing to tri for .stl
    [Fb_PneunetBodyTri,V_stl_temp]=patchCleanUnused(Fb_PneunetBodyTri,V);#removing unused V to stop warning in command window
    
    TR_PnuenetBody=triangulation(Fb_PneunetBodyTri,V_stl_temp);
    stlwrite(TR_PnuenetBody,STL_Pnuenet_Body);#writing .stl file
    
   if SLL_control==1 # create SLL .stl if required
       [F_PneunetSLL,~,~]=element2patch(E2,[],'hex8');
        Fb_PneunetSLL=F_PneunetBody(tesBoundary(F_PneunetSLL),:);#getting boundary of SLL
        Fb_PneunetSLLTri=quad2tri(Fb_PneunetSLL,V);#changing to tri for .stl
        [Fb_PneunetSLLTri,V_stl_temp]=patchCleanUnused(Fb_PneunetSLLTri,V);#removing unused V to stop warning in command window
    
        TR_PnuenetSLL=triangulation(Fb_PneunetSLLTri,V_stl_temp);
        stlwrite(TR_PnuenetSLL,STL_Pnuenet_SLL);#writing .stl file
       
       
   end

# Mould .stl files
elseif stl_control==2
    
### Bottom mould part ###

#sample model to cut into
bottomDimEl=[3 3 2];
stlThickness=3;
[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_STLBase=meshStruct.E;
V_STLBase=meshStruct.V;
F_STLBase=meshStruct.F;
Fb_STLBase=meshStruct.Fb;
Cb_STLBase=meshStruct.faceBoundaryMarker;

#Getting centre for easier control of elements
VE_STLBase=patchCentre(E_STLBase,V_STLBase);

#Shifting data for easier deletion of sink elements
VE_STLBase(:,1)=abs(VE_STLBase(:,1));
VE_STLBase(:,2)=abs(VE_STLBase(:,2));
VE_STLBase(:,3)=VE_STLBase(:,3)-min(VE_STLBase(:,3))+0.5;

#Deleting sink elements for mould pouring
logicDeleteSTLBase= VE_STLBase(:,3)>1 & VE_STLBase(:,2) ==0 & VE_STLBase(:,1) == 0;
E_bottomMould=E_STLBase(~logicDeleteSTLBase,:);
F_bottomMould=element2patch(E_bottomMould,V_STLBase);
Fb_bottomMould=F_bottomMould(tesBoundary(F_bottomMould),:);

#scaling factors for the bottom mould
bottomMouldDimX=[0 stlThickness stlThickness+max(V(:,1)) 2*stlThickness+max(V(:,1))];
bottomMouldDimY=[0 stlThickness stlThickness+max(V(:,2)) 2*stlThickness+max(V(:,2))];
bottomMouldDimZ=[0 stlThickness stlThickness+SLL_thickness+Mat1_base_thickness];

#Moving axes
V_STLBase(:,1)=V_STLBase(:,1)-min(V_STLBase(:,1));
V_STLBase(:,2)=V_STLBase(:,2)-min(V_STLBase(:,2));
V_STLBase(:,3)=V_STLBase(:,3)-min(V_STLBase(:,3));

#scaling
for pp=1:1:size(V_STLBase,1)
    V_STLBase(pp,1)=bottomMouldDimX(V_STLBase(pp,1)+1);
    V_STLBase(pp,2)=bottomMouldDimY(V_STLBase(pp,2)+1);
    V_STLBase(pp,3)=bottomMouldDimZ(V_STLBase(pp,3)+1);
end

#plotting
cFigure;
gpatch(Fb_bottomMould,V_STLBase,'rw')
axisGeom;


### Top mould part A ###

#sample model to cut into
bottomDimEl=[(4+n+n-1) 7 3];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_stlTopA=meshStruct.E;
V_stlTopA=meshStruct.V;
F_stlTopA=meshStruct.F;
Cb_stlTopA=meshStruct.faceBoundaryMarker;

#Getting centre for easier control of elements
VE_stlTopA=patchCentre(E_stlTopA,V_stlTopA);

#Shifting data for easier deletion of sink elements
V_stlTopA(:,1)=V_stlTopA(:,1)-min(V_stlTopA(:,1));
VE_stlTopA_temp=VE_stlTopA(:,1)-min(VE_stlTopA(:,1))+0.5;
VE_stlTopA(:,1)=abs(VE_stlTopA(:,1)); # only direction where changes are symmetric
VE_stlTopA(:,2)=abs(VE_stlTopA(:,2)); # only direction where changes are symmetric

VE_stlTopA(:,3)=VE_stlTopA(:,3)-min(VE_stlTopA(:,3))+0.5;

#Deleting general area  elements for mould pouring
logicDeletestlTopAXY1=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) ==2 & VE_stlTopA(:,1)<=max(VE_stlTopA(:,1))-1;
logicDeletestlTopAXY2=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) <=2 & VE_stlTopA(:,1)<=max(VE_stlTopA(:,1))-1 & VE_stlTopA(:,1)>max(VE_stlTopA(:,1))-2;
logicKeepstlTopASink=logical(~logicDeletestlTopAXY1.*~logicDeletestlTopAXY2);#general sink shape

#Deleting the elements required for channel
VE_stlTopA(:,1)=VE_stlTopA_temp;

DeleteElTopAXCoord=3.5:2:(max(V_stlTopA(:,1))-3.5); #element centres (in x) which must be deleted to make general shape

logicDeletestlTopChannelX=ismember(VE_stlTopA(:,1),DeleteElTopAXCoord);
logicDeletestlTopChannelYZ=VE_stlTopA(:,3)>1 & VE_stlTopA(:,2) <=2 & VE_stlTopA(:,2) > 0;
logicDeletestlTopChannelYZ2=VE_stlTopA(:,3)>2 & VE_stlTopA(:,2) ==0;#second logic for vertical grooves

logicDeletestlTopChannel=logical(logicDeletestlTopChannelX.*logicDeletestlTopChannelYZ);
logicDeletestlTopChannel2=logical(logicDeletestlTopChannelX.*logicDeletestlTopChannelYZ2);
logicKeepstlTopChannel=logical(~logicDeletestlTopChannel.*~logicDeletestlTopChannel2);
logicKeepstlTopA=logical(logicKeepstlTopASink.*logicKeepstlTopChannel);#combining all logic

#contolling faces and element identification
E_TopA=E_stlTopA(logicKeepstlTopA,:);
F_TopA=element2patch(E_TopA,V_stlTopA);
Fb_stlTopA=F_TopA(tesBoundary(F_TopA),:);

#scaling factors for the bottom mould
# bottomMouldDimX=[0 stl_gen_thickness stl_gen_thickness+max(V(:,1)) 2*stl_gen_thickness+max(V(:,1))];
TopADimY=[0 stlThickness ((2*stlThickness)+Side_wall_thickness) ((2*stlThickness)+Side_wall_thickness+((Chamber_width/2)-(Channel_width/2)))...
    ((2*stlThickness)+Side_wall_thickness+((Chamber_width/2)+(Channel_width/2))) ((2*stlThickness)+(Side_wall_thickness)+((Chamber_width))) ...
    ((3*stlThickness)+(2*Side_wall_thickness)+((Chamber_width))) ((4*stlThickness)+(2*Side_wall_thickness)+((Chamber_width)))];
TopADimZEdge=[0 stlThickness stlThickness*1.5 stlThickness*2];
TopADimZCentre=[0 stlThickness stlThickness+Channel_height stlThickness+Chamber_height];

#lenght scaling factor
LengthTopA=1:1:n+n-1; #general length vector - focused on the chamber geomertry
current_length=0;
for i=1:1:n+n-1 #scaling general length vector - focused on the chamber geomertry
   
    if iseven(i)==0
       
        LengthTopA(i)=current_length+Chamber_length;
        
    elseif iseven(i)==1
        
        LengthTopA(i)=current_length+((2*Chamber_wt)+Gap_length);
        
    end
    current_length=LengthTopA(i);
end
LengthTopA=LengthTopA+(2*stlThickness)+Chamber_wt;#shifting along to account for the lips that are added next
LengthTopA=[0 stlThickness (2*stlThickness)+Chamber_wt LengthTopA (max(LengthTopA)+Chamber_wt+stlThickness) (max(LengthTopA)+Chamber_wt+(2*stlThickness))];

#Shifting V nodes
V_stlTopA(:,2)=V_stlTopA(:,2)-min(V_stlTopA(:,2));
V_stlTopA(:,3)=V_stlTopA(:,3)-min(V_stlTopA(:,3));

#applying scaling 
indEdgeTopA=find(V_stlTopA(:,2)<=1 | V_stlTopA(:,2)>=6 | V_stlTopA(:,1)<=1 | V_stlTopA(:,1)>=max(V_stlTopA(:,1)-1));

for pp=1:1:size(V_stlTopA,1)
    #  height (z)
    if ismember(pp,indEdgeTopA)==1
       
        V_stlTopA(pp,3)=TopADimZEdge(V_stlTopA(pp,3)+1);
        
    else 
        
        V_stlTopA(pp,3)=TopADimZCentre(V_stlTopA(pp,3)+1);
        
    end
    
    # x direction
    V_stlTopA(pp,1)=LengthTopA(V_stlTopA(pp,1)+1);
    
    # y direction
    
    V_stlTopA(pp,2)=TopADimY(V_stlTopA(pp,2)+1);
    
end

#plotting
cFigure;
gpatch(Fb_stlTopA,V_stlTopA,'rw')
axisGeom;

### Top Mould part B

#sample model to cut into
bottomDimEl=[(2+n+n-1) 3 2];

[meshStruct]=hexMeshBox(bottomDimEl,bottomDimEl);

E_stlTopB=meshStruct.E;
V_stlTopB=meshStruct.V;
Cb_stlTopB=meshStruct.faceBoundaryMarker;

#Getting centre for easier control of elements
VE_stlTopB=patchCentre(E_stlTopB,V_stlTopB);

#Shifting data for easier deletion of sink elements
V_stlTopB(:,1)=V_stlTopB(:,1)-min(V_stlTopB(:,1));

VE_stlTopB(:,1)=VE_stlTopB(:,1)-min(VE_stlTopB(:,1))+0.5;
VE_stlTopB(:,2)=abs(VE_stlTopB(:,2)); # 
VE_stlTopB(:,3)=VE_stlTopB(:,3)-min(VE_stlTopB(:,3))+0.5;

#Deleting general area  elements for mould pouring
# DeleteElTopBFill=1.5:2:(max(V_stlTopB(:,1))-1.5); #element centres (in x) which must be deleted to make general shape
DeleteElTopBWall=1.5:2:(max(V_stlTopB(:,1))-1.5); #element centres (in x) which must be deleted to make general shape

logicDeleteTopBTrench=VE_stlTopB(:,3)<1 & VE_stlTopB(:,2) <1 & VE_stlTopB(:,1) > 1 & VE_stlTopB(:,1) < max(V_stlTopB(:,1))-1; # creating trench
logicDeleteTopBWallX=ismember(VE_stlTopB(:,1),DeleteElTopBWall);#x direction for deleting scraper elements
logicDeleteTopBWallYZ=VE_stlTopB(:,3)>1 & VE_stlTopB(:,2) <1;#y z direction for deleting scraper elements
logicDeleteTopBWall=logical(logicDeleteTopBWallX.*logicDeleteTopBWallYZ); #combining
logicDeleteTopBWall=logical(logicDeleteTopBWall(~logicDeleteTopBTrench,:)); #accounting for other logic

#controlling faces and element identification & applying logic
E_stlTopB=E_stlTopB(~logicDeleteTopBTrench,:);
E_stlTopB=E_stlTopB(~logicDeleteTopBWall,:);

F_stlTopB=element2patch(E_stlTopB,V_stlTopB);
Fb_stlTopB=F_stlTopB(tesBoundary(F_stlTopB),:);

#scaling
TopBDimY=[0 stlThickness (stlThickness+(2*Side_wall_thickness)+Chamber_width) ((2*stlThickness)+(2*Side_wall_thickness)+Chamber_width)];
TopBDimZ=[0 (Channel_height+Channel_roof_thickness) (Chamber_height+Chamber_roof_thickness)];

LengthTopB=1:1:n+n-1;
current_length=0;
for i=1:1:n+n-1 #scaling general length vector - focused on the chamber geomertry
   
    if iseven(i)==0
       
        LengthTopB(i)=current_length+(2*Chamber_wt)+Chamber_length;
        
    elseif iseven(i)==1
        
        LengthTopB(i)=current_length+Gap_length;
        
    end
    current_length=LengthTopB(i);
end
LengthTopB=LengthTopB+stlThickness;#shifting along to account for the lips that are added next
LengthTopB=[1 stlThickness LengthTopB (max(LengthTopB)+stlThickness)-1];

#shifting
V_stlTopB(:,2)=V_stlTopB(:,2)-min(V_stlTopB(:,2));
V_stlTopB(:,3)=V_stlTopB(:,3)-min(V_stlTopB(:,3));


#assigning scaled values
for pp=1:1:size(V_stlTopB,1)
    #  height (z)
    V_stlTopB(pp,3)=TopBDimZ(V_stlTopB(pp,3)+1);
    
    # x direction
    V_stlTopB(pp,1)=LengthTopB(V_stlTopB(pp,1)+1);
    
    # y direction
    
    V_stlTopB(pp,2)=TopBDimY(V_stlTopB(pp,2)+1);
    
end

#plotting
cFigure;
gpatch(Fb_stlTopB,V_stlTopB,'rw')
axisGeom;


### .stl file creation
# Bottom
Fb_bottomMould=quad2tri(Fb_bottomMould,V_STLBase);#changing to tri for .stl
[Fb_bottomMould,V_STLBase]=patchCleanUnused(Fb_bottomMould,V_STLBase);#removing unused V to stop warning in command window
TR_TopA=triangulation(Fb_bottomMould,V_STLBase);
stlwrite(TR_TopA,STL_BottomMould);#writing .stl file

# Top A
Fb_stlTopA=quad2tri(Fb_stlTopA,V_stlTopA);#changing to tri for .stl
[Fb_stlTopA,V_stlTopA]=patchCleanUnused(Fb_stlTopA,V_stlTopA);#removing unused V to stop warning in command window
TR_TopA=triangulation(Fb_stlTopA,V_stlTopA);
stlwrite(TR_TopA,STL_TopAMould);#writing .stl file

# Top B
Fb_stlTopB=quad2tri(Fb_stlTopB,V_stlTopB);#changing to tri for .stl
[Fb_stlTopB,V_stlTopB]=patchCleanUnused(Fb_stlTopB,V_stlTopB);#removing unused V to stop warning in command window
TR_TopB=triangulation(Fb_stlTopB,V_stlTopB);
stlwrite(TR_TopB,STL_TopBMould);#writing .stl file
       
end


## Defining the FEBio input structure
# See also |febioStructTemplate| and |febioStruct2xml| and the FEBio user
# manual.

#Get a template with default settings 
[febio_spec]=febioStructTemplate;

#febio_spec version 
febio_spec.ATTR.version='3.0'; 

#Module section
febio_spec.Module.ATTR.type='solid'; 

#Create control structure for use by all steps
stepStruct.Control.analysis='STATIC';
stepStruct.Control.time_steps=numTimeSteps;
stepStruct.Control.step_size=1/numTimeSteps;
stepStruct.Control.solver.max_refs=max_refs;
stepStruct.Control.solver.max_ups=max_ups;
stepStruct.Control.solver.symmetric_stiffness=0;
stepStruct.Control.time_stepper.dtmin=dtmin;
stepStruct.Control.time_stepper.dtmax=dtmax; 
stepStruct.Control.time_stepper.max_retries=max_retries;
stepStruct.Control.time_stepper.opt_iter=opt_iter;

#Add template based default settings to proposed control section
[stepStruct.Control]=structComplete(stepStruct.Control,febio_spec.Control,1); #Complement provided with default if missing

#Remove control field (part of template) since step specific control sections are used
febio_spec=rmfield(febio_spec,'Control'); 

febio_spec.Step.step{1}.Control=stepStruct.Control;
febio_spec.Step.step{1}.ATTR.id=1; #inflate to bending pressure
febio_spec.Step.step{2}.Control=stepStruct.Control;
febio_spec.Step.step{2}.ATTR.id=2; #inflate further (used to fix tip)

#Material section
material_number=1; #counter for material to reduce loops, easier to read than indexing size

materialName1='Material1';
febio_spec.Material.material{material_number}.ATTR.name=materialName1;
febio_spec.Material.material{material_number}.ATTR.type='Ogden';
febio_spec.Material.material{material_number}.ATTR.id=material_number;
febio_spec.Material.material{material_number}.c1=c1;
febio_spec.Material.material{material_number}.m1=m1;
febio_spec.Material.material{material_number}.c2=c1;
febio_spec.Material.material{material_number}.m2=-m1;
febio_spec.Material.material{material_number}.k=k1;

if SLL_control == 1
material_number=material_number+1;

materialName2='Material2';
febio_spec.Material.material{material_number}.ATTR.name=materialName2;
febio_spec.Material.material{material_number}.ATTR.type='Ogden';
febio_spec.Material.material{material_number}.ATTR.id=material_number;
febio_spec.Material.material{material_number}.c1=c2;
febio_spec.Material.material{material_number}.m1=m2;
febio_spec.Material.material{material_number}.c2=c2;
febio_spec.Material.material{material_number}.m2=-m2;
febio_spec.Material.material{material_number}.k=k2;

elseif SLL_control >= 2
material_number=material_number+1;
    
materialName2='Material2';
febio_spec.Material.material{material_number}.ATTR.name=materialName2;
febio_spec.Material.material{material_number}.ATTR.type='neo-Hookean';
febio_spec.Material.material{material_number}.ATTR.id=2;
febio_spec.Material.material{material_number}.E=E_material2;
febio_spec.Material.material{material_number}.v=Poissons_material2;
end

if  object_control == 2
    material_number=material_number+1;
    rigidID=material_number;
    materialRigidBody='MaterialRigidBody';
    febio_spec.Material.material{material_number}.ATTR.name=materialRigidBody;
    febio_spec.Material.material{material_number}.ATTR.type='rigid body';
    febio_spec.Material.material{material_number}.ATTR.id=material_number;
    febio_spec.Material.material{material_number}.density=1;
    febio_spec.Material.material{material_number}.center_of_mass=center_of_mass_rigid;

    
elseif object_control==3 # soft object material
    material_number=material_number+1;
    
    MaterialSoftObject='Material3';
    febio_spec.Material.material{material_number}=materialBank.SoftMaterialSample1;
    febio_spec.Material.material{material_number}.ATTR.name=MaterialSoftObject;
    febio_spec.Material.material{material_number}.ATTR.id=material_number;
     
end

#Mesh section
# -> Nodes
febio_spec.Mesh.Nodes{1}.ATTR.name='nodeSet_all'; #The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id=(1:size(V,1))'; #The node id's
febio_spec.Mesh.Nodes{1}.node.VAL=V; #The nodel coordinates 

# -> Elements
element_number=1; #counter for elements, looks cleaner than more if loops or indexing based on size

partBodyPneunet='Part1';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partBodyPneunet; #Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; #Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=(1:1:size(E1,1))'; #Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E1; #The element matrix

if SLL_control == 1 #Thick SLL -> hexa8
    element_number=element_number+1;
    
partSLL='Part2';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partSLL; #Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; #Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; #Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E2; #The element matrix

elseif SLL_control >= 2 #Thin SLL -> quad4
    element_number=element_number+1;
    
partSLL='Part2';
febio_spec.Mesh.Elements{element_number}.ATTR.name=partSLL; #Name of this part
febio_spec.Mesh.Elements{element_number}.ATTR.type='quad4'; #Element type 
febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+(1:1:size(E2,1))'; #Element id's
febio_spec.Mesh.Elements{element_number}.elem.VAL=E2; #The element matrix

end

if SLL_control==0
    E2=[]; #setting empty E2 to remove need for homogenous Pnuenet loop for object control
end

if  object_control ==2
    element_number=element_number+1;
    
    partRigidBody='Part3';
    febio_spec.Mesh.Elements{element_number}.ATTR.name=partRigidBody; #Name of this part
    febio_spec.Mesh.Elements{element_number}.ATTR.type='tri3'; #Element type 
    febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(F_rigid,1))'; #Element id's
    febio_spec.Mesh.Elements{element_number}.elem.VAL=F_rigid; #The element matrix
     
    
elseif object_control ==3 # Assigning part for soft object
    element_number=element_number+1;

    partSoftObject='Part3';
    febio_spec.Mesh.Elements{element_number}.ATTR.name=partSoftObject; #Name of this part
    febio_spec.Mesh.Elements{element_number}.ATTR.type='hex8'; #Element type 
    febio_spec.Mesh.Elements{element_number}.elem.ATTR.id=size(E1,1)+size(E2,1)+(1:1:size(E3,1))'; #Element id's
    febio_spec.Mesh.Elements{element_number}.elem.VAL=E3; #The element matrix
       
end

# -> Surfaces
surfaceName1='LoadedSurface';
febio_spec.Mesh.Surface{1}.ATTR.name=surfaceName1;
febio_spec.Mesh.Surface{1}.quad4.ATTR.id=(1:1:size(F_pressure,1))';
febio_spec.Mesh.Surface{1}.quad4.VAL=F_pressure;

leaderStringSurf='contactSurface';
leaderStringContact='Contact';
if contact_control == 2 # contact surfaces
contactSurfacesValPrimary=2:2:(2*(n-1));# attribute number of contact faces
contactSurfacesValSecondary=3:2:((2*(n-1))+1);# attribute number of contact faces

if n>1 #only need face contact for more than one chamber
    for q=1:1:n-1
        numStringPrimary=num2str(contactSurfacesValPrimary(q));
        numStringSecondary=num2str(contactSurfacesValSecondary(q));
        numStringPair=num2str(q);
              
contactSurfacesStringPrimary=strcat(leaderStringSurf,numStringPrimary);#creates string 'ContactSurfaceINT'
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.ATTR.name=contactSurfacesStringPrimary;
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Primary{q},1))';#size(F_pressure,1)+
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)}.quad4.VAL=ContactPair.Primary{q};

contactSurfacesStringSecondary=strcat(leaderStringSurf,numStringSecondary);#creates string 'ContactSurfaceINT'
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.ATTR.name=contactSurfacesStringSecondary;
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.ATTR.id=(1:1:size(ContactPair.Secondary{q},1))';#size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)}.quad4.VAL=ContactPair.Secondary{q};

# -> Surface pairs
febio_spec.Mesh.SurfacePair{q}.ATTR.name=strcat(leaderStringContact,numStringPair);
febio_spec.Mesh.SurfacePair{q}.primary=contactSurfacesStringPrimary;
febio_spec.Mesh.SurfacePair{q}.secondary=contactSurfacesStringSecondary;
    end
end

else #where there is no chamber contact
q=1;#counter for contact surfaces
contactSurfacesValPrimary(q)=0;# dummy value to allow other contact using same input
contactSurfacesValSecondary(q)=1;# dummy value to allow other contact using same input

end




if object_control>=2 #if there is an object/rigid body
contactSurfacesObjectPrimary='contactObjectPrimary';
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.ATTR.name=contactSurfacesObjectPrimary;
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.quad4.ATTR.id=(1:1:size(contactObjectPrimary,1))';#size(F_pressure,1)+
febio_spec.Mesh.Surface{contactSurfacesValPrimary(q)+2}.quad4.VAL=contactObjectPrimary;

contactSurfacesObjectSecondary='contactObjectSecondary';
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.ATTR.name=contactSurfacesObjectSecondary;

if object_control==3 #if soft object -> quad4
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.quad4.ATTR.id=(1:1:size(contactObjectSecondary,1))';#size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.quad4.VAL=contactObjectSecondary;
elseif object_control==2 #if soft object -> tri3
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.tri3.ATTR.id=(1:1:size(contactObjectSecondary,1))';#size(F_pressure,1)+size(F_contactPrimary,1)+
febio_spec.Mesh.Surface{contactSurfacesValSecondary(q)+2}.tri3.VAL=contactObjectSecondary;
end

# -> Surface pairs
if contact_control==1 # shifting to match up q values if no contact
    q=0;
elseif n==1
    q=0;
end

numStringPair=num2str(q+1);#string for number of the name of pair
febio_spec.Mesh.SurfacePair{q+1}.ATTR.name=strcat(leaderStringContact,numStringPair);
febio_spec.Mesh.SurfacePair{q+1}.primary=contactSurfacesObjectPrimary;
febio_spec.Mesh.SurfacePair{q+1}.secondary=contactSurfacesObjectSecondary;
end


# -> NodeSets
nodeSetName1='bcSupportList';
febio_spec.Mesh.NodeSet{1}.ATTR.name=nodeSetName1;
febio_spec.Mesh.NodeSet{1}.node.ATTR.id=bcSupportList(:);

nodeSetName2='bcTipList';
febio_spec.Mesh.NodeSet{2}.ATTR.name=nodeSetName2;
febio_spec.Mesh.NodeSet{2}.node.ATTR.id=indForceNodes(:);

if object_control==3
nodeSetName3='bcSoftObjectList';
febio_spec.Mesh.NodeSet{3}.ATTR.name=nodeSetName3;
febio_spec.Mesh.NodeSet{3}.node.ATTR.id=bcMidPlane(:);
end

#MeshDomains section
solid_number=1;
shell_number=0;

febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partBodyPneunet;
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=materialName1;

if SLL_control == 1 #SLL domain
    solid_number=solid_number+1;
    
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partSLL;
febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=materialName2;
elseif SLL_control >= 2
    shell_number=shell_number+1;
    
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.name=partSLL;
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.mat=materialName2;
end

if object_control==2 #object domain
    shell_number=shell_number+1;
    
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.name=partRigidBody;
febio_spec.MeshDomains.ShellDomain{shell_number}.ATTR.mat=materialRigidBody;
    
    
elseif object_control==3
    solid_number=solid_number+1;

    febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.name=partSoftObject;
    febio_spec.MeshDomains.SolidDomain{solid_number}.ATTR.mat=MaterialSoftObject;
 
end

#MeshData secion
#-> Element data
if SLL_control >= 2
febio_spec.MeshData.ElementData{1}.ATTR.var='shell thickness';
febio_spec.MeshData.ElementData{1}.ATTR.elem_set=partSLL;
febio_spec.MeshData.ElementData{1}.elem.ATTR.lid=(1:1:size(E2,1))';
febio_spec.MeshData.ElementData{1}.elem.VAL=SLLThicknessShell*ones(size(E2,1),size(E2,2));
end

#Boundary condition section 
# -> Fix boundary conditions
febio_spec.Boundary.bc{1}.ATTR.type='fix';
febio_spec.Boundary.bc{1}.ATTR.node_set=nodeSetName1;
febio_spec.Boundary.bc{1}.dofs='x,y,z';

if object_control==3 #fixing the object nodes
febio_spec.Boundary.bc{2}.ATTR.type='fix';
febio_spec.Boundary.bc{2}.ATTR.node_set=nodeSetName3;
febio_spec.Boundary.bc{2}.dofs='x,y,z';
end

if object_control==1
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.type='prescribe';
febio_spec.Step.step{2}.Boundary.bc{1}.ATTR.node_set=nodeSetName2;
febio_spec.Step.step{2}.Boundary.bc{1}.dof='z';
febio_spec.Step.step{2}.Boundary.bc{1}.scale.ATTR.lc=2;
febio_spec.Step.step{2}.Boundary.bc{1}.scale.VAL=0;
febio_spec.Step.step{2}.Boundary.bc{1}.relative=1;
end

#Loads section
# -> Surface load
febio_spec.Loads.surface_load{1}.ATTR.type='pressure';
febio_spec.Loads.surface_load{1}.ATTR.surface=surfaceName1;
febio_spec.Loads.surface_load{1}.pressure.ATTR.lc=1;
febio_spec.Loads.surface_load{1}.pressure.VAL=designPressureForce;
febio_spec.Loads.surface_load{1}.symmetric_stiffness=1;


#Contact section
q=0;
if contact_control == 2
for q=1:1:n-1
febio_spec.Contact.contact{q}.ATTR.type='sliding-elastic';
febio_spec.Contact.contact{q}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q}.ATTR.name;
febio_spec.Contact.contact{q}.two_pass=1;
febio_spec.Contact.contact{q}.laugon=laugon;
febio_spec.Contact.contact{q}.tolerance=0.2;
febio_spec.Contact.contact{q}.gaptol=0;
febio_spec.Contact.contact{q}.minaug=minaug;
febio_spec.Contact.contact{q}.maxaug=maxaug;
febio_spec.Contact.contact{q}.search_tol=0.01;
febio_spec.Contact.contact{q}.search_radius=0.1*sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));
febio_spec.Contact.contact{q}.symmetric_stiffness=0;
febio_spec.Contact.contact{q}.auto_penalty=1;
febio_spec.Contact.contact{q}.penalty=contactPenalty;
febio_spec.Contact.contact{q}.fric_coeff=fric_coeff;
end

if object_control>=2
#     if object_control==2
#         febio_spec.Contact.contact{q+1}.ATTR.type='sliding-elastic';
#     else
        febio_spec.Contact.contact{q+1}.ATTR.type='sliding-elastic';
#       
#     end
febio_spec.Contact.contact{q+1}.ATTR.surface_pair=febio_spec.Mesh.SurfacePair{q+1}.ATTR.name;
febio_spec.Contact.contact{q+1}.two_pass=0;#1
febio_spec.Contact.contact{q+1}.laugon=laugon;
febio_spec.Contact.contact{q+1}.tolerance=0.2;
febio_spec.Contact.contact{q+1}.gaptol=0;
febio_spec.Contact.contact{q+1}.minaug=minaug;
febio_spec.Contact.contact{q+1}.maxaug=maxaug;
febio_spec.Contact.contact{q+1}.search_tol=0.01;
febio_spec.Contact.contact{q+1}.search_radius=sqrt(sum((max(V,[],1)-min(V,[],1)).^2,2));#0.1*
febio_spec.Contact.contact{q+1}.symmetric_stiffness=0;
febio_spec.Contact.contact{q+1}.auto_penalty=1;
febio_spec.Contact.contact{q+1}.penalty=contactPenalty*2;
febio_spec.Contact.contact{q+1}.fric_coeff=fric_coeff;

end
end

if object_control==2
#Rigid section 
# -> Prescribed rigid body boundary conditions
febio_spec.Rigid.rigid_constraint{1}.ATTR.name='RigidFix_1';
febio_spec.Rigid.rigid_constraint{1}.ATTR.type='fix';
febio_spec.Rigid.rigid_constraint{1}.rb=rigidID;
febio_spec.Rigid.rigid_constraint{1}.dofs='Rx,Ry,Rz,Ru,Rv,Rw';
end


#LoadData section
# -> load_controller
febio_spec.LoadData.load_controller{1}.ATTR.id=1;
febio_spec.LoadData.load_controller{1}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate='LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL=[0 0; 1 (designPressureAngle/designPressureForce); 2 1];

if object_control==1 # Loadcurve for the Z node fix method 
febio_spec.LoadData.load_controller{2}.ATTR.id=2;#loadcurve ID no.
febio_spec.LoadData.load_controller{2}.ATTR.type='loadcurve';
febio_spec.LoadData.load_controller{2}.interpolate='STEP'; 
febio_spec.LoadData.load_controller{2}.points.point.VAL=[1 0; 2 1];
end

#Output section 
# -> log file
febio_spec.Output.logfile.ATTR.file=febioLogFileName;
febio_spec.Output.logfile.node_data{1}.ATTR.file=febioLogFileName_disp;
febio_spec.Output.logfile.node_data{1}.ATTR.data='ux;uy;uz';
febio_spec.Output.logfile.node_data{1}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{1}.VAL=1:size(V,1);

febio_spec.Output.logfile.node_data{2}.ATTR.file=febioLogFileName_force;
febio_spec.Output.logfile.node_data{2}.ATTR.data='Rx;Ry;Rz';
febio_spec.Output.logfile.node_data{2}.ATTR.delim=',';
febio_spec.Output.logfile.node_data{2}.VAL=1:size(V,1);

febio_spec.Output.logfile.element_data{1}.ATTR.file=febioLogFileName_stress;
febio_spec.Output.logfile.element_data{1}.ATTR.data='s1';
febio_spec.Output.logfile.element_data{1}.ATTR.delim=',';
if SLL_control==1
febio_spec.Output.logfile.element_data{1}.VAL=1:size(E,1);
elseif SLL_control>=2
    if object_control==3
febio_spec.Output.logfile.element_data{1}.VAL=1:(size(E{1},1)+size(E{2},1)+size(E3,1));
    else
febio_spec.Output.logfile.element_data{1}.VAL=1:(size(E{1},1)+size(E{2},1));        
    end
end    
## Quick viewing of the FEBio input file structure
# The |febView| function can be used to view the xml structure in a MATLAB
# figure window. 

##
# |febView(febio_spec); #Viewing the febio file|

## Exporting the FEBio input file
# Exporting the febio_spec structure to an FEBio input file is done using
# the |febioStruct2xml| function. 

febioStruct2xml(febio_spec,febioFebFileName); #Exporting to file and domNode
# febView(febioFebFileName); 

## Running the FEBio analysis
# To run the analysis defined by the created FEBio input file the
# |runMonitorFEBio| function is used. The input for this function is a
# structure defining job settings e.g. the FEBio input file name. The
# optional output runFlag informs the user if the analysis was run
# succesfully. 

febioAnalysis.run_filename=febioFebFileName; #The input file name
febioAnalysis.run_logname=febioLogFileName; #The name for the log file
febioAnalysis.disp_on=1; #Display information on the command window
febioAnalysis.runMode=runMode;
febioAnalysis.maxLogCheckTime=100; #Max log file checking time - EDITED 

[runFlag]=runMonitorFEBio(febioAnalysis);#START FEBio NOW!!!!!!!!

## Import FEBio results 

if runFlag==1 #i.e. a succesful run
    
     ## 
    # Importing nodal displacements from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_disp),1,1);
    
    #Access data
    N_disp_mat=dataStruct.data; #Displacement
    timeVec=dataStruct.time; #Time
    
    #Create deformed coordinate set
    V_DEF=N_disp_mat+repmat(V,[1 1 size(N_disp_mat,3)]);
               
    ## 
    # Plotting the simulated results using |anim8| to visualize and animate
    # deformations 
    
    DN_magnitude=sqrt(sum(N_disp_mat(:,:,end).^2,2)); #Current displacement magnitude
        
    # Create basic view and store graphics handle to initiate animation
    plot_number=7;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
    hf=cFigure; #Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('Displacement magnitude [mm]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),DN_magnitude,'k',1); #Add graphics object to animate
#     hp.Marker='.';
#     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([0 max(DN_magnitude)]);    
    axis(axisLim(V_DEF)); #Set axis limits statically    
    camlight headlight;        
        
    # Set up animation features
    animStruct.Time=timeVec; #The time vector    
    for qt=1:1:size(N_disp_mat,3) #Loop over time increments        
        DN_magnitude=sqrt(sum(N_disp_mat(:,:,qt).^2,2)); #Current displacement magnitude
                
        #Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; #Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; #Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),DN_magnitude}; #Property values for to set in order to animate
    end        
    anim8(hf,animStruct); #Initiate animation feature    
    drawnow;
    end        
    ##
    # Importing element stress from a log file
    dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_stress),1,1);
    
    #Access data
    E_stress_mat=dataStruct.data;
    
    E_stress_mat(isnan(E_stress_mat))=0;
    
    ## 
    # Plotting the simulated results using |anim8| to visualize and animate
    # deformations 
     if SLL_control==1
    [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,end));
     elseif SLL_control>=2
    [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
#     CV{1}=faceToVertexMeasure(E1,V,E_stress_mat(:,:,end));
#     CV{2}=faceToVertexMeasure(E2,V,E_stress_mat(:,:,end));
    end
    
    # Create basic view and store graphics handle to initiate animation
    plot_number=8;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
    hf=cFigure; #Open figure  
    gtitle([febioFebFileNamePart,': Press play to animate']);
    title('$\sigma_{1}$ [MPa]','Interpreter','Latex')
    hp=gpatch(Fb,V_DEF(:,:,end),CV,'k',1); #Add graphics object to animate
#     hp.Marker='.';
#     hp.MarkerSize=markerSize2;
    hp.FaceColor='interp';
    
    axisGeom(gca,fontSize); 
    colormap(gjet(250)); colorbar;
    caxis([min(E_stress_mat(:)) max(E_stress_mat(:))]);    
    axis(axisLim(V_DEF)); #Set axis limits statically    
    camlight headlight;        
        
    # Set up animation features
    animStruct.Time=timeVec; #The time vector    
    for qt=1:1:size(N_disp_mat,3) #Loop over time increments        
        if SLL_control==1
        [CV]=faceToVertexMeasure(E,V,E_stress_mat(:,:,qt));
        elseif SLL_control>=2
        [CV]=faceToVertexMeasure(E1,V,E_stress_mat(:,:,qt));
        end    
        #Set entries in animation structure
        animStruct.Handles{qt}=[hp hp]; #Handles of objects to animate
        animStruct.Props{qt}={'Vertices','CData'}; #Properties of objects to animate
        animStruct.Set{qt}={V_DEF(:,:,qt),CV}; #Property values for to set in order to animate
    end        
    anim8(hf,animStruct); #Initiate animation feature    
    drawnow;
    end


## Bending Angle of Each Chamber
Theta=zeros(1,n-2);


if n>4 # results only comparable with 5+ chambers as there are no chambers far enough from an end of the Pneunet
    
    for i=1:1:n-2 
    normLeft=patchNormal(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3)));#normals of left face of chamber wall
    normRight=patchNormal(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));#normals of right face of chamber wall
    
    normLeft=mean(normLeft,1);#average normal for left side
    normRight=mean(normRight,1);#average normal for right side
    
    Theta(i)=rad2deg(acos((dot(normLeft,normRight,2))/((norm(normLeft)*(norm(normRight))))));# angle between chamber wall normals
    
    
    
    end
    
    plot_number=9;
    if plot_control==1 || ismember(plot_number,plot_vector)==1
    cFigure;
    gpatch(Fb,V_DEF(:,:,size(V_DEF,3)),Cb,'k',0.2); hold on;
    title('Contact surface normal vectors')
    for i=1:1:n-2 
        patchNormPlot(ContactPair.Primary{i+1},V_DEF(:,:,size(V_DEF,3)));
        patchNormPlot(ContactPair.Secondary{i},V_DEF(:,:,size(V_DEF,3))); 
    end
    axisGeom;
    disp(Theta)#display chamber angles to command window
    end
    
end


## Bending Angle using Tip Nodes
BendingAngle=zeros(size(V_DEF,3),3);#empty vector of bending angles over time
BendingAngle(:,3)=timeVec;#assigning time values
YZnormVec=[1 0 0];

BendingAngleRefPoint=mean(V(indForceNodes,:),1);

for i=1:1:size(V_DEF,3)
    
    normEnd=mean(patchNormal(F_end,V_DEF(:,:,i)),1);#avergage normal vector of the last face
    nodeEnd=mean(V_DEF(indForceNodes,:,i),1);#averaged end Node tracking over time
    
    BendingAngle(i,1)=rad2deg(acos((dot(normEnd,YZnormVec,2))/((norm(normEnd)*(norm(YZnormVec))))));#angle of end face normal vector using 2D dot product
    BendingAngle(i,2)=rad2deg(atan((abs(nodeEnd(3)))/(nodeEnd(1))));#saving angle of end node
    
end
##
plot_number=10;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
ha(1)=plot(BendingAngle(:,3),BendingAngle(:,1),'r');
hold on;
ha(2)=plot(BendingAngle(:,3),BendingAngle(:,2),'b');
legend(ha,{'Normal Vector Angle', 'End Point Angle'});#'SLLsurface'
xlabel('Time (s)'); ylabel('Angle (°)');
#title('Bending Angle');
end
## Force

# Importing nodal displacements from a log file
dataStruct=importFEBio_logfile(fullfile(savePath,febioLogFileName_force),1,1);

#Access data
Force_mat=dataStruct.data; #All force values
timeVec=dataStruct.time; #Time

# Extract Forces at Restricted nodes
Force_mean=zeros(1,size(Force_mat,3));
Force_total=zeros(1,size(Force_mat,3));

for i=1:1:size(Force_mat,3)
    Force_nodal=Force_mat(indForceNodes,3,i); #Force only at restricted nodes
    Force_mean(i)=mean(Force_nodal);#mean force at the end nodes due to Z restriction
    Force_total(i)=sum(Force_nodal);# sum of all the end node forces due to Z restriction
end


## Plotting
plot_number=11;
if plot_control==1 || ismember(plot_number,plot_vector)==1
cFigure;
#hforce(1)=plot(timeVec,Force_mean,'r');
hforce=plot(timeVec,Force_total,'b');
#title('Nodal Force ()');
xlabel('Time (s)'); ylabel('Total End-Node Force in z-direction (N)');
end

end # successful run check



=#