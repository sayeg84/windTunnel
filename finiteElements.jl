using Plots, Triangulate, BenchmarkTools, NearestNeighbors, LinearAlgebra, Printf
gr()

@doc """
   triangleArea(vertex)

    Calculates the are of triangle with vertices 
    `vertex`. `vertex` must be a (2,3) Float64 array, 
    with each column representing a 2D vertex.
"""
function triangleArea(vertex::Array{Float64,2})::Float64
    if size(vertex) != (2,3)
        error("Vertex size must be (2,3)")
    end
    # explicit formula using cross product.
    return 1/2*abs((vertex[1,1]-vertex[1,3])*(vertex[2,2]-vertex[2,1])-(vertex[1,1]-vertex[1,2])*(vertex[2,3]-vertex[2,1]))
end

test = [0.0 0.0 2.0 ;
0.0 2.0 0.0
]

println(triangleArea(test))
@doc """
   triangleAreaArray(triangulation::TriangulateIO)

    Calculates the areas of all triangles in a triangulation and returns it
    as a (1,ntriangles) Float64 Array
"""
function triangleAreaArray(triangulation::TriangulateIO)::Array{Float64,2}
    #npoints = size(triangulation.pointlist)[2]
    nt = size(triangulation.trianglelist)[2]
    res = Array{Float64,2}(undef,(1,nt))
    for i in 1:nt
        tlist = triangulation.trianglelist[:,i]
        res[1,i] = triangleArea(triangulation.pointlist[:,tlist])
    end
    return res
end
@doc """
    addTriangleArea!(triangulation::TriangulateIO)

    Adds the array of triangle areas to the `trianguleattributelist`
    field of `triangulation`.
"""
function addTriangleArea!(triangulation::TriangulateIO)
    areas = triangleAreaArray(triangulation)
    if size(triangulation.triangleattributelist) == (0,0)
        triangulation.triangleattributelist = areas
    else
        println("Warning: $(size(triangulation.triangleattributelist)[1]) attributes already exist")
        triangulation.triangleattributelist = vcat(triangulation.triangleattributelist,areas)
    end
    return nothing
end

@doc """
   barycenter(vertex)

    Calculates the barycenter of triangle with vertices 
    `vertex`. `vertex` must be a (2,3) Float64 array, 
    with each column representing a 2D vertex.
"""
function barycenter(vertex::Array{Float64,2})::Array{Float64,1}
    if size(vertex) != (2,3)
        error("Vertex size must be (2,3)")
    end
    # explicit formula using cross product.
    return 1/3*reshape(sum(vertex,dims=2),2)
end

@doc """
   triangleAreaArray(triangulation::TriangulateIO)

    Calculates the areas of all triangles in a triangulation and returns it
    as a (1,ntriangles) Float64 Array
"""
function barycenterArray(triangulation::TriangulateIO)::Array{Float64,2}
    #npoints = size(triangulation.pointlist)[2]
    ntri = size(triangulation.trianglelist)[2]
    res = Array{Float64,2}(undef,(2,ntri))
    for i in 1:ntri
        tlist = triangulation.trianglelist[:,i]
        res[:,i] = barycenter(triangulation.pointlist[:,tlist])
    end
    return res
end


@doc """
    plotTriangle(triangulation::TriangulateIO;edgColor="black",verColor="blue")

    Plots a the `triangulation` using Plots library. The plots are
    appended to an existing plot, so a new figure must be created 
    before calling.
"""
function plotTriangle(triangulation::TriangulateIO;edgColor="black",verColor="blue")
    xs = triangulation.pointlist[1,:]
    ys = triangulation.pointlist[2,:]
    Plots.scatter!(xs,ys,color=verColor,legend=false)
    for edg in eachcol(triangulation.edgelist)
        if !(-1 in edg)
            xs = triangulation.pointlist[1,edg]
            ys = triangulation.pointlist[2,edg]
            Plots.plot!(xs,ys,color=edgColor,legend=false)
        end
    end
end


@doc """
    planeInterpolationFromMatrix(f::Function,vertex::Array{Float64,2})
    
    Given a vector function `f` and a (2,3) array of points `vertex`, it returns
    an array `res` such that the plane g(x) = res[1] + x[1]*res[2] + x[2]*res[3] 
    is the equation of the plane that passes through f(vertex[:,1]),
    f(vertex[:,2]) and f(vertex[:,3])
    
    It solves the problem by directly solving the linear system of the problem using
    Julia `\` operator. 

    Will return an error if vertex are colineal
"""
function planeInterpolationFromMatrix(f::Function,vertex::Array{Float64,2})::Array{Float64,1}
    if size(vertex) != (2,3)
        error("vertex size must be (2,3)")
    end
    farr = [f(x) for x in eachcol(vertex)]
    # creating matrix of equation
    A = hcat(ones(3,1),transpose(vertex))
    return A \ farr
end

@doc """
    planeInterpolationDirect(f::Function,vertex::Array{Float64,2})
    
    Given a vector function `f` and a (2,3) array of points `vertex`, it returns
    an array `res` such that the plane g(x) = res[1] + x[1]*res[2] + x[2]*res[3] 
    is the equation of the plane that passes through f(vertex[:,1]),
    f(vertex[:,2]) and f(vertex[:,3])
    
    It solves the problem by directly solving the using the formula for a plane given
    three of its points. Slower than `planeInterpolationFromMatrix`

    Will return an error if vertex are colineal
"""
function planeInterpolationDirect(f::Function,vertex::Array{Float64,2})::Array{Float64,1}
    if size(vertex) != (2,3)
        error("vertex size must be (2,3)")
    end
    farr = [f(x) for x in eachcol(vertex)]
    points = vcat(vertex,transpose(farr))
    normal = cross(points[:,2]-points[:,1],points[:,3]-points[:,1])
    return [dot(normal,points[:,1])/normal[3],-normal[1]/normal[3],-normal[2]/normal[3]]
end


@doc """
interpolateFunction(func::Function,triangulation::TriangulateIO)::Array{Float64,2}
    
    Returns an array `res` of size `(3,numberoftriangles(triangulate))` such that `res[:,i]` is 
    the `planeInterpolationFromMatrix` of function `func` for triangle `triangulate.trianglelist[:,i]`
"""
function interpolateFunction(func::Function,triangulation::TriangulateIO)::Array{Float64,2}
    ntri = size(triangulation.trianglelist)[2]
    res = Array{Float64,2}(undef,(3,ntri))
    for i in 1:ntri
        tlist = triangulation.trianglelist[:,i]
        res[:,i] = planeInterpolationFromMatrix(func,triangulation.pointlist[:,tlist])
    end
    return res
end
function evalInterpolationInTriangulation(triangulation::TriangulateIO,interpolation::Array{Float64,2},voronoiTree)::Array{Float64,1}
    return [evalInterpolation(triout.pointlist[:,j],interpolation,voronoiTree) for j in 1:size(triout.pointlist)[2]]
end

@doc """
    graphEdgeList(triangulation::TriangulateIO)::Array{Array{Int32,1},1}

    returns the edgeList of a given triangulation
"""
function graphEdgeList(triangulation::TriangulateIO)::Array{Array{Int32,1},1}
    final = []
    for n in 1:size(triangulation.pointlist)[2]
        res = []
        for vec in eachcol(triout.edgelist)
            if vec[1] == n
                push!(res,vec[2])
            elseif vec[2] == n
                push!(res,vec[1])
            end
        end
        push!(final,res)
    end
    return final
end

@doc """
    graphEdgeList(triangulation::TriangulateIO)::Array{Array{Int32,1},1}

    returns the edgeList of a given triangulation
"""
function triangleList(triangulation::TriangulateIO)::Array{Array{Int32,1},1}
    final = []
    for n in 1:size(triangulation.pointlist)[2]
        res = []
        for j in 1:size(triangulation.trianglelist)[2]
            if triangulation.trianglelist[1,j] == n || triangulation.trianglelist[2,j] == n || triangulation.trianglelist[3,j] == n
                push!(res,j)
            end
        end
        push!(final,res)
    end
    return final
end

@doc """
    sign(point1::Array{Float64,1},point2::Array{Float64,1},point3::Array{Float64,1})::Float64

    Auxiliary function for checking whether a point is inside a triangle
"""
function sign(point1::Array{Float64,1},point2::Array{Float64,1},point3::Array{Float64,1})::Float64
    return (point1[1]-point3[1])*(point2[2]-point3[2]) - (point2[1]-point3[1])*(point1[2]-point3[2])
end

@doc """
    pointInTriangle(point::Array{Float64,1},vertex::Array{Float64,2})::Bool

    function for veryfing if `point` is inside the triangle defined by `vertex`
"""
function pointInTriangle(point::Array{Float64,1},vertex::Array{Float64,2})::Bool
    if size(vertex)!= (2,3)
        error("bad vertex shape")
    end
    d1 = sign(point,vertex[:,1],vertex[:,2])
    d2 = sign(point,vertex[:,2],vertex[:,3])
    d3 = sign(point,vertex[:,3],vertex[:,1])
    P = (d1<0) || (d2 < 0) || (d3 < 0)
    Q = (d1>0) || (d2 > 0) || (d3 > 0)
    return !(P && Q)
end

@doc """
    findTriangleContaining(point::Array{Float64,1},triangulate::TriangulateIO,triangleList,voronoiTree)::Int32

    Given a `triangleList` object, such that `triangleList[i]` is an array of the indexes of triangles of which 
    `triangulate.pointlist[i]` is a vertex, and a NearestNeighbors tree `voronoiTree` that searches for closest 
    vertex in triangulation, the function returns the index of the triangle in `triangulate` closest to `point`.

    It returns `-1` if no such triangle is found.
"""
function findTriangleContaining(point::Array{Float64,1},triangulate::TriangulateIO,triangleList,voronoiTree)::Int32
    closestVertex = NearestNeighbors.knn(voronoiTree,point,1,false)[1][1]
    for tri in triangleList[closestVertex]
        trianglePoints = triangulate.pointlist[:,triout.trianglelist[:,tri]]
        if pointInTriangle(point,trianglePoints)
            return tri
        end
    end
    print(point)
    return -1
end

@doc """
    evalInterpolation(point::Array{Float64,1},triangulate::TriangulateIO,triangleList,interpolation::Array{Float64,2},voronoiTree)::Float64

    Returns the value of a function interpolated by array `interpolation` over a triangulateion `triangulate` with 
    triangle neighbor list `triangleList` on  `point`. 

    Will raise error if point is not inside of triangulation.
"""
function evalInterpolation(point::Array{Float64,1},triangulate::TriangulateIO,triangleList,interpolation::Array{Float64,2},voronoiTree)::Float64
    triangleindex = findTriangleContaining(point,triangulate,triangleList,voronoiTree)
    coefs = interpolation[:,triangleindex]
    return coefs[1] + coefs[2]*point[1] + coefs[3]*point[2]
end

@doc """
    massMatrix(triangulate::TriangulateIO)::Array{Float64,2}

    If builds the mass matrix using normal quadrature for a triangulation `triangulate`.
"""
function massMatrix(triangulate::TriangulateIO)::Array{Float64,2}
    np = size(triangulate.pointlist)[2]
    nt = size(triangulate.trianglelist)[2]
    M = zeros(np,np)
    for k in 1:nt
        Mtemp = 1/12*[2.0 1.0 1.0; 1.0 2.0 1.0; 1.0 1.0 2.0]* triangulate.triangleattributelist[1,k]
        M[triangulate.trianglelist[:,k],triangulate.trianglelist[:,k]] += Mtemp
    end
    return M
end

@doc """
    hatGradientInTriangle(k::Integer,triangulate::TriangulateIO)::Array{Float64,2}

    Calculates the gradient of the hat function for triangle `k` of triangulation `triangulate`
"""
function hatGradientInTriangle(k::Integer,triangulate::TriangulateIO)::Array{Float64,2}
    points = triangulate.pointlist[:,triangulate.trianglelist[:,k]]
    H =  [points[2,2]-points[2,3] points[2,3]-points[2,1] points[2,1]-points[2,2] ; -points[1,2]+points[1,3] -points[1,3]+points[1,1] -points[1,1]+points[1,2]]
    H *= 1/(2*triangulate.triangleattributelist[1,k])
end

@doc """
    localStiffness(k::Integer,triangulate::TriangulateIO)::Array{Float64,2}

    Builds the local stiffness Matrix for the constant unity scalar field over triangle `k` 
    in triangulation `triangulate`
"""

function localStiffness(k::Integer,triangulate::TriangulateIO)::Array{Float64,2}
    G = hatGradientInTriangle(k,triangulate)
    A = zeros(3,3)
    for i in 1:3
        for j in i:3
            A[i,j] = G[1,i]*G[1,j] + G[2,i]*G[2,j]
            A[j,i] = A[i,j]
        end
    end
    A *= triangulate.triangleattributelist[1,k]
    return A
end

@doc """
    stiffnessMatrix(triangulate::TriangulateIO)::Array{Float64,2}

    Uses `localStiffness` to build the stiffness matrix for an unity scalar field over triangulation `triangulate`
"""
function stiffnessMatrix(triangulate::TriangulateIO)::Array{Float64,2}
    np = size(triangulate.pointlist)[2]
    nt = size(triangulate.trianglelist)[2]
    M = zeros(np,np)
    for k in 1:nt
        Mtemp = localStiffness(k,triangulate)
        M[triangulate.trianglelist[:,k],triangulate.trianglelist[:,k]] += Mtemp
    end
    return M
end

@doc """
convectionMatrix(triangulate::TriangulateIO,bx::Array{Float64,1},by::Array{Float64,1})::Array{Float64,2}

    Builds the convection matrix over triangulation `triangulate` for the vector field defined by `bx` and `by`
"""
function convectionMatrix(triangulate::TriangulateIO,bx::Array{Float64,1},by::Array{Float64,1})::Array{Float64,2}
    np = size(triangulate.pointlist)[2]
    nt = size(triangulate.trianglelist)[2]
    C = zeros(np,np)
    for k in 1:nt
        bxmid = mean(bx[triangulate.trianglelist[:,k]])
        bymid = mean(by[triangulate.trianglelist[:,k]])
        aux = [bxmid bymid] * hatGradientInTriangle(k,triangulate)
        Ctemp = vcat(aux,aux,aux)*triangulate.triangleattributelist[1,k]/3
        C[triangulate.trianglelist[:,k],triangulate.trianglelist[:,k]] += Ctemp
    end
    return C
end

@doc """
    DFGbenchmarkMesh(nc::Integer;maxarea::Float64=0.05)

    Builds the normal and voronoi mesh for the DFG benchmark: a [0,2.2]x[0,0.41] rectagnular section with a 
    circle centered at [0.2,0.2] and radius 0.5. The circle is descripted by `nc` points and the 
    maximum area of each triangle is `maxarea`.

    `maxarea` is the parameter that changes the number of triangles in the grid.
"""
function DFGbenchmarkMesh(nc::Integer;maxarea::Float64=0.05)
    angles = range(0,stop=2*pi,length=nc)
    circle = hcat([0.2+0.05*cos(t) for t in angles],[0.2+0.05*sin(t) for t in angles])[1:(end-1),:]
    circleSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ])
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Array{Float64,2}(transpose([0.0 0.0 ; 
            0.0 0.41 ; 
            2.2  0.41 ; 
            2.2 0.0;
    ]))
    triin.pointlist = hcat(triin.pointlist,transpose(circle))
    triin.segmentlist=Array{Int64,2}(transpose([1 2 ; 2 3 ; 3 4 ; 4 1 ]))
    triin.segmentlist = hcat(triin.segmentlist,transpose(circleSegments))
    triin.segmentmarkerlist=vcat([1,1,1,1],[2 for i in 1:(nc-1)])    
    triin.holelist = transpose([0.2 0.2 ])
    area=@sprintf("%.15f",maxarea) # Don't use exponential format!
    (triout, vorout)=triangulate("pea$(area)DQv", triin)
    #fig = Plots.plot(aspect_ratio=1)
    #plotTriangle(triout)
    #display(fig)
    return triout,vorout
end

@doc """
    ellipse(nc::Integer;maxarea::Float64=0.05)

    Builds the normal and voronoi mesh for the DFG benchmark: a [0,2.2]x[0,0.41] rectagnular section with a 
    circle centered at [0.2,0.2] and radius 0.5. The circle is descripted by `nc` points and the 
    maximum area of each triangle is `maxarea`.

    `maxarea` is the parameter that changes the number of triangles in the grid.
"""
function ellipse(nc::Integer;maxarea::Float64=0.05)
    angles = range(0,stop=2*pi,length=nc)
    circle = hcat([0.2+0.05*cos(t) for t in angles],[0.2+0.08*sin(t) for t in angles])[1:(end-1),:]
    circleSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ])
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Array{Float64,2}(transpose([0.0 0.0 ; 
            0.0 0.41 ; 
            2.2  0.41 ; 
            2.2 0.0;
    ]))
    triin.pointlist = hcat(triin.pointlist,transpose(circle))
    triin.segmentlist=Array{Int64,2}(transpose([1 2 ; 2 3 ; 3 4 ; 4 1 ]))
    triin.segmentlist = hcat(triin.segmentlist,transpose(circleSegments))
    triin.segmentmarkerlist=vcat([1,1,1,1],[2 for i in 1:(nc-1)])    
    triin.holelist = transpose([0.2 0.2 ])
    area=@sprintf("%.15f",maxarea) # Don't use exponential format!
    (triout, vorout)=triangulate("pea$(area)DQv", triin)
    #fig = Plots.plot(aspect_ratio=1)
    #plotTriangle(triout)
    #display(fig)
    return triout,vorout
end

@doc """
    joukowskyTrans(C)

    Computes the conformal transformation A(C) = C + 1/C for a vector `C` 
"""
function joukowskyTrans(C)
    a = C[1]^2 + C[2]^2
    return [C[1]*(a+1)/a C[2]*(a-1)/a]
end

@doc """
    joukowskiWing(nc::Integer;maxarea::Float64=0.05)

    Builds the normal and voronoi mesh for a [0,2.2]x[0,0.41] rectagnular section with a 
    joukowski Wing centered around [0.2,0.2]. The wing is descripted by `nc` points and the 
    maximum area of each triangle is `maxarea`.

    `maxarea` is the parameter that changes the number of triangles in the grid.
"""
function joukowskiWing(nc::Integer;maxarea::Float64=0.05)
    r = 1.25
    c = [-0.2,0.2]
    angles = range(0,stop=2*pi,length=nc)
    circlex = [c[1]+r*cos(t) for t in angles]
    circley = [c[2]+r*sin(t) for t in angles]
    circle = hcat(circlex,circley)[1:(end-1),:]
    wing = vcat([joukowskyTrans(c) for c in eachrow(circle)]...)
    wing = 0.1*wing
    wing = vcat([transpose(c + [0.4, 0.2]) for c in eachrow(wing)]...)
    wingSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ]);
    triin=Triangulate.TriangulateIO()
    triin.pointlist=Array{Float64,2}(transpose([0.0 0.0 ; 
                0.0 0.41 ; 
                2.2  0.41 ; 
                2.2 0.0;
                ]))
    triin.pointlist = hcat(triin.pointlist,transpose(wing))
    triin.segmentlist=Array{Int64,2}(transpose([1 2 ; 2 3 ; 3 4 ; 4 1 ]))
    triin.segmentlist = hcat(triin.segmentlist,transpose(wingSegments))
    triin.segmentmarkerlist=vcat([1,1,1,1],[2 for i in 1:(nc-1)])
    triin.holelist = transpose([0.4 0.2 ])
    area=@sprintf("%.15f",maxarea) # Don't use exponential format!
    (triout, vorout)=triangulate("pea$(area)DQv", triin)
    return triout, vorout
end

@doc """
    DFGsolver(titer::Integer=100;nu = 0.001, nc::Integer=20,maxarea::Float64=0.001,meshFunc=DFGbenchmarkMesh)

    Using the mesh built by function `meshFunc(nc,maxarea)` for a fluid with kinematic viscosity `nu`. 
    For stability reasons, the default time step is set to `0.01` and the number of time iterations is `titer`.

    It returns a dictionary with the velocity field and the pressure.
"""
function DFGsolver(titer::Integer=100;nu = 0.001, nc::Integer=20,maxarea::Float64=0.001,meshFunc=DFGbenchmarkMesh)
    mesh = meshFunc(nc,maxarea=maxarea)
    addTriangleArea!(mesh[1])
    fig = Plots.plot(aspect_ratio=1)
    plotTriangle(mesh[1])
    Plots.xlims!(0,2.2)
    Plots.ylims!(0,0.41)
    Plots.savefig("mesh.pdf")
    np = size(mesh[1].pointlist)[2]
    outnodes = Array{Int32,1}()
    innodes = Array{Int32,1}()
    noslipnodes = Array{Int32,1}()
    obstaclenodes = [j for j in 1:numberofpoints(mesh[1]) if mesh[1].pointmarkerlist[j]==2]
    boundary = [j for j in 1:numberofpoints(mesh[1]) if mesh[1].pointmarkerlist[j]==1 || mesh[1].pointmarkerlist[j]==2]
    for j in boundary
        if mesh[1].pointlist[1,j] <= 0.01
            push!(innodes,j)
        elseif mesh[1].pointlist[1,j] >= 2.2
            push!(outnodes,j)
        else
            push!(noslipnodes,j)
        end
    end
    R = zeros(np,np)
    for node in outnodes
        R[node,node] = 10^6
    end
    println("$(length(innodes)) innodes: ")
    println(innodes)
    println("$(length(outnodes)) outnodes: ")
    println(outnodes)
    println("$(length(noslipnodes)) noslipnodes: ")
    println(noslipnodes)
    mask = ones(np)
    inSpeed = zeros(np)
    umax = 0.3
    for node in innodes
        y = mesh[1].pointlist[2,node]
        inSpeed[node] = 4*umax*y*(0.41-y)/(0.41^2)
        mask[node] = 0.0
    end
    for node in noslipnodes
        mask[node] = 0.0
    end
    dt = 0.01
    println("Reynolds Number = $(0.2*0.1/nu)")
    A = stiffnessMatrix(mesh[1])
    M = massMatrix(mesh[1])
    Minv = inv(M)
    #Minv = [1/M[i,i] for i in 1:size(M)[1]]
    Bx = convectionMatrix(mesh[1],ones(np),zeros(np))
    By = convectionMatrix(mesh[1],zeros(np),ones(np))
    U = zeros(np)
    V = zeros(np)
    P = zeros(np)
    resU = Array{Float64,2}(undef,(titer,np))
    resV = Array{Float64,2}(undef,(titer,np))
    resP = Array{Float64,2}(undef,(titer,np))
    for i in 1:titer
        println("Timestep: $i")
        C = convectionMatrix(mesh[1],U,V)
        #proposal step
        Uaux = U - Minv*(dt*(nu*A+C)*U)
        Vaux = V - Minv*(dt*(nu*A+C)*V)
        # enforce boundary conditions
        Uaux = Uaux.*mask + inSpeed
        Vaux = Vaux.*mask
        # solve PPE
        aux1 = (A+R)
        aux2 = -1/dt*(Bx*Uaux + By*Vaux)
        P = aux1 \ aux2
        # update velocity
        U = Uaux-Minv*(dt*(Bx*P))
        V = Vaux-Minv*(dt*(By*P))
        resU[i,:] = copy(U)
        resV[i,:] = copy(V)
        resP[i,:] = copy(P)
    end
    times = [i*dt for i in 1:titer]
    res = Dict("innodes"=>innodes,"outnodes"=>outnodes,"obstaclenodes"=>obstaclenodes,"mesh"=>mesh[1],"ts"=>times,"U"=>resU,"V"=>resV,"P"=>resP)
    return res
end
