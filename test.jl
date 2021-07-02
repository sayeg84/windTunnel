include("./finiteElements.jl")

points = rand(2,60)
triin = Triangulate.TriangulateIO(pointlist = points)
triout, vorout = Triangulate.triangulate("veQ",triin);
println(size(triout.trianglelist))
println(size(vorout.pointlist))
trilist = triangleList(triout)
kdtree = NearestNeighbors.KDTree(triout.pointlist,NearestNeighbors.Euclidean())
testFunction(x) = sin(pi*x[1])*cos(pi*x[2])
interpolation = interpolateFunction(testFunction,triout)
n = 20
bx = 0.3
by = 0.6
xs = [x for x in range(bx,stop=by,length=n) for y in range(bx,stop=by,length=n)]
ys = [y for x in range(bx,stop=by,length=n) for y in range(bx,stop=by,length=n)]
zs = [testFunction([xs[i],ys[i]]) for i in 1:length(xs)]
fig = Plots.plot()
Plots.surface!(xs,ys,zs,alpha=0.5,camera=(0,90))
Plots.scatter!(xs,ys,zs,alpha=0.5)
#Plots.savefig(fig,"3Doriginal.pdf")
display(fig)


trianList = triangleList(triout)
disInter = [evalInterpolation(triout.pointlist[:,j],trianList,interpolation,kdtree) for j in 1:size(triout.pointlist)[2]]
fig = Plots.plot()
Plots.surface!(triout.pointlist[1,:],triout.pointlist[2,:],disInter,alpha=0.5,camera=(0,90))
Plots.scatter!(triout.pointlist[1,:],triout.pointlist[2,:],disInter,alpha=0.5)
display(fig)

zs= [evalInterpolation([xs[i],ys[i]],trianList,interpolation,kdtree) for i in 1:length(xs)]
fig = Plots.plot()
Plots.surface!(xs,ys,zs,alpha=0.5,camera=(0,90))
Plots.scatter!(xs,ys,zs,alpha=0.5)
display(fig)

findTriangleContaining([0.37894736842105264, 0.3],triout,trianList,kdtree)

nc = 22
angles = range(0,stop=2*pi,length=nc)
circlex = [0.2*cos(t) for t in angles]
circley = [0.2*sin(t) for t in angles]
circle = hcat(circlex,circley)[1:(end-1),:]
circleSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ])
triin=Triangulate.TriangulateIO()
triin.pointlist=Array{Float64,2}(transpose([-1.0 -1.0 ; 
            1.0 -1.0 ; 
            1.0  1.0 ; 
            -1.0 1.0;
            ]))
triin.pointlist = hcat(triin.pointlist,transpose(circle))
triin.segmentlist=Array{Int64,2}(transpose([1 2 ; 2 3 ; 3 4 ; 4 1 ]))
triin.segmentlist = hcat(triin.segmentlist,transpose(circleSegments))
triin.segmentmarkerlist=vcat([1,1,1,1],[2 for i in 1:(nc-1)])
triin.holelist = transpose([0.0 0.0 ])
area=@sprintf("%.15f",0.01) # Don't use exponential format!
(triout, vorout)=triangulate("pea$(area)DQ", triin)
#(triout, vorout)=triangulate("peDQ", triin)
fig = Plots.plot(aspect_ratio=1)
plotTriangle(triout)
display(fig)
nc = 40
r = 1.1
c = [-0.05,0.05]
angles = range(0,stop=2*pi,length=nc)
circlex = [c[1]+r*cos(t) for t in angles]
circley = [c[2]+r*sin(t) for t in angles]
circle = hcat(circlex,circley)[1:(end-1),:]
wing = vcat([joukowskyTrans(c) for c in eachrow(circle)]...)
#circleSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ])
fig = Plots.plot()
Plots.scatter!(wing[:,1],wing[:,2])
Plots.ylims!(-1,1)
Plots.xlims!(-4,4)
display(fig)
wingSegments = hcat([i+4 for i in 1:(nc-1)],[mod1(i+1,nc-1)+4 for i in 1:(nc-1) ]);

triin=Triangulate.TriangulateIO()
triin.pointlist=Array{Float64,2}(transpose([-4.0 -1.0 ; 
            -4.0 1.0 ; 
            4.0  1.0 ; 
            4.0 -1.0;
            ]))
triin.pointlist = hcat(triin.pointlist,transpose(wing))
triin.segmentlist=Array{Int64,2}(transpose([1 2 ; 2 3 ; 3 4 ; 4 1 ]))
triin.segmentlist = hcat(triin.segmentlist,transpose(wingSegments))
triin.segmentmarkerlist=vcat([1,1,1,1],[2 for i in 1:(nc-1)])
triin.holelist = transpose([0.0 0.0 ])
area=@sprintf("%.15f",0.03) # Don't use exponential format!
(triout, vorout)=triangulate("pea$(area)DQv", triin)
#(triout, vorout)=triangulate("peDQ", triin)
fig = Plots.plot()
Plots.scatter!(wing[:,1],wing[:,2])
plotTriangle(triout)
Plots.ylims!(-1,1)
Plots.xlims!(-4,4)
display(fig)


#=
points = rand(2,40)
triin = Triangulate.TriangulateIO(pointlist = points)
triout, vorout = Triangulate.triangulate("veQ",triin);
println(size(triout.trianglelist))
println(size(vorout.pointlist))
#barycenters = barycenterArray(triout)
#vorout.pointlist = barycenters
fig = Plots.plot(aspect_ratio=1)
plotTriangle(triout,edgColor="black",verColor="blue")
plotTriangle(vorout,edgColor="green",verColor="orange")
Plots.xlims!(0,1)
Plots.ylims!(0,1)
display(fig)

xs = [x for x in range(0,stop=1,length=201) for y in  range(0,stop=1,length=201)]
ys = [y for x in range(0,stop=1,length=201) for y in  range(0,stop=1,length=201)]
xs = [xs[i]+0.5*ys[i] for i in 1:length(xs)]
points = transpose(hcat(xs,ys))
triin = Triangulate.TriangulateIO(pointlist = points)
triout, vorout = Triangulate.triangulate("veQ",triin);
println(size(triout.trianglelist))
println(size(vorout.pointlist))
#fig = Plots.plot(aspect_ratio=1)
#plotTriangle(triout,edgColor="black",verColor="blue")
#plotTriangle(vorout,edgColor="green",verColor="orange")
#Plots.scatter!([0.23],[0.05])
#Plots.xlims!(0,1.5)
#Plots.ylims!(0,1)
#display(fig)

kdtree = NearestNeighbors.KDTree(vorout.pointlist,NearestNeighbors.Euclidean())
balltree = NearestNeighbors.BallTree(vorout.pointlist,NearestNeighbors.Euclidean())
brutetree = NearestNeighbors.BruteTree(vorout.pointlist,NearestNeighbors.Euclidean())
=#

# @benchmark  knn(balltree,[0.23,0.05],1,false)
