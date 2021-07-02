include("./finiteElements.jl")
using ArgParse
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mesh"
            help = "geometry of the problem"
            arg_type = String
            default = "circle"
        "--titer"
            help = "Number of time iterations"
            arg_type = Int
            default = 100
        "--nu"
            help = "Kinematic viscosity"
            arg_type = Float64
            default = 1e-3
        "--nc"
            help = "Number of points used to describe obstable"
            arg_type = Int
            default = 20
        "--maxarea"
            help = "Maximum area for a triangle in the mesh. Do not use exponential format for this value"
            arg_type = Float64
            default = 0.001
        "--gifres"
            help = "Frecuency for gif"
            arg_type = Int
            default = 5
    end

    return parse_args(s)
end


println()
println("Initializing")
println()

args = parse_commandline()
if args["mesh"]=="circle"
    meshFunc = DFGbenchmarkMesh
elseif args["mesh"] == "wing"
    meshFunc = joukowskiWing
elseif args["mesh"] == "ellipse"
    meshFunc = ellipse
else
    error("No mesh for $(args["mesh"]) description")
end
res = DFGsolver(args["titer"],nu=args["nu"],nc=args["nc"],maxarea=args["maxarea"],meshFunc=meshFunc)
gifres = args["gifres"]

println()
println("Done")
println()

println()
println("Making pressure gif")
println()
anim = @animate for i in 1:length(res["ts"])
    fig = Plots.plot(aspect_ratio=1)
    Plots.surface!(res["mesh"].pointlist[1,:],res["mesh"].pointlist[2,:],res["P"][i,:],alpha=0.5,camera=(0,90),clim=(0,1),color=:viridis)
    Plots.title!("Pressure, t = $(round(res["ts"][i],digits=3))")
    Plots.xlims!(0,2.2)
    Plots.ylims!(0,0.41)
end every gifres

Plots.gif(anim,"pressure.gif",fps=30)

println()
println("Making U gif")
println()
anim = @animate for i in 1:length(res["ts"])
    fig = Plots.plot(aspect_ratio=1)
    Plots.surface!(res["mesh"].pointlist[1,:],res["mesh"].pointlist[2,:],res["U"][i,:],alpha=0.5,camera=(0,90),clim=(0.0,0.4),color=:curl)
    Plots.title!("U, t = $(round(res["ts"][i],digits=3))")
    Plots.xlims!(0,2.2)
    Plots.ylims!(0,0.41)
end every gifres

Plots.gif(anim,"Uspeed.gif",fps=30)

println()
println("Making V gif")
println()
anim = @animate for i in 1:length(res["ts"])
    fig = Plots.plot(aspect_ratio=1)
    Plots.surface!(res["mesh"].pointlist[1,:],res["mesh"].pointlist[2,:],res["V"][i,:],alpha=0.5,camera=(0,90),clim=(0.0,0.4),color=:curl)
    Plots.title!("V, t = $(round(res["ts"][i],digits=3))")
    Plots.xlims!(0,2.2)
    Plots.ylims!(0,0.41)
end every gifres

Plots.gif(anim,"Vspeed.gif",fps=30)

println()
println("Making Speed Norm gif")
println()
anim = @animate for i in 1:length(res["ts"])
    fig = Plots.plot(aspect_ratio=1)
    new = [res["U"][i,j]^2 +  res["V"][i,j]^2 for j in 1:size(res["U"])[2]]
    Plots.surface!(res["mesh"].pointlist[1,:],res["mesh"].pointlist[2,:],new,alpha=0.5,camera=(0,90),clim=(0.0,0.4),color=:curl)		

    Plots.title!("Speed norm, t = $(round(res["ts"][i],digits=3))")
    Plots.xlims!(0,2.2)
    Plots.ylims!(0,0.41)
end every gifres

Plots.gif(anim,"speedNorm.gif",fps=30)



