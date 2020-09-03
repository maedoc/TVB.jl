# not sure how Julia packages work yet..
# using BrainNetworkModels
using StochasticDelayDiffEq
# using StochasticDiffEq
using Statistics
using Plots
GR()
using DelimitedFiles
using ZipFile
# https://diffeq.sciml.ai/stable/features/progress_bar/#Using-Progress-Bars-Outside-Juno
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())
using DiffEqDevTools

function get_tvb_connectome()
    fname = "conn76.zip"
    if ~isfile(fname)
        url = "https://github.com/the-virtual-brain/tvb-data/raw/master/tvb_data/connectivity/connectivity_76.zip"
        download(url, fname)
    end
    zf = ZipFile.Reader(fname)
    # @show zf.files
    l = readdlm(zf.files[6])
    w = readdlm(zf.files[7])
    close(zf)
    w, l
end
w,l=get_tvb_connectome()

struct kp
    G::Float64
    σ::Float64
    w::Array{Float64,2}
    l::Array{Float64,2}
end

# Kuramoto order parameter
kop(sol) = [abs(mean(exp.(1.0im .* u_t))) for u_t in sol.u];

function kf(du,u,h,p,t)
    G = p.G/length(u)
    du .= 2*2*π/1000
    for j=1:length(u)
        for i=1:length(u)
            if i!=j
                uⱼ = h(p, t - p.l[i,j]; idxs=j)[1]
                du[i] += G * p.w[i,j] * sin(uⱼ - u[i])
            end
        end
    end
end

function kg(du,u,h,p,t)
    du .= sqrt(2*p.σ)
end

function kh(p,t;idxs=nothing)
    nn = size(p.w,1)
    x = ((1:nn)/nn*2*pi .+ t) .* 2 .* pi
    collect(idxs==nothing ? x : x[idxs])
end

function make_prob()
    # w, l = get_tvb_connectome()
    p = kp(0.1, 0.01, w, l)
    ic = kh(p, 0.0)
    nn = size(w,1)
    prob = SDDEProblem(kf, kg, ic, kh, (0.,1000.), p);
    prob
end


function compute_work_precision()
    # doesn't work yet, cache problem
    prob = remake(make_prob(),tspan=(0.0,1.0))
    reltols = 1.0 ./ 10.0 .^ (1:5)
    abstols = reltols#[0.0 for i in eachindex(reltols)]
    setups = [#Dict(:alg=>SRIW1())
            Dict(:alg=>EM(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1))
            #Dict(:alg=>RKMil(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
            #Dict(:alg=>SRIW1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
            #Dict(:alg=>SRA1(),:dts=>1.0./5.0.^((1:length(reltols)) .+ 1),:adaptive=>false)
            #Dict(:alg=>SRA1())
            ]
    names = ["EM"]#["SRIW1","EM"]#,"RKMil","SRIW1 Fixed","SRA1 Fixed","SRA1"]
    wp = WorkPrecisionSet(prob,abstols,reltols,setups;numruns=10,names=names,maxiters=1e7,error_estimate=:l2)
    plot(wp)
end
compute_work_precision()


function sweep_alg()
    # would be good to focus on strong order, since
    # e.g. in case of seizure, we want to strongly converge
    # to correct solution in case of specific noise path
    prob = make_prob()
    for alg=[LambaEM(), SRIW1(), SOSRA(), SOSRI(), SRA3(), DRI1(), SKenCarp()]
        @show alg
        solve(prob, alg, dt=1.0);
        for tol=[1e-3, 1e-4, 1e-5, 1e-6]
            @time sol = solve(prob, alg, abstol=tol);
        end
    end
    for alg=[EM(), EulerHeun(), SROCK1()]
        @show alg
        solve(prob, alg, dt=1.0);
        for dt=[1.0, 0.1]
            @time sol = solve(prob, alg, dt=dt);
        end
    end
end
sweep_alg()

# we can sweep but how to judge the accuracy of the solution?
# how to handle fixed step size?