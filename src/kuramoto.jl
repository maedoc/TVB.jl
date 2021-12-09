"Kuramoto model parameter structure."
struct kp
    G::Float64
    σ::Float64
    w::Array{Float64,2}
    l::Array{Float64,2}
end

"Kuramoto order parameter"
kop(sol) = [abs(mean(exp.(1.0im .* u_t))) for u_t in sol.u];

"Kuramoto dfun."
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

"Kuramoto noise function."
function kg(du,u,h,p,t)
    du .= sqrt(2*p.σ)
end

"Kuramoto history function."
function kh(p,t;idxs=nothing)
    nn = size(p.w,1)
    x = ((1:nn)/nn*2*pi .+ t) .* 2 .* pi
    collect(idxs==nothing ? x : x[idxs])
end
