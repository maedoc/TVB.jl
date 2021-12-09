
function mpr_delays(nt, dt, r, V, weights,idelays,g,Delta,tau,eta,J,I)
    nl, nt_, nn, nk = size(r)
    rtau = 1 / tau
    Delta_rpitau = Delta / (pi * tau)
    Threads.@threads for k=1:nk
        for t=1:nt
            for i=1:nn
                @inbounds for l=1:nl
                    acc = 0.0
                    for j=1:nn
                        acc += weights[j, i]*r[l, 256 + t - idelays[j, i], j, k]
                    end
                    r_cl = g[l, k] * acc
                    rl = r[l, 256 + t, i, k]
                    Vl = V[l, 256 + t, i, k]
                    r_noisel = r[l, 256 + t + 1, i, k]
                    V_noisel = V[l, 256 + t + 1, i, k]
                    drl = rtau * (Delta_rpitau + 2 * Vl * rl)
                    dVl = 1/tau * ( Vl^2 - pi^2 * tau^2 * rl^2 + eta + J * tau * rl + I + r_cl ) 
                    r[l, 256 + t + 1, i, k] = rl + dt*drl + r_noisel
                    V[l, 256 + t + 1, i, k] = Vl + dt*dVl + V_noisel
                end
            end
        end
    end
end

function mpr_benchmark()
    using BenchmarkTools
    dt = 1e-4
    nt = 1000
    k = 16
    nl = 8
    g = reshape(range(0,stop=1,length=k*nl),(nl,k))

    r = randn(Float32, nl, 2000, 76, k) * 1e-5;
    V = randn(Float32, nl, 2000, 76, k) * 1e-5;
    @. V -= 2;

    pars = [1.0,1.0,-5.0,15.0,0.0]

    weights, lengths = VirtualBrain.conn76()


    weights = Float32.(weights);
    idelays = trunc.(Int32, lengths / 10.0 / 0.1);
    maximum(idelays)

    delays(1000, dt, r, V, weights, idelays, g, pars...)

    @time delays(1000, dt, r, V, weights, idelays, g, pars...)
    @time delays(1000, dt, r, V, weights, idelays, g, pars...)
    @time delays(1000, dt, r, V, weights, idelays, g, pars...)

    1000 * k * nl / 0.138 * 1e-6
    # matches numbers for numba until boundscheck=False 
end