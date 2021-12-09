
function make_epileptor(weights)
    a = 1.0f0
    b = 3.0f0
    c = 1.0f0
    d = 5.0f0
    r = 0.00015f0
    s = 4.0f0
    Iext = 3.1f0
    slope = 0.0f0
    Iext2 = 0.45f0
    tau = 10.0f0
    aa = 6.0f0
    bb = 2.0f0
    Kvf = 0.0f0
    Kf = 0.01f0
    Ks = 0.01f0
    tt = 1.0f0
    modification = 0.0f0
    x0 = -2f0
    weights = Float32.(weights)

    function node(i,ydot,y,p,t,c_pop1,c_pop2)
        if y[1] < 0.0
            ydot[i,1] = - a * y[1] ^ 2 + b * y[1]
        else
            ydot[i,1] = slope - y[3] + 0.6 * (y[3] - 4.0) ^ 2
        end
        ydot[i,1] = tt * (y[2] - y[3] + Iext + Kvf * c_pop1 + ydot[i,1] * y[1])
        ydot[i,2] = tt * (c - d * y[1] ^ 2 - y[2])

        # energy
        if y[3] < 0.0
            ydot[i,3] = - 0.1 * y[3] ^ 7
        else
            ydot[i,3] = 0.0
        end
        if modification>0
            h =  x0 + 3/(1 + exp(-(y[1]+0.5)/0.1))
        else
            h = 4 * (y[1] - x0) + ydot[i,3]
        end
        ydot[i,3] = tt * (r * (h - y[3]  + Ks * c_pop1))

        # population 2
        ydot[i,4] = tt * (-y[5] + y[4] - y[4] ^ 3 + Iext2 + bb * y[6] - 0.3 * (y[3] - 3.5) + Kf * c_pop2)
        if y[4] < -0.25
            ydot[i,5] = 0.0
        else
            ydot[i,5] = aa * (y[4] + 0.25)
        end
        ydot[i,5] = tt * ((-y[5] + ydot[i,5]) / tau)

        # filter
        ydot[i,6] = tt * (-0.01 * (y[6] - 0.1 * y[1]))
    end

    function network(ydot,y,p,t)
        c_pop1 = dropdims(sum(weights .* (y[:,1] .- y[:,1]'); dims=2); dims=2)
        c_pop2 = dropdims(sum(weights .* (y[:,4] .- y[:,4]'); dims=2); dims=2)
        for i=1:size(weights,1)
            node(i,ydot,y[i,:],p,t,c_pop1[i],c_pop2[i])
        end
    end

    return network
end