include("PhasePlot.jl")

# Function for a dual repression system without cooperativity

# x1: range of x1 values (i.e. X values)

# x2: range of x2 values (i.e. Z values)

# We use `@.` to apply the calculations across all rows.

# Note that model parameters are specified within the function

# Returns computed (dx1/dt, dx2/dt) over the range of (x1, x2)

function toggleMono(x1, x2)

    ax = 1.5            #constant production rate of x

    bx = 5.0            #scalar production rate of x

    nzx = 2.7           #hill coefficient of z binding to x

    nxz = 2.7           #hill coefficient of x binding to z

    zx = 0.4            #constant that controls repression of x by z

    xz = 1.5            #constant that controls repression of z by x

    del_z = 1           #degradation rate of z

    S = 105.0           #Signal concentration

    u = @. (ax+bx*S)/(1+S+((x2/zx)^nzx))-x1

    v = @. 1/(1+((x1/xz)^nxz))-del_z*x2  #eqn 12

    return (u,v)

end

#Range of x1, x2 values

x1range = (0,8,50)          #Has the form (min, max, num points)

x2range = (0,8,50)          #Has the form (min, max, num points)



x₀ = ([0.0,0.0],)  #initial state vectors; a common must be included after the first array

tspan=(0.0,100.0)             #time span

#Call the phaseplot functon to construct the phase portrait

phaseplot(toggleMono, x1range, x2range, xinit=x₀, t=tspan, clines=true,

        norm=true, scale=0.5)
