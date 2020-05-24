using DifferentialEquations

using Plots

#This code is used to modify the S value, last value of p1, to find points near
#either the saddle point bifurcation or the hopf bifurcation
function toggleMono3(du,u,p,t)

ax,ay,bx,by,del_y,del_z,zx,yz,xz,xy,nzx,nxy,nxz,nyz,S = p

    du[1] = @. (ax+bx*S)/(1+S+((u[3]/zx)^nzx))-u[1] #eqn 11

    du[2] = @. (ay+by*S)/(1+S+((u[1]/xy)^nxy))-del_y*u[2]

    du[3] = @. 1/(1+((u[1]/xz)^nxz) + ((u[2]/yz)^nyz))-del_z*u[3]   #eqn 12
end

tspan = (0.0,4975.0)
u0 = [0.0,0.0,0.0]
p1 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 35500)


prob1 = ODEProblem(toggleMono3,u0,tspan,p1)
sol1 = solve(prob1)

print(sol1[end])

plot(sol1,vars = (0,3), label = "Z",xaxis=("time"), yaxis=("Expression"))
