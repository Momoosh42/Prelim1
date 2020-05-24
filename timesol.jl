#Code for 2d.

using DifferentialEquations

using Plots


function toggleMono3(du,u,p,t)

ax,ay,bx,by,del_y,del_z,zx,yz,xz,xy,nzx,nxy,nxz,nyz,S = p

    du[1] = @. (ax+bx*S)/(1+S+((u[3]/zx)^nzx))-u[1] #eqn 11

    du[2] = @. (ay+by*S)/(1+S+((u[1]/xy)^nxy))-del_y*u[2]

    du[3] = @. 1/(1+((u[1]/xz)^nxz) + ((u[2]/yz)^nyz))-del_z*u[3]   #eqn 12
end

tspan = (0.0,50.0)
u0 = [0.0,0.0,0.0]
p1 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 0.02) #S = 0.02
p2 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 10.0)    #S = 10.0
p3 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 100000.0)    #S = 100000

#Solve the 3 ODEs
prob1 = ODEProblem(toggleMono3,u0,tspan,p1)
sol1 = solve(prob1)
prob2 = ODEProblem(toggleMono3,u0,tspan,p2)
sol2 = solve(prob2)
prob3 = ODEProblem(toggleMono3,u0,tspan,p3)
sol3 = solve(prob3)

#Plot the 3 ODEs
plot(sol1,vars = (0,1), title = "X vs. time", label = "S = 0.02",xlims = (0,50), ylims = (0,6),xaxis=("time"), yaxis=("Expression"))
plot!(sol2, vars = (0,1), label = "S = 10.0",xlims = (0,50),ylims = (0,6), xaxis=("time"), yaxis=("Expression"))
plot!(sol3, vars = (0,1),  label = "S = 10^5",xlims = (0,50),ylims = (0,6), xaxis=("time"), yaxis=("Expression"))

#Uncomment below to see graph's of X,Y and Z vs. time for S = 0.02, 10.0, 10^5

#plot(sol1, title = "X vs. time", label = ["X(0.02)" "Y(0.02)" "Z(0.02)"],xaxis=("time"), yaxis=("Expression"))
#plot!(sol2, label = ["X(10)" "Y(10)" "Z(10)"])
#plot!(sol3,  label = ["X(5)" "Y(5)" "Z(5)"])
