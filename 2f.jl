using DifferentialEquations
using Plots


#S = 35500, slightly below Hopf Bifurcation point

function toggleMono3(du,u,p,t)

ax,ay,bx,by,del_y,del_z,zx,yz,xz,xy,nzx,nxy,nxz,nyz,S = p

    du[1] = @. (ax+bx*S)/(1+S+((u[3]/zx)^nzx))-u[1] #eqn 11

    du[2] = @. (ay+by*S)/(1+S+((u[1]/xy)^nxy))-del_y*u[2]

    du[3] = @. 1/(1+((u[1]/xz)^nxz) + ((u[2]/yz)^nyz))-del_z*u[3]   #eqn 12
end

tspan1 = (0.0,20.0)

tspan2 = (20.00001,100)

u0 = [0.014524873383241137, 1.376946930379048, 0.0009157064354738723] #Cell at all SS values
u1 = u0*1.25                        #Cell at 25% above SS
u2 = u0*0.75                        #Cell at 25% below SS
#p1 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 105.0)
p2 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 100.0)


prob2 = ODEProblem(toggleMono3,u0,tspan2,p2)
sol2 = solve(prob2)


prob4 = ODEProblem(toggleMono3,u1,tspan2,p2)
sol4 = solve(prob4)


prob6 = ODEProblem(toggleMono3,u0,tspan2,p2)
sol6 = solve(prob6)

#plot(sol1,vars = (0,3), title = "Z vs. time", label = "Normal SS",xlims = (0,50),xaxis=("time"), yaxis=("Expression"))
plot(sol2, vars = (0,3), label = "Normal SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
#plot!(sol3, vars = (0,3), label = "+25% SS", xlims=(0,50),xaxis=("time"), yaxis=("Expression"))
plot!(sol4, vars = (0,3), label = "+25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
#plot!(sol5, vars = (0,3), label = "-25% SS", xlims=(0,50),xaxis=("time"), yaxis=("Expression"))
plot!(sol6, vars = (0,3), label = "-25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
#Uncomment below to see graph's of X,Y and Z vs. time for S = 0.02, 10.0, 10^5

#plot(sol1, title = "X vs. time", label = ["X(0.02)" "Y(0.02)" "Z(0.02)"],xaxis=("time"), yaxis=("Expression"))
#plot!(sol2, label = ["X(10)" "Y(10)" "Z(10)"])
#plot!(sol3,  label = ["X(5)" "Y(5)" "Z(5)"])
