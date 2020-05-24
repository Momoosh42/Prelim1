using DifferentialEquations
using Plots


#S = 0.5, slightly below Hopf Bifurcation point
#X,Y,Z = [0.001437877053527658, 0.5724866675720363, 0.00035486318379280115]
function toggleMono3(du,u,p,t)

ax,ay,bx,by,del_y,del_z,zx,yz,xz,xy,nzx,nxy,nxz,nyz,S = p

#Equation 2

    du[1] = @. (ax+bx*S)/(1+S+((u[3]/zx)^nzx))-u[1]

    du[2] = @. (ay+by*S)/(1+S+((u[1]/xy)^nxy))-del_y*u[2]

    du[3] = @. 1/(1+((u[1]/xz)^nxz) + ((u[2]/yz)^nyz))-del_z*u[3]
end

#Timespan when values are at initial steady state
tspan1 = (0.0,20.0)

#Timespan at which S has been changed to 100
tspan2 = (20.00001,500)

u0 = [0.0014378770535276562, 0.5724866675720368, 0.00035486318379280104] #Cell at all SS values
u1 = u0*1.25                        #Cell at 25% above SS
u2 = u0*0.75                        #Cell at 25% below SS
p1 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 0.5) #S = 0.5
p2 = (0.039, 0.043, 6.1, 5.7, 1.05, 1.04, 0.000013, 0.011, 0.12, 0.00079, 2.32, 2.0, 2.0, 2.0, 100.0) #S = 100



prob2 = ODEProblem(toggleMono3,u0,tspan2,p2)
sol2 = solve(prob2)


prob4 = ODEProblem(toggleMono3,u1,tspan2,p2)
sol4 = solve(prob4)


prob6 = ODEProblem(toggleMono3,u2,tspan2,p2)
sol6 = solve(prob6)

z_i1 = fill(u0[3], (20,1))
z_i2 = fill(u1[3], (20,1))
z_i3 = fill(u2[3], (20,1))

plot(z_i1,y=0:1:19, title = "Z vs. time", label = "Normal SS",xlims = (0,100),xaxis=("time"), yaxis=("Expression"))
plot!(sol2, vars = (0,3), label = "Normal SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
plot!(z_i2,y=0:1:19, label = "+25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
plot!(sol4, vars = (0,3), label = "+25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
plot!(z_i3,y=0:1:19, label = "-25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"))
plot!(sol6, vars = (0,3), label = "-25% SS", xlims=(0,100),xaxis=("time"), yaxis=("Expression"), legend = :best)
#Uncomment below to see graph's of X,Y and Z vs. time for S = 0.02, 10.0, 10^5

#plot(sol1, title = "X vs. time", label = ["X(0.02)" "Y(0.02)" "Z(0.02)"],xaxis=("time"), yaxis=("Expression"))
#plot!(sol2, label = ["X(10)" "Y(10)" "Z(10)"])
#plot!(sol3,  label = ["X(5)" "Y(5)" "Z(5)"])
