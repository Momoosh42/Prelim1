using DifferentialEquations
using Plots

#Function - toggleMono3
#du: derivatives of state variables; [dX/dt, dZ/dt]
#u: state variables; [X, Z]
#p: parameter value; S
#t: time span; (start_t, end_t)
#Here, du will be updated with newly computed
#values of dX/dt and dZ/dt
function toggleMono3(du,u,p,t)

ax,bx,del_z,zx,xz,nzx,nxz,S = p

#bistable switch equations
    du[1] = @. (ax+bx*S)/(1+S+((u[2]/zx)^nzx))-u[1]

    du[2] = @. 1/(1+((u[1]/xz)^nxz))-del_z*u[2]
end

#time range
tspan = (0.0,4975.0)

#Through the use of a Phase portrain producing script,
#Phaseplot.jl and ToggleMonoPrelim.jl
#(which is similar to PhasePortrait.jl and ToggleMono.jl, provided in HW #5)
#and visual inspection it was determined that there is an instability somewhere
#between S = 0.5 and S = 1.3.
#The first two loops solve the ODE from S = 0 to S = 2 with each one using a
#different initial condition. One initial condition lies far to the right, while
#the other lies far upwards. These IC's should converge when there is only one
#stable point, but should produce different steady state values when there are
#3 fixed points in the system.
#To determine the steady state value, the last point of the time dependent solve
#was utilized
#Each initial condition produces a distinct jump. These were found to be the
#limits of the 3 fixed point system, changing the initial values did not change
#the outer limits of where the system gains an unstable fixed point
#The limits where found by subtracting one solved value from the previous solved
#value and when the difference was greater than 1 the jump point was recorded
#Plots are produced for these two loops, with two jumps one indicating the
#beginning of the 3-fixed point system while the other indicating the end

#initial condition
u0 = [8.0,0.0]

SS = zeros(2560) #array of steady state values
x_val = zeros(2560) #array of x-axis values
j = 1 #counter
first_jump = 0 #jump point
for i in 0.001:0.001:2.56
p1 = (1.5, 5.0, 1.0, 0.4, 1.5, 2.7, 2.7, i)
prob1 = ODEProblem(toggleMono3,u0,tspan,p1)
sol1 = solve(prob1)

SS[j]=sol1[1,end] #records only the X values
if j>1 && (SS[j]-SS[j-1])>1
    global first_jump = i
end
x_val[j] = i #records the x-axis values
global j += 1
end


u1 = [0.0,8.0]

SS1 = zeros(2560) #array of steady state values
x_val = zeros(2560)
j = 1
second_jump = 0
for i in 0.001:0.001:2.56
p1 = (1.5, 5.0, 1.0, 0.4, 1.5, 2.7, 2.7, i)
prob1 = ODEProblem(toggleMono3,u1,tspan,p1)
sol1 = solve(prob1)

SS1[j]=sol1[1,end]
if j>1 && (SS1[j]-SS1[j-1])>1
    global second_jump = i-0.001
end
x_val[j] = i
global j += 1
end

println(first_jump)
print(second_jump)
plot(x_val,SS, label = "(8.0, 0.0) IC")
p1 = plot!(x_val,SS1, title = "X vs. S with jumps shown", xlabel = "S", ylabel = "X", label = "(0.0, 8.0) IC")


#Now that the S's that produces the jumps is known, the stable points can be plotted
#The first loop plots the stable points after the first jump
#The second loop plots the stable points before the second jump

u0_final = [8.0,0.0]

h = 0.001
last = 2.56
space = Int((last-first_jump+h)/h) #calculate number of elements needed in array
SS0_final = zeros(space) #array of steady state values
x_val0_final = zeros(space)
j = 1 #counter
for i in first_jump:h:last
p1 = (1.5, 5.0, 1.0, 0.4, 1.5, 2.7, 2.7, i)
prob1 = ODEProblem(toggleMono3,u0_final,tspan,p1)
sol1 = solve(prob1)
SS0_final[j]=sol1[1,end]

x_val0_final[j] = i
global j += 1
end


u1_final = [0.0,8.0]

space = Int(second_jump/h)
SS1_final = zeros(space) #array of steady state values
x_val1_final = zeros(space)
j = 1
for i in 0.001:h:second_jump
p1 = (1.5, 5.0, 1.0, 0.4, 1.5, 2.7, 2.7, i)
prob1 = ODEProblem(toggleMono3,u1_final,tspan,p1)
sol1 = solve(prob1)

SS1_final[j]=sol1[1,end]

x_val1_final[j] = i
global j += 1
end

plot(x_val0_final,SS0_final)
p2 = plot!(x_val1_final,SS1_final, title = "X vs. S only stable points shown", xlabel = "S", ylabel = "X")
plot(p1,p2, size = (800,800))
