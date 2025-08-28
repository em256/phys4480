"""
    evolve(x0,dxdt,timerange,stepper)

integrates a differential equation.

`x0` - the dynamical variables at the original time -- can be a number, a vector, a tuple, etc
`dxdt(x,t)` -- a function of `x` and `t`,  Should return an object of the same type as `x0`
            -- evaluates to the right hand side of the differential equation
`timerange` -- a tuple `(ti,tf,dt)` listing the initial and final time, as well as the timestep
`stepper`

    stepper(x,dxdt,t,deltat))

gives you the value of `x` at the next timestep
"""
function evolve(x0,dxdt,timerange,stepper=eulerstep)
    ti,tf,dt=timerange
    numsteps=floor(Int,(tf-ti)/dt)
    x=collect(x0) #converts to a vector
                  #works fine without "collect" if you restrict input to vectors
    result=Array{typeof(x)}(undef,numsteps+1)
    result[1]=x
    t=ti
    for j in 1:numsteps
        x=stepper(x=x,dxdt=dxdt,t=t,deltat=dt)
        result[j+1]=x
        t+=dt
    end
    return result
end

evolve(;x0,dxdt,timerange,stepper=eulerstep)=evolve(x0,dxdt,timerange,stepper)

function rk4step(x,dxdt,t,deltat)
    k1=dxdt(x,t)
    k2=dxdt(x+k1*(deltat/2),t+deltat/2)
    k3=dxdt(x+k2*(deltat/2),t+deltat/2)
    k4=dxdt(x+k3*deltat,t+deltat)
    return x+(k1+2*k2+2*k3+k4)*(deltat/6)
end

rk4step(;x,dxdt,t,deltat)=rk4step(x,dxdt,t,deltat)

function pendulum_dxdt(x,t)
    (theta,v)=x
    return [v,-sin(theta)]
end

"""
    projectile_trajectory(xv0,dxdt,dt)
takes an initial vector `xv0=[x,y,vx,vy]`, 
a function `dxdt` which gives the derivitive,
and a timestep `dt`, and uses Runge Kutta
to integrate the differential equation until `y` becomes negative

    dxdt(xv,t)
returns the vector [dx/dt,dy/dt,dvx/dt,dvy/dt]
"""
function projectile_trajectory(xv0,dxdt,dt)
    xv=collect(xv0)
    trajectory=[xv]
    times=[zero(dt)]
    t=dt
    while xv[2]>=0
        xv=rk4step(xv,dxdt,t,dt)
        push!(trajectory,xv)
        push!(times,t)
        t+=dt
    end
    return (times,trajectory)
end


"""
    projectile_dxdt(m,gamma,delta,g)

creates a function-like object:

    dxdt1= projectile_dxdt(1,0,1,9.8)

and

    dxdt1(xv,t)

is the right hand side of our pendulum differential equation.
Here `xv=(x,y,vx,vy)`  and `dxdt1` returns `(dxdt,dydt,dvxdt,dvydx(`
"""
struct projectile_dxdt
    m
    gamma
    delta
    g
end

projectile_dxdt(;m,gamma,delta,g)=projectile_dxdt(m,gamma,delta,g)

function (dxdt::projectile_dxdt)(xv,t)
    x,y,vx,vy=xv
    m=dxdt.m
    gamma=dxdt.gamma
    delta=dxdt.delta
    g=dxdt.g
    dragprefactor=(gamma/m)*(vx^2+vy^2)^((delta-1)/2)
    dvx=-dragprefactor*vx
    dvy=-g-dragprefactor*vy
    return [vx,vy,dvx,dvy]
end

(dxdt::projectile_dxdt)(;xv,t)=dxdt(xv,t)

"""
    projectile_trajectory(xv0,dxdt,dt)
takes an initial vector `xv0=[x,y,vx,vy]`, 
a function `dxdt` which gives the derivitive,
and a timestep `dt`, and uses Runge Kutta
to integrate the differential equation until `y` becomes negative

    dxdt(xv,t)
returns the vector [dx/dt,dy/dt,dvx/dt,dvy/dt]
"""
function projectile_trajectory(xv0,dxdt,dt)
    xv=collect(xv0)
    trajectory=[xv]
    times=[zero(dt)]
    t=dt
    while xv[2]>=0
        xv=rk4step(xv,dxdt,t,dt)
        push!(trajectory,xv)
        push!(times,t)
        t+=dt
    end
    return (times,trajectory)
end

function interpx0(xvp,xvf)
    xp,yp,vxp,vyp=xvp
    xf,yf,vxf,vyf=xvf
    return (xf*yp-xp*yf)/(yp-yf)
end

function distance(initialv,m,g,gamma,delta,dt)
    dxdt=projectile_dxdt(m=m,gamma=gamma,delta=delta,g=g)
    times,traj=projectile_trajectory([0,0,initialv...],dxdt,dt)
    interpx0(traj[end],traj[end-1])
end

distance(;initialv,m,g,gamma,delta,dt)=distance(initialv,m,g,gamma,delta,dt)

