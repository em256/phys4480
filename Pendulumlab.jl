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
function evolve(x0,dxdt,timerange,stepper=rk4step)
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

"""
    rk4step(x,dxdt,t,deltat)

Implements the RK4 stepping algorithm

`x` = the current variables at time `t`
`dxdt(x,t)` is a function of `x` and `t`
"""
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


