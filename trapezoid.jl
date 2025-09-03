using LinearAlgebra

"""
    trap_data(f,range)

stores the intermediate information for our trapezoid rule. 

Fields:
    f       # store the function we are integrating
    range   # store the range (a,b)
    val     # trapezoid rule integral at current depth, m
    m       # depth -- number of trapezoids = 2ᵐ

To increase the depth use command `subdivide!` or `subdivide`.

"""
mutable struct trap_data
    f       # store the function we are integrating
    range   # store the range (a,b)
    val     # trapezoid rule integral at current depth, m
    m       # depth -- number of trapezoids = 2ᵐ
    function trap_data(f,range)
        new(f,range,(f(range[1])+f(range[2]))*(range[2]-range[1])/2,0)
    end
end

function subdivide!(data ::trap_data)
    # Extract needed info from `data`
    f=data.f
    m=data.m+1
    range=data.range
    nextresult=sum_even_points(f,range,2^m)
    result=nextresult+data.val/2
    # Update data
    data.m=m
    data.val=result
end

function sum_even_points(f,range,n)
    (a,b)=range
    dx=(b-a)/n
    result=f(a+dx)*dx
    for j in 3:2:n-1
        result+=f(a+j*dx)*dx
    end
    return result
end

"""
`nrtrap(f,range;namedargs...)`

Uses the trapezoid rule algorithm from Numerical Recipies to calculate ``\\int_a^b f(x) dx``
where `range=(a,b)`.

`namedargs`
- `eps=1E-10`: algorithm terminates when relative change in integral falls below `eps`
- `max_divisions=20`: make no more than `max_divisions` subdivisions
- `debug=false`: set to `true` to see debug informations
"""
function nrtrap(f,range;eps=1E-10,max_divisions=20,debug=false)
    data=trap_data(f,range)
    oldval=zero(data.val)
    val=zero(data.val)
    for j in 1:max_divisions
        oldval=data.val
        subdivide!(data)
        val=data.val
        if norm(val-oldval)<eps*norm(oldval) || (val==0 && oldval==0)
            if debug
                print("converged after "*string(j)*" iterations")
            end
            return val
        end
    end
    error("Did not converge to precision "*string(eps)*
          "after "*string(max_divisions)*" subdivisions."*
          "previous 2 results: ("*string(oldval)*
          ","*string(val)*")")
end