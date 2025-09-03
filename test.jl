using QuadGK
PeriodIntegrand(x,thm)=(sqrt(2)/pi)*(thm/sqrt(cos(x*thm)-cos(thm)))
result=quadgk(x->PeriodIntegrand(x,pi/2),0,1)
println(result)