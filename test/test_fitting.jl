#using PyPlot

f(p,x)=p[1]/sqrt(2*pi)/p[3]*exp.(-(x-p[2]).^2/2/p[3]^2)+p[4];
xv=collect(linspace(-10,10,3e2));
yv=f([10,3,0.5,1],xv)+2rand(size(xv)...);

fit=nlfit([10,3,0.5,1],xv,yv,f)

#figure();plot(xv,f(fit.param,xv),color="blue");plot(xv,yv);
