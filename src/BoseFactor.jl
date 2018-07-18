kᴮ=8.6173303e-2 # meV/K

InverseBoseFactor(energy,temperature)=1-exp( -energy/(kᴮ*temperature) )
BoseFactor(energy,temperature)=1/InverseBoseFactor(energy,temperature)

function removeBoseFactor!(a::Neutrond;yname=a.y,ename="en",tname="tt")
    @assert isCol(a,yname)
    @assert isCol(a,ename)
    @assert isCol(a,tname)
    ibf=InverseBoseFactor.(a[ename],a[tname])
    setDat(a,yname,a[yname].*ibf)
    return nothing
end
