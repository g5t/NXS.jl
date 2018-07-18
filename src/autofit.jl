export autofit!,getValDrift,peakintensity
function autofit!(a::Scatterd{0};p::Array=[],peaksat=[],func::Function=gaussian,maxpeaks=1,sortpositions=true)
    # every Scatterd{0} object is the data from a point detector.
    # often scans are performed across a peak which can be approximated
    # by a gaussian.
    x=getVal(a,a.x); y=getVal(a,a.y);
    if !isempty(peaksat)
        if !isempty(p)
            pinp=[any(p[2:3:end-1].==z) for z in peaksat];
            peaksat=peaksat[.!pinp]
        end
        pae=vcat([estimate_gaussian_parameters(x,y,z) for z in peaksat]...)
        p = isempty(p)? [pae;minimum(y)] : [p[1:end-1];pae;p[end]]
    end
    if !isempty(p) # try an initial fitting
        fr=nlfit!(func,p,a;fitinstrument=false)
        p0=fr.p
        yfit=func(p0,x);
        ykep= (yfit-p0[end]) .< 1e-3*(maximum(yfit)-minimum(yfit))
    else
        p0=[minimum(y)]
        yfit=p0[1]+zeros(y);
        ykep=trues(size(yfit))
    end
    while cld(length(p0)-1,3)<maxpeaks
        ym=y[ykep]-yfit[ykep]; # mask known peaks and remove the fit
        xm=x[ykep]
        length(ym) > length(p0) || break;
        imax = findfirst(ym.==maximum(ym))
        # add the new peak onto p0 (after any existing peaks)
        p0=[p0[1:end-1];estimate_gaussian_parameters(xm,ym,xm[imax]);p0[end]]
        if sortpositions
            # p0 ~ [a1,p1,w1,a2,p2,w2,...,bkg] -- sort p0 on [p1,p2,p3...]
            pp=p0[2:3:end-1]
            sp=[vec(3*(sortperm(pp)-1)' .+ [1,2,3]);length(p0)]
            p0=p0[sp]
        end
        fr=nlfit!(func,p0,a;fitinstrument=false)
        p0=fr.p # pull out just the fit values
        cld(length(p0)-1,3)<maxpeaks || break;
        yfit=func(p0,x)
        ykep= (yfit-p0[end]) .< 1e-3*(maximum(yfit)-minimum(yfit))
    end
    return fr
end

function estimate_gaussian_parameters(x::Vector,y::Vector,x0::Number)
    @assert length(x)==length(y)
    s=sortperm(x)
    x=x[s]; y=y[s]-minimum(y);
    ipk = indmin( abs.(x-x0) )
    il=ipk; ir=ipk;
    imin=1; imax=length(x)
    yh=y[ipk]/2
    while y[il]>yh && il>imin+1; il-=1; end
    while y[ir]>yh && ir<imax-1; ir+=1; end
    fwhm=abs(x[ir]-x[il])
    return [fwhm*2*yh,x[ipk],fwhm]
end

function getValDrift(a::Scatterd{0},b)
    z=getVal(a,b);
    (n,x)=extrema(z);
    m=mean(z);
    Measured(m,(m-x)^2,(m-n)^2)
end

tdiff(a::NTuple{2,T}) where T = a[2]-a[1]
function peakintensity(x::Scatterd{0},d=tdiff(extrema(getVal(x,x.y))))
    @assert hasbeenfit(x) "peakintensity is only defined for fit scans"
    @assert x.fitres.fitfun.SF == gaussian "peakintensity is only defined for gaussian fit scans"
    p=2sqrt(log(2)/pi)*x.fitres.pM[1]/x.fitres.pM[3]
    uncertainty(p)>value(p)?Measured(0,d^2,0):p
end
function peakarea(x::Scatterd{0},d=tdiff(extrema(getVal(x,x.y))))
    @assert hasbeenfit(x) "peakintensity is only defined for fit scans"
    @assert x.fitres.fitfun.SF == gaussian "peakintensity is only defined for gaussian fit scans"
    p=x.fitres.pM[1]
    uncertainty(p)>value(p)?Measured(0,d^2,0):p
end
