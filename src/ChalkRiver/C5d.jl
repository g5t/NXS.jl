#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
"""
C5 is a single-detector triple-axis neutron spectrometer at Chalk River.
It can operate in a polarized mode, allowing for the data block to be two or higher dimensions.
The mutable type `C5d{N}` is a subtype of `TripleAxisd` with `N` designating the dimensionality
of the data block -- `N=2` is a "standard" experiment, `N>2` is likely a polarized experiment
of some sort. As with other `TripleAxisd<:Neutrond` instrument types, the first dimension
of the data block is the points in a scan and the last dimension is the columns recorded during
the scan.
"""
type C5d{N} <: TripleAxisd{0} #Fixme? Move this to C5d{0} and C5d{1} somehow?
# C5 is a single-detector triple-axis neutron spectrometer at Chalk River
# it can operate in polarized mode so the data block can be > 2D
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    data::Array{Measured,N}
    mask::BitArray
    instrument::Array{TripleAxis} # should be Array{TripleAxis,N-1} but we can't do math here :(
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString

#    C5d()=new()
#    C5d(d::C5d)=deepcopy(d)
end

function loadC5d(crdata::CR,run::Integer,definst::Function=defineC5)
    genericoutput=CRTASLoad(crdata,run)
    isa(genericoutput,Void)&&(return)
    (filename, command, header, columns, dataval, parameters)=genericoutput
    columnnames=name.(columns) # columns is of type Array{Column,1}
    #(Note, map(Measured,dataval) would return Array{MeasuredSymmetric} and we need the abstract Array{Measured})
    thedata=Array{Measured}(size(dataval)...); thedata[:]=dataval[:] # promote dataval to Array{Measured}
    # give columns their uncertanties and units based on type [c]outners, [t]imers, [r].l.u., [e]nergies, [a]ngles
    alist=matchBA(columnnames,[r"^2t[am]$",r"^p[sh]i$"]) # angles (2TM,PSI,PHI,2TA)
    clist=matchBA(columnnames,[r"^mon=",r"^sig="]) # counter and monitors
    dlist=matchBA(columnnames,r"^[s,m,a]t[u,l]$") # translations
    νlist=matchBA(columnnames,"nu") # C5 at least uses frequency (THz) instead of energy transfer :(
    elist=matchBA(columnnames,r"^e[i,f,n]?$") # matches ei,ef,en,e
    flist=matchBA(columnnames,r"[hv]f$") # magnetic fields all end in HF or VF
    ilist=matchBA(columnnames,r"^i[hv]f.+$") # all currents have IHFx or IVFx
    klist=matchBA(columnnames,r"^[a-z]temp$") # non-exhaustive list of possible temperature "motors"
    rlist=matchBA(columnnames,r"^z?eta$") # matches eta and zeta
    tlist=matchBA(columnnames,"time")
    wlist=matchBA(columnnames,r"^k[i,f]$") # ki, kf
    cln=repeat([Colon()],outer=ndims(dataval)-1) # to enable viewing for 2D, 3D (or ND) data
    # countingerror! inserts variance equal to the counts; replaceerror! inserts symmetric variance equal to the passed error squared
    any(clist) && (countingerror!(view(thedata,cln...,clist));      addunit!(view(columns,clist),"counts") )
    any(tlist) && ( replaceerror!(view(thedata,cln...,tlist),1e-2); addunit!(view(columns,tlist),"seconds") )
    any(rlist) && ( replaceerror!(view(thedata,cln...,rlist),1e-4); addunit!(view(columns,rlist),"rlu") )
    any(elist) && ( replaceerror!(view(thedata,cln...,elist),1e-4); addunit!(view(columns,elist),"eV";power=-3) )
    any(alist) && ( replaceerror!(view(thedata,cln...,alist),1e-3); addunit!(view(columns,alist),"°") )
    any(dlist) && ( replaceerror!(view(thedata,cln...,dlist),1e-2); addunit!(view(columns,dlist),"m";power=-3) )
    any(νlist) && ( replaceerror!(view(thedata,cln...,νlist),1e-4); addunit!(view(columns,νlist),"Hz";power=12) )
    any(klist) && ( replaceerror!(view(thedata,cln...,klist),1e-3); addunit!(view(columns,klist),"K") )
    any(flist) && ( replaceerror!(view(thedata,cln...,flist),1e-4); addunit!(view(columns,flist),"T") )
    any(ilist) && ( replaceerror!(view(thedata,cln...,ilist),1e-4); addunit!(view(columns,ilist),"A") )
    any(wlist) && ( replaceerror!(view(thedata,cln...,wlist),1e-4); addunit!(view(columns,wlist),"Å⁻¹") )
    data=thedata
    counters=columnnames[clist]
    timers=columnnames[tlist]
    motors=columnnames[.!(clist.|tlist)]
    # To start, all points should be unmasked (for fitting/plotting)
    mask=falses(size(data,1))
    # Determining default x and y columns is somewhat straightforward.
    y="sig=sig" # Default y should be "Sig=SIG" on C5.

    # XXX you left off here.
    # TODO? replace (2tm,psi,phi,2ta) by (2tm/2,2tm,psi,phi,2ta/2,2ta) and rename them a1-a6
    # TODO alternately, keep (2tm,psi.phi,2ta) but use fillinstrument/TASdef to correctly pull together the six defining angles

    N=ndims(data)
    intdims=(2:N-1...) # for averaging and squeezing higher-order data blocks
    intdata=squeeze(mean(value(data),intdims),intdims)
    mode=get(parameters,"mode",:none)
    if mode==:q # this is a Q-scan? look to eta, zeta and nu for default x
        cmpr=matchBA(columnnames,[r"^z?eta",r"^nu"])
    elseif mode==:r # guessing that r is rocking. maybe it's :a for angle?
        cmpr=alist
    else
        info("unknown C5 scan mode $mode; selecting biggest step from angles and pseudo-Q")
        cmpr=matchBA(columnnames,[r"^z?eta",r"^nu",r"^2t[am]$",r"^p[sh]i$"])
    end
    biggest=indmax(vec(mean(abs.(diff(intdata[:,cmpr],1)),1)))
    x=columnnames[cmpr][biggest]

    instrument=Array{TripleAxis}(0) # fillinstrument!(this)
    fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
    fits=FitResult[] # initialize an empy array for the fits

    out=C5d{N}(filename,command,header,columns,motors,counters,timers,parameters,data,mask,instrument,definst,fitres,fits,x,y)
    fillinstrument!(out)
    return out
end
