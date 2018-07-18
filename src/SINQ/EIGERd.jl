EIGERload(o...)=EIGERd(SINQdataPath("eiger",o...))

#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
#abstract SINQd <: TripleAxisd
type EIGERd <: SINQd{0}
# TASP is a single-detector triple-axis neutron spectrometer at SINQ
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{Column,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    variables::Dict{AbstractString,Any}
    zeros::Dict{AbstractString,Any}
    data::Array{Measured,2}
    mask::BitArray
    instrument::Array{TripleAxis,1}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::AbstractString
    y::AbstractString
    #t::Int

    EIGERd()=new()
    EIGERd(d::EIGERd)=deepcopy(d)

    function EIGERd(filename::AbstractString,TASdef::Function=defineEIGER)
        this=new()
        this.filename=filename
        isfile(filename) || (warn("The file $filename does not exist; returning uninitialized structure");return)
        (this.command, this.header, this.columns, dataval, flatcone, this.variables, this.parameters, this.zeros)=ILLTASLoad(filename)
        isempty(flatcone) || warn("EIGER flatcone information ignored. Please extend my code to match new hardware capabilities.")
        columnnames=name.(this.columns) # this.columns is of type Array{Column,1}
        #(Note, map(Measured,dataval) would return Array{MeasuredSymmetric} and we need the abstract Array{Measured})
        thedata=Array{Measured}(size(dataval)...); thedata[:]=dataval[:] # promote dataval to Array{Measured}
        # give columns their uncertanties and units based on type [c]outners, [t]imers, [r].l.u., [e]nergies, [a]ngles
        alist=matchBA(columnnames,[r"^a[1-6]$",r"^[s,m,a]g[u,l]$"]) # angles (a1-a6,tilts)
        clist=matchBA(columnnames,"cnts").|matchBA(columnnames,r"^m[1-4]$") # counter and monitors
        dlist=matchBA(columnnames,r"^[s,m,a]t[u,l]$") # translations
        elist=matchBA(columnnames,r"^e[i,f,n]?$") # matches ei,ef,en,e
        flist=matchBA(columnnames,"mf") # magnetic field is always(?) "mf"
        klist=matchBA(columnnames,[r"^t[t,s]$",r"^temp?$"]) # non-exhaustive list of possible temperature "motors"
        rlist=matchBA(columnnames,r"^q?[h,k,l]$") # matches qh,qk,ql,h,k,l
        tlist=matchBA(columnnames,"time")
        wlist=matchBA(columnnames,r"^k[i,f]$") # ki, kf
        # countingerror! inserts variance equal to the counts; replaceerror! inserts symmetric variance equal to the passed error squared
        any(clist) && (countingerror!(view(thedata,Colon(),clist));      addunit!(view(this.columns,clist),"counts") )
        any(tlist) && ( replaceerror!(view(thedata,Colon(),tlist),1e-2); addunit!(view(this.columns,tlist),"seconds") )
        any(rlist) && ( replaceerror!(view(thedata,Colon(),rlist),1e-4); addunit!(view(this.columns,rlist),"rlu") )
        any(elist) && ( replaceerror!(view(thedata,Colon(),elist),5e-4); addunit!(view(this.columns,elist),"eV";power=-3) )
        any(alist) && ( replaceerror!(view(thedata,Colon(),alist),1e-2); addunit!(view(this.columns,alist),"°") )
        any(dlist) && ( replaceerror!(view(thedata,Colon(),dlist),1e-2); addunit!(view(this.columns,dlist),"m";power=-3) )
        any(klist) && ( replaceerror!(view(thedata,Colon(),klist),1e-3); addunit!(view(this.columns,klist),"K") )
        any(flist) && ( replaceerror!(view(thedata,Colon(),flist),1e-4); addunit!(view(this.columns,flist),"T") )
        any(wlist) && ( replaceerror!(view(thedata,Colon(),wlist),1e-4); addunit!(view(this.columns,wlist),"Å⁻¹") )
        this.data=thedata
        this.counters=columnnames[clist]
        this.timers=columnnames[tlist]
        this.motors=columnnames[.!(clist.|tlist)]
        # To start, all points should be unmasked (for fitting/plotting)
        this.mask=falses(size(this.data,1))
        # Determining default x and y columns is somewhat straightforward.
        this.y="cnts" # Default y should be "cnts" on TASP.
        this.x=guessscanvariableILLTAS(this.header,this.command,this.parameters,columnnames,this.data) # Default x can be determined from the entered command:
        this.definst=TASdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empy array for the fits
        return this
    end
end
EIGERd{T<:AbstractString}(fn::AbstractArray{T},TASdef::Function=defineEIGER)=map(x->EIGERd(x,TASdef),fn)

include("EIGERd2TripleAxis.jl")
