#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
type IN22d <: TripleAxisd{0}
# IN22 is a single-detector triple-axis neutron spectrometer at ILL
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

    function IN22d(fn::AbstractString,cmd::AbstractString,hd::Vector{S1},cl::Vector{C},mt::Vector{S2},
                    cnt::Vector{S3},tmr::Vector{S4},par::Dict{S5},var::Dict{S6},zer::Dict{S7},dat::Array{M,2},
                    msk::BitArray,x::AbstractString,y::AbstractString,
                    instdef::Function=defineIN22) where {S1<:AbstractString,
                                                          S2<:AbstractString,
                                                          S3<:AbstractString,
                                                          S4<:AbstractString,
                                                          S5<:AbstractString,
                                                          S6<:AbstractString,
                                                          S7<:AbstractString,
                                                          C<:Column,M<:Measured}
        this=new();
        this.filename=fn
        this.command=cmd
        this.header=hd
        this.columns=cl
        this.motors=mt
        this.counters=cnt
        this.timers=tmr
        this.parameters=par
        this.variables=var
        this.zeros=zer
        this.data=dat
        this.mask=msk
        this.x=x
        this.y=y
        this.definst=instdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empy array for the fits
        return this
    end

    IN22d()=new()
    IN22d(d::IN22d)=deepcopy(d)
end
