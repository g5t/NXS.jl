#abstract Scatterd
#abstract Neutrond <: Scatterd
#abstract TripleAxisd <: Neutrond
type CAMEAd <: TripleAxisd{2}
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,Any}
    variables::Dict{AbstractString,Any}
    zeros::Dict{AbstractString,Any}
    data::NArray{Measured,2,Column}
    detector::NArray{Measured,4,Column}
    mask::BitArray
    instrument::Array{TripleAxis,3}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::Vector{AbstractString}
    y::AbstractString

    function CAMEAd(fn::AbstractString,cmd::AbstractString,hd::Vector{S1},
                       mt::Vector{S2},cnt::Vector{S3},tmr::Vector{S4},
                       par::Dict{S5},var::Dict{S6},zer::Dict{S7},
                       dat::NArray{M,2,C},det::NArray{M,4,C},msk::BitArray,
                       x::Vector{S8},y::AbstractString,
                       instdef::Function=defineRITA2) where {S1<:AbstractString,
                                                             S2<:AbstractString,
                                                             S3<:AbstractString,
                                                             S4<:AbstractString,
                                                             S5<:AbstractString,
                                                             S6<:AbstractString,
                                                             S7<:AbstractString,
                                                             S8<:AbstractString,
                                                             C<:Column,M<:Measured}
        this=new();
        this.filename=fn
        this.command=cmd
        this.header=hd
        this.motors=mt
        this.counters=cnt
        this.timers=tmr
        this.parameters=par
        this.variables=var
        this.zeros=zer
        this.data=dat
        this.detector=det
        this.mask=msk
        this.x=x
        this.y=y
        this.definst=instdef
        fillinstrument!(this)
        this.fitres=FitResult{Float16,Float16,Float16,Measured,Int8}() # as a small(ish) placeholder
        this.fits=FitResult[] # initialize an empy array for the fits

        # With the instrumens array defined ensure that Qh,Qk,Ql,En are in the detector block
        QE=getQE.(this.instrument)
        hkl=[gethkl(x[1]) for x in QE] # an Array{Array{T,1}(3)}(size(instrumet))
        avec=getMultiple(this.parameters,["ax","ay","az"],0)
        bvec=getMultiple(this.parameters,["bx","by","bz"],0)
        avecname="["*fancy_hkl(avec...)*"]"
        bvecname="["*fancy_hkl(bvec...)*"]"
        all(iszero,avec) || addname!(this.detector,dot.([avec],hkl),Column(avecname,unit="rlu"))
        all(iszero,bvec) || addname!(this.detector,dot.([bvec],hkl),Column(bvecname,unit="rlu"))
        (all(iszero,avec)||all(iszero,bvec)) || (this.x=[avecname,bvecname]) # set the orienting vectors as default plotting axes
        #isCol(this,"qh") || addname!(this.detector,[geth(x[1]) for x in QE],Column("qh",unit="rlu"));
        #isCol(this,"qk") || addname!(this.detector,[getk(x[1]) for x in QE],Column("qk",unit="rlu"));
        #isCol(this,"ql") || addname!(this.detector,[getl(x[1]) for x in QE],Column("ql",unit="rlu"));
        isCol(this,"en") || addname!(this.detector,[     x[2]  for x in QE],Column("en",unit="meV"));
        return this
    end

    CAMEAd()=new()
    CAMEAd(d::CAMEAd)=deepcopy(d)
end

function applycorrections!(a::CAMEAd;int_corr_csv::AbstractString="",weights_csv::AbstractString="",alive_csv::AbstractString="")
    raw=getDat(a,"intensity")
    (npt,na4,nen)=size(raw) # if CAMEA na4==31, nen==5

    if isfile(int_corr_csv)
        int_corr=readdlm(int_corr_csv,',')
        @assert size(int_corr)==(na4,nen)
        cormat=permutedims(reshape(int_corr,(na4,nen,1)),(3,1,2))
    else
	isempty(int_corr_csv) || status(:info,"Requested intensity correction file $int_corr_csv does not exist, no correction performed.")
        corrmat=ones(1,na4,nen)
    end
    if isfile(weights_csv)
        weights=readdlm(weights_csv,',')
        @assert size(weights)==(na4,1) # the weights are defined per a4 channel only
        wghmat=permutedims(reshape(weights,(na4,1,1)),(3,1,2))
    else
	isempty(weights_csv) || status(:info,"Requested intensity weight file $weight_csv does not exist, no correction performed.")
        wghmat=ones(1,na4,1)
    end
    if isfile(alive_csv)
        dead=readdlm(alive_csv,' ').==0 # why is alive.csv space delineated? it's a CSV!!!
        @assert size(dead)==(na4,nen)
        dedmat=repeat( permutedims(reshape(dead,(na4,nen,1)),(3,1,2)), outer=[npt,1,1] ) # BitArray .& requires equal size arrays
    else
	isempty(alive_csv) || status(:info,"Requested detector status file $alive_csv does not exist, no masking performed.")
        dedmat=falses(npt,na4,nen)
    end

    # from PlotMF, corrected_intensity = ( (raw_intensity - 0.6*background)/int_corr ) * weights
    # and the background field is "empty" by default
    setDat(a,"intensity", (raw ./ cormat) .* wghmat )

    # The dead pixels can be handled by the mask, which is true for bad pixels
    a.mask .|= dedmat

    return nothing
end

"""
    flat2CAMEA(columns,data,flat)
The ILL TAS data format was designed for 0D detector instruments while CAMEA has 5 1D detectors.
To shoehorn CAMEA data into a ILL TAS data file some (constant) information is discarded and all detector pixels are stored in a single "point".
This works for data storage but is undesirable for data visualization and manipulation.
This function transforms the 2D ILL TAS data block into a 3D CAMEA data block restoring information which was stripped.
"""
function flat2CAMEA(cols::Vector{C},data::Matrix{T},flat::Matrix{T}) where {T,C}
    na4=31
    nen=5
    @assert size(flat,2) == (na4+1)*nen "each row should have $((na4+1)*nen) values for CAMEA"
    npt=size(flat,1)
    if size(data,1)!=npt
        status(:info,"unequal data and flat arrays")
        npt=min(npt,size(data,1))
        data=data[1:npt,:]
        flat=flat[1:npt,:]
    end
    @assert size(data,1) == npt "the data and flat arrays must contain the same number of rows"
    cn=name.(cols)

    # TODO: move [a5,a6] calculations to defineCAMEA since they're constants
    τs=gettau("PG(002)")
    Efs=(5:9)/2
    a6s=2asind.(τs./(2sqrt.(Efs/ħ²2mₙ)))
    a6mat=repeat(reshape(a6s,(1,1,nen)),outer=(npt,na4,1)) # (npt,na4,nen)
    a5mat=a6mat/2

    a4offset=(-15:15)*2.5
    a4flag=matchBA(cn,"a4")
    a4mat=repeat(a4offset' .+ data[:,a4flag], outer=(1,1,nen)) # (npt,na4,nen)

    keepasis=.!(a4flag)

    detectordata=cat(4,a4mat,a5mat,a6mat)
    detectorcols=["a4","a5","a6"]

    #keepasis.&=.!(matchBA(cn,"cnts")) # discard total counts
    #for ct in cnttmr[cnttmr.!="cnts"] # and ensure it's left off here too
    #    flag=matchBA(cn,ct)
    #    mat=repeat( data[:,flag], outer=(1,na4,nen) )
    #    detectordata=cat(4,detectordata,mat)
    #    detectorcols=vcat(detectorcols,ct)
    #    keepasis.&=.!(flag)
    #end

    detector=permutedims(reshape(flat,(npt,nen,na4+1))[:,:,1:end-1],(1,3,2)) # (npt,na4,nen)
    detectordata=cat(4,detectordata,detector)
    detectorcols=cat(1,detectorcols,"intensity")

    return (cols[keepasis],data[:,keepasis],Column(detectorcols),detectordata)
end

# CAMEA is a CAMEA-detector triple-axis neutron spectrometer at SINQ
# CAMEAd0 holds data from a CAMEAd object in a single data block, this
# requires replicating 0-D information for each point and is very memory
# hungry as a result. This type should probably only be used by Slice, which
# removes as much unneccessary 0-D data as possible
type CAMEAd0 <: TripleAxisd{0}
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
    data::Array{Measured,4}
    #data::NArray{Measured,3,Column} # -- hopefully this works
    mask::BitArray
    instrument::Array{TripleAxis,3}
    definst::Function
    fitres::FitResult
    fits::Vector{FitResult}
    x::Vector{AbstractString}
    y::AbstractString
    #t::Int

    CAMEAd0()=new()
    CAMEAd0(d::CAMEAd0)=deepcopy(d)
end
zerotype(::CAMEAd)=CAMEAd0
