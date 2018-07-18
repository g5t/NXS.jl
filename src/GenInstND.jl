type GenTAS0d <: TripleAxisd
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,AbstractString}
    #data::Array{Union{MeasuredSymmetric,MeasuredAsymmetric},2} ### Measured_asym.jl
    data::Array{Measured,2}                                     ### Measured.jl
    mask::BitArray
    instrument::Array{TripleAxis,1} # (npts,)
    definst::Function
    fitres::FitResult
    x::AbstractString
    y::AbstractString
end
type GenTAS1d <: TripleAxisd
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,AbstractString}
    #data::Array{Union{MeasuredSymmetric,MeasuredAsymmetric},3} ### Measured_asym.jl
    data::Array{Measured,3}                                     ### Measured.jl
    mask::BitArray
    instrument::Array{TripleAxis,2}# (npts,ndets)
    definst::Function
    fitres::FitResult
    x::Vector{AbstractString}
    y::AbstractString
end
type GenTAS2d <: TripleAxisd
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,AbstractString}
    #data::Array{Union{MeasuredSymmetric,MeasuredAsymmetric},4} ### Measured_asym.jl
    data::Array{Measured,4}                                     ### Measured.jl
    mask::BitArray
    instrument::Array{TripleAxis,3}# (npts,nxdets,nydets)
    definst::Function
    fitres::FitResult
    x::Vector{AbstractString}
    y::AbstractString
end
type GenTAS3d <: TripleAxisd
    filename::AbstractString
    command::AbstractString
    header::Array{AbstractString,1}
    columns::Array{AbstractString,1}
    motors::Array{AbstractString,1}
    counters::Array{AbstractString,1}
    timers::Array{AbstractString,1}
    parameters::Dict{AbstractString,AbstractString}
    #data::Array{Union{MeasuredSymmetric,MeasuredAsymmetric},5} ### Measured_asym.jl
    data::Array{Measured,5}                                     ### Measured.jl
    mask::BitArray
    instrument::Array{TripleAxis,4}# (npts,nxdets,nydets,nzdets)
    definst::Function
    fitres::FitResult
    y::AbstractString
end
