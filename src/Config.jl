# Sample information, collimations, and orientation can be read from the datafile.
# A more accurate TripleAxis object can be created when loading data files by first
# setting appropriate keys in the TASPsetup dictionary, e.g., from your script files:
#
#         using NXS
#         NXS.config["tasp/analyzer/curvature/horizontal"]=NXS.flatBragg
#         scanobject=TASPd(filename)
#
# would define TripleAxis objects for the points in the scan with flat analyzers (since
# the TASP analyzer is only horizontally focusing to begin with)
#
# If one forgets to update alf[1-4] in SICS but changes the collimators, the correct
# value(s) can be forced into the TripleAxis objects by, e.g.,
#
#        using NXS
#        NXS.config["tasp/collimators/monochromator/sample"]=40; # arc-minutes
#        NXS.config["tasp/collimators/sample/analyzer"]=80; # arc-minutes
#        scanobject=TASPd(filename)
#
# if the collimators in place were -/40'/80'/- but the SICSserver had [any other values]
#
# A last alternative is to create your own equivalent defineTASP() function and pass it
# to the TASPd constructor, e.g.,
#
# using NXS
# function mydefineTASP()
# ...
# end
# defaultTASPscanobject = TASPd(filename)
# modifiedTASPscanobject= TASPd(filename,mydefineTASP)

#TODO move config to a Structs.Struct, then, e.g., config["/tasp"] would return the TASP Structs.Struct
config=Dict{String,Any}() # user-defined setup parameters can go here
# accessible via NXS.config

function getsetup(setup::Dict,key::AbstractString)
    haskey(setup,key) || error("key not present!")
    setup[key]
end
function getsetup(setup::Dict,key::AbstractString,def)
    get(setup,key,def)
end
function getsetup(setup::Dict,base::AbstractString,key::AbstractString,def)
    get(setup,base*key,def)
end
function getsetup{T<:AbstractString}(setup::Dict,base::AbstractString,keys::Array{T},def)
    fullkeys=base.*keys
    foundkey=findfirst(map(x->haskey(setup,x),fullkeys))
    return foundkey>0?setup[fullkeys[foundkey]]:def
end
# further getsetups that add a unit? so (base,key[s],def,unit) that look for, e.g., base*key.*["/cm","/mm","/m"]

checksetup(s::Dict,k::AbstractString)=haskey(s,k)
checksetup(s::Dict,k::AbstractString,T)=haskey(s,k)&&isa(s[k],T)
checksetup(s::Dict,k::AbstractString,T,L::Integer)=checksetup(s,k,T)&&length(s[k])>=L

# instrument-specific variants of get and check setup
for z in (:get,:check), (x,y) in zip((:TASP,:C5,:DMC,:EIGER,:WAND,:FLEXX,:RITA2,:IN22),("tasp","c5","dmc","eiger","wand","flexx","rita2","in22"))
  @eval $(Symbol(z,x,:setup))(b::AbstractString,o...)=$(Symbol(z,:setup))(config,$y*"/"*b,o...)
end
