
config=Struct() # instrument and user-defined setup parameters can go here
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
  @eval $(Symbol(z,x,:setup))(b::AbstractString,o...)=$(Symbol(z,:setup))(config,"/"*$y*"/"*b,o...)
end
