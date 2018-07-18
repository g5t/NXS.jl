abstract type SampleEnvironment end
abstract type TemperatureControl <: SampleEnvironment end
abstract type FieldControl <: SampleEnvironment end
abstract type TemperatureFieldControl <: SampleEnvironment end

type Cryostat <: TemperatureControl
    temperature::Real
end
type Magnet <: FieldControl
    field::Real
end
type VectorMagnet <: FieldControl
    field::Real
    direction::LatticeVector #?
end
type Cryomagnet <: TemperatureFieldControl
    temperature::Real
    field::Real
end
type VectorCryoMagnet <: TemperatureFieldControl
    temperature::Real
    field::Real
    direction::LatticeVector #?
end
#type MuPAD <: SampleEnvironment
#    # MuPAD stuff
#end
getT(a::Union{TemperatureControl,TemperatureFieldControl})=a.temperature
getT(a::FieldControl)=300 #K, room temperature
setT!(a::Union{TemperatureControl,TemperatureFieldControl},t::Real)=a.temperature=t;

getH(a::TemperatureControl)=0 
getH(a::Union{TemperatureFieldControl,FieldControl})=a.field
setH!(a::Union{TemperatureFieldControl,FieldControl},h::Real)=a.field=h;

# Aliases for various sample environments
const ILL5 = Cryostat
const ORI3 = Cryostat
const ILL2 = Cryostat

const MA06 = CryoMagnet
const MA09 = CryoMagnet
const MA15 = CryoMagnet

const MA02 = VectorCryoMagnet
const MA07 = VectorCryoMagnet
const MA11 = VectorCryoMagnet 

