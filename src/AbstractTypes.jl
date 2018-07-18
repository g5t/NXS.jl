"""`Scatterd` is the highest echelon scattering data type."""
abstract type Scatterd{N} end
"""`Neutrond` is a subtype of `Scatterd` and should be the parent type for all
neutron scattering data types."""
abstract type Neutrond{N} <: Scatterd{N} end
"""`TripleAxisd` is a subtype of `Neutrond` and should be a parent type of all
triple-axis neutron scattering data types."""
abstract type TripleAxisd{N} <: Neutrond{N} end
"""`SINQd` is a subtype of `TripleAxisd` and is an abstract type type for all end
data from SINQ triple-axis neutron scattering instruments."""
abstract type SINQd{N} <: TripleAxisd{N} end
"""`Diffractometerd` is a subtype of `Neutrond` and should be a parent type of
all neutron diffractometers."""
abstract type Diffractometerd{N} <: Neutrond{N} end


abstract type Source end
abstract type Detector end
