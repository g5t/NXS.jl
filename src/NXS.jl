# With precompile turned on, the first pre-compilation of NXS takes a long time
# since every module used gets individualy precompiled as well.
# If any changes are then made to the NXS code the precompiled file becomes "stale"
# and julia will automatically recompile when the module is next used.
# The time to recompile is probably variable, but in simple testing was <24 seconds.
# This time was similar to the "first" compilation time if the NXS.ji file was deleted
# without touching the precompiled files for the submodules.
# The precompilation saves time on later "using NXS" statements as the time to load 
# the precompiled information is <9 seconds, compared to <17 seconds if not precompiled.
# A significant fraction of this time seems to be loading PyPlot, which can't be precompiled
# due to it being a thin wrapper for python's Matplotlib.pyplot.
__precompile__(false) # XXX !!! precompile causes a Seg. Fault for *some* plotting routines?!?!?
"""
NXS is a module to implement functionality relevant for the analysis
of data from neutron (and x-ray) scattering instruments. In its current
form it is mostly limited to triple-axis neutron scattering instruments
but there is functionality (and defined instrument types) relevant to 
diffractometers and multi-analyzer/multi-detector spectrometers. Further
instrument types can be introduced without too much effort as needed.
"""
module NXS
    import Base: (==),(+),(-),(.+),(.-),(/),(*),(^),(./),(.*),(%),(<),(>),(<=),(>=),
    sqrt,cos,sin,tan,acos,asin,atan,atan2,log,log10,zero,one,prod,convert,
    sort,abs,isless,isapprox,isnan,sum,getindex,setindex!,cross,dot,norm,
    round,ceil,floor,size,ndims,convert,promote_rule,sum,bin,var,get,find

    using Measureds     # A package that defines measured values (with uncertainty) and operations on such values
    using Lattices, Lattices.SecondMoments      # contains all Lattice-related functions, plus Crystal and Sample definitions
    importall PyPlot2TikZ
    importall Lattices
    import Lattices: _multMv!, _multMv # why does this need to be separate from "importall"?

    export FittingFunction
    export FitResult
    export SINQafsPath,SINQload,TASPload,ILLTASkeyLines2Dict,ILLTASLoad,
           Scatterd,Neutrond,TripleAxisd,SINQd,TASPd,isCol,getIdx,getDat,getVal,
           getVar,getErr,setDat,setVal,setVar,setErr,getAvg,getSum,
           delCol,delCol!,addCol,addCol!,exchCol,moveCol,combineRepeatedPoints,
           bin,bin!,mask,mask!,plot
    export Source,CircSource,ElliSource,RectSource,Detector,CircDetector,
           ElliDetector,RectDetector,getradius,getwidth,getheight,setradius,
           setwidth,setheight
    export hasbeenfit,gaussian,absgaussian,lorentzian,compatible,ndgrid,ndgridvecs
    export Bragg,flatBragg,verticallyCurvedBragg,horizontalyCurvedBragg,
           setHorzC!,setVertC!,setHorzR!,setVertR!,getHorzC,getVertC,getHR,
           getVR,gettau,tauname,resense!,setsense,getsense,getTh2Th,getwidth,getheight,
           getshape,TripleAxis,calck,getki,getkf,calcEs,calcEi,calcEf,getphi,
           iskiFixed,iskfFixed,getQ,getE,getQE,getangles,calcAngles,gotoQE!,
           gotoQE,getmonoHR,getmonoVR,getanaHR,getanaVR,getαs,getβs,
           geometricαs,geometricβs,geteffictiveαs,geteffectiveβs,
           getηs,getηh,getηv,getx,gety,getz,get3vector,get4vector,
           lab2sample,lab2sample!,sample2lab
    export ResFloat,ResMatrix,ResRM,coopernathansR0M,coopernathansRmM,
           popoviciR0M,popoviciRmM,convertcoopernathansR02Rm,convertpopoviciR02Rm
    export voigt,convres,gausshermite,
           convolute,evaluate,convoluteδE,convoluteδQω
    export FitResult,checkcorrelations,nlfitconvolution,nlfit,nlfitconvolution!,nlfit!

    using Colors        # for manipulation of plotting colors
    import OptimBase,LsqFit  # for fitting -- import does not merge namespaces (contrary to its name, it does not include exported functions)
    import Distributions # for TDist (used in FitResult.jl)
    using FastGaussQuadrature # for gausshermite (used in Convolution.jl)
    using SpecialFunctions # for, e.g., erfcx (the Faddeeva function) -- used for accurately approximating the Voigt function

    using SimpleNamedArrays # Provides NArray and Column types
    using Structs # Provides Struct for NXS.config

    # maybe it's useful to define some physical constants?
    # TODO Consider implementing these following Unitful.jl. Then, with the definition
    #   @unit Å "Å" Angstrom 1//10^10*u"m" false
    # we can do things like:
    #   const _speed_of_light = u"c0" |> Å/u"s" # u"c" is itself a unit?
    #   const meV = u"meV"
    #   const mₙ = u"mn" |> meV*u"s"^2/Å^2
    #   const kᴮ = u"k" |> meV/u"K"
    #   const ħ = u"ħ" |> meV*u"s"
    #   const ħ²2mₙ = ustrip( ħ^2/2/mₙ )
    # one could then, e.g., utilize unit-analysis to prevent stupid mistakes
    const _electron_radius =Measured(2.8179403227,0.0000000019^2)*1e-15 # m
    const _neutron_mass_kg =Measured(1.674927471,0.000000021^2)*1e-27 # kg
    const _speed_of_light  =299792458*1e10                            # Å/s
    const _neutron_mass    =Measured(939.5654133,0.0000058^2)*1e9/_speed_of_light^2     # meV s²/Å²
    const _neutron_g_factor=Measured(-3.82608545,0.00000090^2)        # unitless
    const _neutron_γ= -_neutron_g_factor/2 # This is Shirane, Shapiro, and Tranquada's "gyromagetic ratio"
    const _planck=Measured(4.135667662,0.000000025^2)*1e-12 # meV s
    const _planck_reduced=_planck/2pi
    const _boltzmann=Measured(8.6173303,0.0000050^2)*1e-2 # meV K⁻¹
    # the following has its error stripped because it's used in multiple places where Measureds aren't supported
    const ħ²2mₙ=value(_planck_reduced^2/(2*_neutron_mass))
    #const ħ²2mₙ=2.0721484 # ħ²/(2mₙ) = 2.0721484(3) meV Å² (ħ,² and ₙ used to avoid naming conflicts)

    global DEBUG=false
    debug(newdebug::Bool=DEBUG)=(global DEBUG=newdebug)
    testfield(test,field,a,b)=test(getfield(a,field),getfield(b,field))||(DEBUG&&(status(:debug,"$field mismatch");false))
    include("StatusMessages.jl")
    modfiles=["Config_new",
              "AbstractTypes",
              "SimpleDistributions",
              "FittingFunction",
              "FitResult",
              "SourceDetector",
              "Bragg",
              "TripleAxis",
              "TwoAxis",
              "TwoAxisBanana",
              "Resolution",
              "Convolution",
              "Fitting",
              "Utilities",
              "Scatterd",
              "ILLFormat",
              "SINQ/SINQ",
              "NCNR/MACSd",
              "ChalkRiver/CRDATA",
              "ChalkRiver/C5d",
              "GetSetParameters",
              "libnxs",
              "BinND",
              "Slice",
              "SPICE/SPICE",
              "SLS/X04SAd",
              "HZB/FLEXX",
              "ILL/IN22",
              "Plotting",
              "OperatorOverloading",
              "NCNR/MACSd2TripleAxis",
              "ChalkRiver/C5d2TripleAxis",
              "FunctionOnData",
              "autofit",
              "FormFactor",
              "Spurions",
              "Absorption"].*".jl"
    for file in modfiles
        DEBUG && status(:debug,"including $file")
        include(file)
    end
end
