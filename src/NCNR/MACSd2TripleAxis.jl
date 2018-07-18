# TODO move defineMACS* to using the config dictionary

# A general function to create a TripleAxis object representative of a typical MACS setup.
# takes no input, returns a MACS TripleAxis object
function defineMACS(a::MACSd)
    MACS=TripleAxis(defineMACSsource(a),     defineMACSmonochromator(a),
                    defineMACSmonitor(a),    defineMACSsample(a),
                    defineMACSanalyzer(a),   defineMACSdetector(a),
                    defineMACSdistances(a),  defineMACScollimators(a)...,
                    defineMACSfixedenergy(a),defineMACSmode(a))
end
# FIXME!!! check all distances, mosaics, etc.
function defineMACSsource(a::MACSd)
    Source(30.,120.) # the cross-section of the incident beam, w×h in mm²
end
function defineMACSmonochromator(a::MACSd)
    tau=2pi/parse(a.parameters["monospacing"])
    extent=[3*50,5*25+4*1,2.] # overall w×h×t mm³; will be converted to second moments automatically
    mosaic=[pi,pi]/270 # 40'×40' = pi/270 × pi/270 radian
    #                                    vv this +1 is the monochromator scattering sense TODO check this
    monochromator=Bragg(tau,tauname(tau),+1,mosaic,extent,[horizontalyCurvedBragg,verticallyCurvedBragg]) ## TODO find a way to check focusing state?
end
function defineMACSmonitor(a::MACSd)
    Detector(50.,160.)# the cross-section of the post-monochromator beam monitor, w×h in mm² #TODO
end
function defineMACSsample(a::MACSd)
    # setup a default sample by reading values from the parameters Dict
    uv=eval(parse( "Float64["*replace(get(a.parameters,"orient","1 0 0 0 1 0")," ",";")*"]"))
    @assert 6==length(uv)
    #o1=uv[1:3]; o2=uv[4:6]
    o2=-uv[1:3]; o1=-uv[4:6] # FIXME this gives the same [h,k] orientation as Dave. Were the two orienting peaks somehow left-handed? Does the ng0 file record this?
    lat=eval(parse( "["*replace(get(a.parameters,"lattice","2pi 2pi 2pi 90 90 90")," ",";")*"]"))
    @assert 6==length(lat)
    ld=lat[1:3]; la=lat[4:6]
    shape=defineSampleShape(a,o1,o2,ld,la)
    mosaic=defineSampleMosaic(a)
    sample=Sample(Crystal(ld...,la...),o1,o2,shape,mosaic)
    resense!(sample,-1) # FIXME get(a.parameters,"ss",1)) # Does the MACS datafile contain the sample scattering sense?
    return sample
end
defineSampleShape(a::MACSd,o...)=diagm(SecondMoment{Float64}[100/12.,100/12.,100/12.]) # defaulting to 10×10×10 mm³
defineSampleMosaic(a::MACSd)=[0,0.] # defaulting to 0'×0' mosaic
function defineMACSanalyzer(a::MACSd)
    tau=2pi/parse(a.parameters["anaspacing"])
    extent=[3*50,5*25+4*1,2.] # overall w×h×t mm³; will be converted to second moments automatically FIXME lookup real values
    mosaic=[pi,pi]/270 # 40'×40' = pi/270 × pi/270 radian
    # check analyzer scattering sense. add double Bragg reflection mono/analyzer type? TODO
    analyzer=Bragg(tau,tauname(tau),+1,mosaic,extent,[flatBragg,verticallyCurvedBragg])
end
function defineMACSdetector(a::MACSd)
    Detector(30.,160.) # the cross-section of the detector, w×h in mm² FIXME
end
function defineMACSdistances(a::MACSd)
    [800,2e3,1170,750,1e3] # FIXME
end
function defineMACScollimators(a::MACSd)
    alphas=[60,90,90,120]
    betas=[600,600,600,600]
    #alphas[1]=-2 # m=2 guide before the monochromator
    #betas[1]=-2 # m=2 guide before the monochromator
    (alphas,betas)
end
function defineMACSfixedenergy(a::MACSd)
    fixede=eval(parse(get(a.parameters,"fixede","")))
    isa(fixede,Void)?3.7:fixede # FIXME
end
function defineMACSmode(a::MACSd)
#    ki_is_fixed=get(a.parameters,"fx",2)==2?false:true ## false==kf fixed, true==ki fixed
    false # FIXME
end

function fillinstrument!(a::MACSd)
    oneMACS=a.definst(a)
    (nopts,nodet,nocol)=size(a.data) # nodet should be 20 for MACS
    As=Array{Measured}(nopts,nodet,6)
    for aN in ("a1","a2","a3","a4","a5","a6")
      @assert isCol(a,aN) "Required Column $aN is missing!"
    end
    As[:,:,1]=getDat(a,"a1")/180*pi
    As[:,:,2]=getDat(a,"a2")/180*pi
    As[:,:,3]=getDat(a,"a3")/180*pi
    As[:,:,4]=getDat(a,"a4")/180*pi
    As[:,:,5]=getDat(a,"a5")/180*pi
    As[:,:,6]=getDat(a,"a6")/180*pi
    MACS=setangles(oneMACS,permutedims(value.(As),(3,1,2))) # setangles takes a 6xNxM Array of angles and returns a NxM array of TripleAxis objects
    a.instrument=MACS
end
