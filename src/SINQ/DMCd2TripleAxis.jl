function defaultconfigDMC()
    src=Struct(); mono=Struct(); moni=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); dmc=Struct()
    src[:width]=30.
    src[:height]=120.
    #
    mono[:mosaic]=pi/270. # 40' in radian
    mono[:width]=3*50. # 3×50mm crystals wide
    mono[:height]=5*25.+4*1. # 5×25mm crystals + 4×1mm gaps tall
    mono[:depth]=2. # 2mm thick crystals
    mono["/curvature/horizontal"]=flatBragg
    mono["/curvature/vertical"]=verticallyCurvedBragg
    #
    moni[:width]=50.
    moni[:height]=160.
    #
    det[:width] =  1. #mm,
    det[:height]=160. #mm, height of standard detector mask
    #
    dist[:source_monochromator] = 800.
    dist[:monochromator_sample] =1400.
    dist[:sample_detector]      =1150.
    dist[:monochromator_monitor]=1000.
    #
    coll[:source_monochromator]= -2.
    coll[:monochromator_sample]=coll[:sample_detector]=800.
    #
    dmc[:source]=src; dmc[:monochromator]=mono; dmc[:monitor]=moni;
    dmc[:detector]=det;       dmc[:distances]=dist;
    dmc["/collimators/horizontal"]=coll; dmc["/collimators/vertical"]=coll
    return dmc
end
reset_config_DMC()=config[:dmc]=defaultconfigDMC()
reset_config_DMC()

# A general function to create a TripleAxis object representative of a typical DMC setup.
function defineDMC(a::DMCd)
    deft=defaultconfigDMC()
    dmc=config[:dmc]
    #
    src=get(dmc,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(dmc,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        tau   =get(mono,:tau   ,gettau("pg002"))
        named =get(mono,:name  ,tauname(tau))
        sense =get(mono,:sense ,get(a.parameters,"sm",+1))
        mosaic=get(mono,:mosaic,def[:mosaic])
        if isa(mosaic,Struct)
            horz=get(mosaic,"/horizontal",def[:mosaic])
            vert=get(mosaic,"/vertical"  ,def[:mosaic])
            mosaic=[horz,vert]
        end
        isa(mosaic,Array) || (mosaic=[mosaic,mosaic])
        width =get(mono,:width,def[:width])
        height=get(mono,:height,def[:height])
        depth =get(mono,:depth,def[:depth])
        extent=[width,height,depth] # overall w×h×t mm³; will be converted to second moments automatically
        cfh=get(mono,"/curvature/horizontal",def["/curvature/horizontal"])
        cfv=get(mono,"/curvature/vertical"  ,def["/curvature/vertical"]  )
        mono=Bragg(tau,named,sense,mosaic,extent,[cfh,cfv])
    end
    moni=get(dmc,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(dmc,:sample,Struct())
    if !isa(sample,Sample)
        o1=get(sample,"/orient/x",[1,0,0.])
        o2=get(sample,"/orient/y",[0,1,0.])
        if haskey(sample,:crystal) && isa(sample[:crystal],Crystal)
            crystal=sample[:crystal]
        else
            ld=get(sample,"/lattice/lengths",getMultiple(a.parameters,["as","bs","cs"],2pi))
            la=get(sample,"/lattice/angles", getMultiple(a.parameters,["aa","bb","cc"],90.))
            crystal=Crystal(ld...,la...)
        end
        mosaic=get(sample,:mosaic,Struct())
        isa(mosaic,Struct) && (mosaic=[get(mosaic,:horizontal,0.),get(mosaic,:vertical,0.)])
        isa(mosaic,Number) && (mosaic=[mosaic,mosaic])
        sense=get(sample,:sense,get(a.parameters,"ss",1))
        if haskey(sample,:geometry)&&isa(sample[:geometry],Lattices.SampleMesh)
            invAbsLen=get(sample,:inverse_absorption_length,0.)
            sample=Sample(crystal,o1,o2,sample[:geometry],mosaic,sense,invAbsLen)
        else
            shape=get(sample,:shape,Struct())
            if isa(shape,Struct)
                xx=get(shape,:x,100/12.) # <x*x> for a 10 mm cube
                yy=get(shape,:y,100/12.) # <y*y> for a 10 mm cube
                zz=get(shape,:z,100/12.) # <z*z> for a 10 mm cube
                shape=[xx,yy,zz]
            end
            isa(shape,Vector) && (shape=diagm(shape))
            shape=convert.(SecondMoment,shape) # a non-op if eltype(shape)==SecondMoment already
            sample=Sample(crystal,o1,o2,shape,mosaic,sense)
        end
    end
    det=get(dmc,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(dmc,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2d=get(dist,:sample_detector,def[:sample_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2d,m2m]
    end
    isa(dist,Array) || error("config[\"/dmc/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(dmc,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_detector,def[:sample_detector])
        hcoll=[α1,α2,α3]
    end
    isa(hcoll,Array)||error("config[\"/dmc/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(dmc,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_detector,def[:sample_detector])
        vcoll=[β1,β2,β3]
    end
    isa(vcoll,Array)||error("config[\"/dmc/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))

    f=get(dmc,:fixed,Struct())
    lambda=get(f,:lambda,get(a.parameters,"lambda",4.5))
    fixedk=get(f,:k,get(a.parameters,"k",2pi/lambda))
    energy=get(f,:E,get(a.parameters,"E",ħ²2mₙ*fixedk^2))

    tad=TwoAxis(src,mono,moni,sample,det,dist,hcoll,vcoll,energy)
    return tad
end

# There are multiple problems with using Measured values here. Drop them for the time being.
function fillinstrument!(a::DMCd)
    oneDMC=a.definst(a)
    nopts=size(a.data,1)
    anglecols=["a1","a2","a3","a4"]
    headAs=getMultiple(a.parameters,anglecols,0.)
    As=zeros(nopt,length(anglecols))
    for i=1:length(anglecols)
        As[:,i]=isCol(a,anglecols[i])?getVal(a,anglecols[i]):headAs[i]
    end
    # setangles takes a 4xN Array of angles and returns a Nx0 vector of TripleAxis objects
    dmc=setangles(oneDMC,permutedims(As/180*pi,(2,1))) # make sure all angles are in radian!
    a.instrument=dmc
end
