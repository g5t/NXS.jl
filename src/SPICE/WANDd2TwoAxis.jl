function defaultconfigWAND()
    src=Struct(); mono=Struct(); moni=Struct();
    det=Struct(); dist=Struct(); coll=Struct(); wand=Struct()
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
    wand[:source]=src; wand[:monochromator]=mono; wand[:monitor]=moni;
    wand[:detector]=det;       wand[:distances]=dist;
    wand["/collimators/horizontal"]=coll; wand["/collimators/vertical"]=coll
    return wand
end
reset_config_WAND()=config[:wand]=defaultconfigWAND()
reset_config_WAND()

# A general function to create a TripleAxis object representative of a typical WAND setup.
function defineWAND(a::WANDd)
    scansense=get(a.parameters,"sense","+++")

    deft=defaultconfigWAND()
    wand=config[:wand]
    #
    src=get(wand,:source,Struct())
    if !isa(src,Source)
        def=deft[:source]
        src=Source(get(src,:width,def[:width]),get(src,:height,def[:height]))
    end
    mono=get(wand,:monochromator,Struct())
    if !isa(mono,Bragg)
        def=deft[:monochromator]
        named =get(mono,:name  ,get(a.parameters,"monochromator","ge113"))
        tau   =get(mono,:tau   ,gettau(named))
        sense =get(mono,:sense ,'+'==scatsense[1]?+1:-1)
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
    moni=get(wand,:monitor,Struct())
    if !isa(moni,Detector)
        def=deft[:monitor]
        moni=Detector(get(moni,:width,def[:width]),get(moni,:height,def[:height]))
    end
    sample=get(wand,:sample,Struct())
    if !isa(sample,Sample)
        o1=get(sample,"/orient/x",get(a.parameters,"orientx",[1,0,0.]))
        o2=get(sample,"/orient/y",get(a.parameters,"orienty",[0,1,0.]))
        if haskey(sample,:crystal) && isa(sample[:crystal],Crystal)
            crystal=sample[:crystal]
        else
            lc=get(a.parameters,"latticeconstants",[2pi,2pi,2pi,90,90,90])
            ld=get(sample,"/lattice/lengths",lc[1:3])
            la=get(sample,"/lattice/angles", lc[4:6])
            crystal=Crystal(ld...,la...)
        end
        mosaic=get(sample,:mosaic,Struct())
        isa(mosaic,Struct) && (mosaic=[get(mosaic,:horizontal,0.),get(mosaic,:vertical,0.)])
        isa(mosaic,Number) && (mosaic=[mosaic,mosaic])
        sense=get(sample,:sense,'+'==scatsense[2]?+1:-1)
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
    det=get(wand,:detector,Struct())
    if !isa(det,Detector)
        def=deft[:monitor]
        det=Detector(get(det,:width,def[:width]),get(det,:height,def[:height]))
    end
    dist=get(wand,:distances,Struct())
    if isa(dist,Struct)
        def=deft[:distances]
        s2m=get(dist,:source_monochromator,def[:source_monochromator])
        m2s=get(dist,:monochromator_sample,def[:monochromator_sample])
        s2d=get(dist,:sample_detector,def[:sample_detector])
        m2m=get(dist,:monochromator_monitor,def[:monochromator_monitor])
        dist=[s2m,m2s,s2d,m2m]
    end
    isa(dist,Array) || error("config[\"/wand/distances\"] should be either a Struct or an array of the five distances not $(typeof(dist))")
    isa(dist,Vector)||(dist=vec(dist))
    hcoll=get(wand,"/collimators/hoizontal",Struct())
    if isa(hcoll,Struct)
        def=deft["/collimators/horizontal"]
        α1=get(hcoll,:source_monochromator,def[:source_monochromator])
        α2=get(hcoll,:monochromator_sample,def[:monochromator_sample])
        α3=get(hcoll,:sample_detector,def[:sample_detector])
        hcoll=[α1,α2,α3]
    end
    isa(hcoll,Array)||error("config[\"/wand/collimators/horizontal\"] should be either a struct or an array of the four collimations not $(typeof(hcoll))")
    isa(hcoll,Vector)||(hcoll=vec(hcoll))
    vcoll=get(wand,"/collimators/vertical",Struct())
    if isa(vcoll,Struct)
        def=deft["/collimators/vertical"]
        β1=get(vcoll,:source_monochromator,def[:source_monochromator])
        β2=get(vcoll,:monochromator_sample,def[:monochromator_sample])
        β3=get(vcoll,:sample_detector,def[:sample_detector])
        vcoll=[β1,β2,β3]
    end
    isa(vcoll,Array)||error("config[\"/wand/collimators/vertical\"] should be either a struct or an array of the four collimations not $(typeof(vcoll))")
    isa(vcoll,Vector)||(vcoll=vec(vcoll))

    f=get(wand,:fixed,Struct())
    lambda=get(f,:lambda,get(a.parameters,"lambda",1.5))
    fixedk=get(f,:k,get(a.parameters,"k",2pi/lambda))
    energy=get(f,:E,get(a.parameters,"E",ħ²2mₙ*fixedk^2))

    tad=TwoAxisBanana(0.,0.,0.,[0.],src,mono,moni,sample,det,dist,hcoll,vcoll,energy)
    return tad
end


# There are multiple problems with using Measured values here. Drop them for the time being.
function fillinstrument!(a::WANDd)
    oneWAND=a.definst(a)
    (nopts,nodatcol)=size(a.data)
    (nopts,nodet,nodetcol)=size(a.detector)
    As=Array{Measured}(nopts,3+nodet)
    for aN in ("m1","s1","s2")
        @assert isCol(a,aN) "Required Column $aN is missing!"
    end
    # using getDat would scale-up the two single-columns to the size of the detector, which we don't wants
    # this is precisely what get_data_or_detector_column is for
    As[:,1]    =get_data_or_detector_column(a,"m1")/180*pi
    As[:,2]    =51.5/180*pi # the monochromator take-off angle is fixed at HB-2C
    As[:,3]    =get_data_or_detector_column(a,"s1")/180*pi
    As[:,4:end]=get_data_or_detector_column(a,"s2")/180*pi
    # setangles takes a 4xN Array of angles and returns a Nx0 vector of TripleAxis objects
    wand=setangles(oneWAND,permutedims(value.(As),(2,1))) # XXX make sure all angles are in radian! XXX
    a.instrument=wand
end
