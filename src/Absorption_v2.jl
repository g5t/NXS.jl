struct Corner{T}
    data::SVector{3,T}
end
# a volume can be defined by its corners
# Move all of this to, e.g., the implementation of GeometricTypes?

# A 3-D object can be defined by its vertices, with three or more vertices
# defining a facet of the object. In general, a vertex will be used in
# defining more than one facet (e.g., a tetrahedron has 4 vertices and 4 3-vertex facets)
# so storing the positions of the vertices separate from the indicies of the
# facets is typically worthwhile.
# The GeometricTypes module(s) can do some of the heavy lifting, so we should
# utilize that work here. To make some work easier, restrict ourselves to the case
# of triangular facets
const GT=GeometryTypes
const GTp=GeometryTypes.Point{3,Float32}
const GTn=GeometryTypes.Normal{3,Float32}
const no_intersection= (false, 0., GTp(0.,0.,0.))
function ray_face_intersect(ri::GT.Point{3,T},rf::GT.Point{3,T},a::GT.Point{3,R},b::GT.Point{3,R},c::GT.Point{3,R}) where {T,R}
    # ri,rf -- coordinates for the start and end point of a ray
    # a,b,c -- coordinates for the three corners of the facet
    # n     -- normal direction for the facet (if a,b,c defined in the right order this is (b-a)×(c-a) )
    v1=b-a
    v2=c-a
    n=cross(v1,v2) # calculate the normal instead of relying on the mesh-stored normal
    n̂=n/norm(n)
    r=rf-ri # the ray vector
    x=( dot(a,n̂) - dot(ri,n̂) )/dot(r,n̂)
    0<=x<=1 || return no_intersection
    ip=ri+x*r # vector to the intersection point (from the origin)
    p=ip-a-dot(ip-a,n̂)*n̂ # in-facet-plane vector pointing from a to the intersection

    v̂1=v1/norm(v1)
    v̂2=v2/norm(v2)
    v̂1v̂2=dot(v̂1,v̂2)
    den=1-v̂1v̂2^2
    pv̂1=dot(p,v̂1)
    pv̂2=dot(p,v̂2)
    #sv1= ( pv̂1 - v̂1v̂2*pv̂2 )/den * v̂1
    #tv2= ( pv̂2 - v̂1v̂2*pv̂1 )/den * v̂2
    #s=norm(sv1)/norm(v1)
    #t=norm(tv2)/norm(v2)
    s= ( pv̂1 - v̂1v̂2*pv̂2 )/den/norm(v1)
    t= ( pv̂2 - v̂1v̂2*pv̂1 )/den/norm(v2)
    0<s<1 && 0<t<1 && s+t<1 || return no_intersection
    #return (true, x*norm(r)*sign(dot(r,n)), x*r+ri, a+sv1+tv2) # (flag, ray-length before hitting facet + out/in == (+,-), intersection point calculated two ways)
    return (true, x*norm(r)*sign(dot(r,n)), ip)
end
function ray_mesh_intersect(ri::GT.Point{3,T},rf::GT.Point{3,T},mesh::GT.HomogenousMesh) where T
    # a HomogenousMesh has fields :faces, :vertices, :normals
    nf=length(mesh.faces)
    all_intersections = Array{typeof(no_intersection)}(nf)
    for i=1:nf
        all_intersections[i]=ray_face_intersect(ri,rf,mesh.vertices[mesh.faces[i]]...)
    end
    return all_intersections
    #
    #return (does_intersect,all_intersections)
end
function ray_mesh_pathlength(ri::GT.Point{3,T},rf::GT.Point{3,T},mesh::GT.HomogenousMesh) where T
    all_intersections=ray_mesh_intersect(ri,rf,mesh)
    does_intersect = [x[1] for x in all_intersections]
    any(does_intersect) ? sum(x[2] for x in all_intersections[does_intersect]) : zero(T)
end

function transform(r::SMatrix{3,3},mesh::GT.HomogenousMesh)
    verts=similar(mesh.vertices)
    norms=similar(mesh.normals)
    for i=1:length(verts); verts[i]=r*mesh.vertices[i]; end
    for i=1:length(norms); norms[i]=r*mesh.normals[i];  end
    HomogenousMesh(verts,mesh.faces,norms,mesh.texturecoordinates,mesh.color,mesh.attributes,mesh.attribute_id)
    # or, non-type-preserving but cleaner:
    #HomogenousMesh([r].*mesh.vertices,mesh.faces,[r].*mesh.normals,mesh.texturecoordinates,mesh.color,mesh.attributes,mesh.attribute_id)
end
transform(r::Matrix,f)=transform(SMatrix{size(r)...}(r),f)

function axis_limits(m::GT.HomogenousMesh,v::GTp)
    extrema(dot.([v/norm(v)],m.vertices))
end
function axis_limits(m::GT.HomogenousMesh)
    v=m.vertices
    [ extrema(dot.([q],v)) for q in ( GTp(1,0,0), GTp(0,1,0), GTp(0,0,1) ) ]
end
function line_integration_limits(ri::GT.Point{3,T},rf::GT.Point{3,T},mesh::HomogenousMesh) where T
    all_intersections=ray_mesh_intersect(ri,rf,mesh)
    does_intersect = [x[1] for x in all_intersections]
    any(does_intersect) ? maximum(x[2] for x in all_intersections[does_intersect]) : zero(T)
end

# This absorption correction fuction requires that the facets defining the crystal
# have already been rotated to the correct orientation for the required Q point
# The incident beam is along the x-axis, and the outgoing beam is scattered by
# the angle tθ.
function absorption_correction(tθ::Real,μ::Real,mesh::HomogenousMesh;method::Symbol=:gausslegendre,points::Integer=10)
    o=GTp(0,0,0)
    x=GTp(1,0,0)
    y=GTp(0,1,0)
    z=GTp(0,0,1)
    din=GTp(-1,0,0)
    dout=GTp(cos(tθ),sin(tθ),0)
    (p,w)=gausslegendre(points)
    result=0.
    volume=0.
    mesh_limits=axis_limits(mesh)
    xmin,xmax=mesh_limits[1]
    ylim=mesh_limits[2]; zlim=mesh_limits[3]
    for i=1:points
        xi=(p[i]*(xmax-xmin)+(xmax+xmin))*x/2
        ymin=line_integration_limits(o+xi,o+xi+2*ylim[1]*y,mesh)*sign(ylim[1]) # since sign(ylim) determines the search direction
        ymax=line_integration_limits(o+xi,o+xi+2*ylim[2]*y,mesh)*sign(ylim[2])
        for j=1:points
            yj= (p[j]*(ymax-ymin) + (ymax+ymin))*y/2
            zmin=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[1]*z,mesh)*sign(zlim[1])
            zmax=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[2]*z,mesh)*sign(zlim[2])
            for k=1:points
                zk=(p[k]*(zmax-zmin)+(zmax+zmin))*z/2
                v=o+xi+yj+zk
                Ti=ray_mesh_pathlength(v,v+din,mesh)
                To=ray_mesh_pathlength(v,v+dout,mesh)
                if Ti>=0 || To>=0 # otherwise this point isn't in the crystal
                    vijk = w[i]*(xmax-xmin)/2 * w[j]*(ymax-ymin)/2 * w[k]*(zmax-zmin)/2
                    result += vijk*exp(-μ*(Ti+To))
                    volume += vijk # the volume is the integral of f(x)=1
                end
            end
        end
    end
    #return (result,volume)
    invA= result>0 ? volume/result : 0 # A = (∭e^{-μT}dV)/V. if V≡0 A=∞, which we might want to protect againts.
    return invA # since we care about 1/A, V≡0 is OK, but we need to protect against result≡0
end
absorption_correction(tθ::Real,s::GeometrySample;kwds...)=absorption_correction(tθ,s.μ,s.geometry;kwds...)

_float_type(::GT.HomogenousMesh{GT.Point{N,T}}) where {N,T} = T

function determine_second_moments(mesh::HomogenousMesh;method::Symbol=:gausslegendre,points::Integer=10)
    o=GTp(0,0,0)
    x=GTp(1,0,0)
    y=GTp(0,1,0)
    z=GTp(0,0,1)
    (p,w)=gausslegendre(points)
    result=zeros(3,3)
    volume=0.
    mesh_limits=axis_limits(mesh)
    xmin,xmax=mesh_limits[1]
    ylim=mesh_limits[2]; zlim=mesh_limits[3]
    for i=1:points
        xi=(p[i]*(xmax-xmin)+(xmax+xmin))*x/2
        ymin=line_integration_limits(o+xi,o+xi+2*ylim[1]*y,mesh)*sign(ylim[1]) # since sign(ylim) determines the search direction
        ymax=line_integration_limits(o+xi,o+xi+2*ylim[2]*y,mesh)*sign(ylim[2])
        for j=1:points
            yj= (p[j]*(ymax-ymin) + (ymax+ymin))*y/2
            zmin=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[1]*z,mesh)*sign(zlim[1])
            zmax=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[2]*z,mesh)*sign(zlim[2])
            for k=1:points
                zk=(p[k]*(zmax-zmin)+(zmax+zmin))*z/2
                pijk=o+xi+yj+zk
                vijk = w[i]*(xmax-xmin)/2 * w[j]*(ymax-ymin)/2 * w[k]*(zmax-zmin)/2
                result+= vijk*(pijk*pijk.')
                volume+= vijk
            end
        end
    end
    sm=result/volume
    sm.*(abs.(sm).>eps(_float_type(mesh))) # set absolute values smaller than 2e-16 to +/-0
end

function make_circle_z(d::Float32,N::Integer=25;center::GTp=GTp(0,0,0))
    @assert N>2
    θ=linspace(0,2π,N+1)
    r=d/2
    ring=[GTp(r*cos(x),r*sin(x),0) for x in θ[1:end-1]] # 1:end-1 to avoid the doubled point
    verts=vcat([center],[center].+ring)
    idxs=1+hcat( zeros(Int,N), 1:N, circshift(1:N,-1) ) # 1-based indexing
    return (verts,idxs)
end

function make_cylinder_z(l::Float32,d::Float32,N::Integer=25;center::GTp=GTp(0,0,0))
    tvrt,tidx=make_circle_z(d,N;center=GTp(0,0, l/2))
    bvrt,bidx=make_circle_z(d,N;center=GTp(0,0,-l/2))
    verts=vcat(tvrt,bvrt).+[center] # 2(N+1) vertices spanning 0:2N+1
    bidx = bidx[:,[1,3,2]] + N+1 # the bottom face indicies must be offset and flipped
    #bnrm *= -1 # and their normals must be reversed
    i1 = 1:N; i2 = N+2:2N+1; i3=circshift(i2,-1); i4=circshift(i1,-1);
    midx = 1+reshape(hcat( i1,i2,i3,  i1,i3,i4 ).',(3,2N)).' # 1-based indexing
    return (verts, vcat(tidx,midx,bidx))
end

function indices_to_faces(i::Array{I,2}) where I<:Integer
    (a,b)=size(i)
    [GT.Face{b,GT.OffsetInteger{-1,UInt32}}(i[j,:]...) for j=1:a]
end

function mesh_cylinder_z(l::Float32,d::Float32,N::Integer=25;center::GTp=GTp(0,0,0))
    v,i = make_cylinder_z(l,d,N;center=center)
    GT.GLNormalMesh(GT.PlainMesh(vertices=v,faces=indices_to_faces(i)))
end
#
# function triangle_normal(v::Vector,i::Vector{I}) where I<:Integer
#     @assert minimum(i)>0 && maximum(i)<=length(v)
#     @assert length(i)==3
#     n=cross( v[i[2]]-v[i[1]], v[i[3]]-v[i[1]] )
#     GTn(n./norm(n))
# end
# function triangle_normal(v::Vector,i::Array{I,2}) where I<:Integer
#     @assert minimum(i)>0 && maximum(i)<=length(v) # 1 indexing
#     (l,order)=size(i)
#     @assert order==3 # triangular faces
#     n=[ cross( v[i[j,2]]-v[i[j,1]] , v[i[j,3]]-v[i[j,1]] ) for j=1:l]
#     return GTn.(n./norm.(n))
# end
