module TriMeshes

using StaticArrays
using LinearAlgebra
using Polynomials

const Vertex = SVector{3,Float64}

include("Splines.jl")

export scale, transformVerts, Edge, Face, Net, vertex!, face!, stl, wrap, flatRing!, quadRSlice, quadVSlice, outerSlice, innerSlice, apply, cylinder!, flatDisk!, tube!

struct Edge
	from::Integer
	to::Integer
end

struct Face
	AB::Edge
	BC::Edge
	CA::Edge
end
	
struct Net
	vertices::Vector{Vertex}
	vhashes::Dict{UInt64, Int64}
	faces::Vector{Face}
	edges::Vector{Edge}
	Net() = new(Vector{Vertex}(), Dict{UInt64, Int64}(), Vector{Face}(), Vector{Edge}())
end

vertex!(n::Net, x, y, z) = vertex!(n, Vertex(x, y, z))
vertex!(v::Vector{Vertex}, x, y, z) = vertex!(v, Vertex(x, y, z))
vertex!(n::Net, v::Vertex) = vertex!(n.vertices, v)

function vertex!(vs::Vector{Vertex}, v::Vertex)
	h = hash(v)
	i = get(vhashes, h, 0)
	if i == 0
		push!(vs, v)
		i = length(vs)
		n.vhashes[h] = i
	end
	i
end

scale(n::Net, x, y, z) = map!(v->Vertex(v.x*x, v.y*y, v.z*z), n.vertices, n.vertices)

face!(n::Net, v1, v2, v3) = push!(n.faces, Face(Edge(v1, v2), Edge(v2, v3), Edge(v3, v1)))

face!(f::Vector{Face}, v1, v2, v3) = push!(f, Face(Edge(v1, v2), Edge(v2, v3), Edge(v3, v1)))

apply(net, vs, fn) = foreach(v->net.vertices[v] = fn(net.vertices[v]), vs)

function wrap(n::Net, az2r)

	function trans(v::Vertex)
		a = 2pi * v.x / 360
		r = az2r(a, v.z) + v.y
		Vertex(r * cos(a), r * sin(a), v.z)
	end

	map!(trans, n.vertices, n.vertices)
end

function quadRSlice(az2r, a, astep, zmin, zmax) # (r1L, r2L, r1H, r2H)
	return [az2r(a, zmin), az2r(a+astep, zmin), az2r(a, zmax), az2r(a+astep, zmax)]
end

function quadVSlice(net::Net, a::Float64, astep, zmin, zmax, rs)

	v1L = vertex!(net, rs[1] * cos(a), rs[1] * sin(a), zmin)
	v2L = vertex!(net, rs[2] * cos(a+astep), rs[2] * sin(a+astep), zmin)
	v1H = vertex!(net, rs[3] * cos(a), rs[3] * sin(a), zmax)
	v2H = vertex!(net, rs[4] * cos(a+astep), rs[4] * sin(a+astep), zmax)

	return [v1L, v2L, v1H, v2H]
end

function outerSlice(net::Net, az2r::Function, asteps, zmin, zmax)
	vs = zeros(Int, 4asteps)
	astep = 2pi / asteps
	vi = 1
	for a in astep:astep:2pi
		vs[vi:vi+3] = quadVSlice(net, a, astep, zmin, zmax, quadRSlice(az2r, a, astep, zmin, zmax))
		face!(net, vs[vi], vs[vi+1], vs[vi+2])
		face!(net, vs[vi+3], vs[vi+2], vs[vi+1])
		vi += 4
	end
	vs
end

function quadVSlice(net::Net, o::Vertex, a::Float64, astep, zstep, rs)

	v1L = vertex!(net, o.x+rs[1] * cos(a), o.y+rs[1] * sin(a), o.z)
	v2L = vertex!(net, o.x+rs[2] * cos(a+astep), o.y+rs[2] * sin(a+astep), o.z)
	v1H = vertex!(net, o.x+rs[3] * cos(a), o.y+rs[3] * sin(a), o.z+zstep)
	v2H = vertex!(net, o.x+rs[4] * cos(a+astep), o.y+rs[4] * sin(a+astep), o.z+zstep)

	return [v1L, v2L, v1H, v2H]
end


function outerSlice(net::Net, o::Vertex, az2r::Function, asteps, zstep)
	vs = zeros(Int, 4asteps)
	astep = 2pi / asteps
	vi = 1
	for a in astep:astep:2pi
		vs[vi:vi+3] = quadVSlice(net, o, a, astep, zstep, quadRSlice(az2r, a, astep, o.z, o.z+zstep))
		face!(net, vs[vi], vs[vi+1], vs[vi+2])
		face!(net, vs[vi+3], vs[vi+2], vs[vi+1])
		vi += 4
	end
	vs
end

function quadVSlice(net, o, t, a, astep, rs)

	v1L = vertex!(net, o.x+rs[1] * cos(a), o.y+rs[1] * sin(a), o.z)
	v2L = vertex!(net, o.x+rs[2] * cos(a+astep), o.y+rs[2] * sin(a+astep), o.z)
	v1H = vertex!(net, t.x+rs[3] * cos(a), t.y+rs[3] * sin(a), t.z)
	v2H = vertex!(net, t.x+rs[4] * cos(a+astep), t.y+rs[4] * sin(a+astep), t.z)

	return [v1L, v2L, v1H, v2H]
end

function outerSlice(net, o, t, az2r, astart, asteps)
	vs = zeros(Int, 4asteps)
	astep = 2pi / asteps
	vi = 1
	for a in astep+astart:astep:2pi+astart
		vs[vi:vi+3] = quadVSlice(net, o, t, a, astep, quadRSlice(az2r, a, astep, o.z, t.z))
		face!(net, vs[vi], vs[vi+1], vs[vi+2])
		face!(net, vs[vi+3], vs[vi+2], vs[vi+1])
		vi += 4
	end
	vs
end

function innerSlice(net, az2r, astep, zmin, zmax) 
	vs = zeros(Int, 4asteps)
	astep = 2pi / asteps
	vi = 1
	for a in astep:astep:2pi
		vs[vi:vi+3] = quadVSlice(net, a, astep, zmin, zmax, quadRSlice(az2r, a, astep, zmin, zmax))
		face!(net, vs[vi+2], vs[vi+1], vs[vi])
		face!(net, vs[vi+1], vs[vi+2], vs[vi+3])
		vi += 4
	end
	vs
end

function cylinder!(net::Net, centreL, centreH, asteps, az2r)
	flat(net, centreL, az2r, v.asteps, "down")
	outerSlice(net, az2r, v.asteps, centreL.z, centreH.z)
	flat(net, centreH, az2r, v.asteps, "up")
end

function flatDisk!(net::Net, centre, az2r, astart, asteps, upordown)
	astep = 2pi/asteps
	for a in astep+astart:astep:2pi+astart
		v1 = vertex!(net, centre.x + az2r(a, centre.z) * cos(a), centre.y + az2r(a, centre.z) * sin(a), centre.z)
		v2 = vertex!(net, centre.x, centre.y, centre.z)
		v3 = vertex!(net, centre.x + az2r(a+astep, centre.z) * cos(a+astep), centre.y + az2r(a+astep, centre.z) * sin(a+astep), centre.z)
		if upordown == "up"
			face!(net, v3, v2, v1)
		else
			face!(net, v1, v2, v3)
		end
	end
end

function flatRing!(n::Net, originL, astart, asteps, router, rinner, thickness)
	originH = originL + Vertex(0,0,thickness)
	astep = 2pi / asteps
	for a in astep+astart:astep:2pi+astart

		vo = map(p->vertex!(n, p[1] + Vertex(router*cos(p[2]), router*sin(p[2]), 0)), [(originL, a), (originL, a+astep), (originH, a), (originH, a+astep)])
		vi = map(p->vertex!(n, p[1] + Vertex(rinner*cos(p[2]), rinner*sin(p[2]), 0)), [(originL, a), (originL, a+astep), (originH, a), (originH, a+astep)])

		face!(n, vo[1], vo[2], vo[3])
		face!(n, vo[4], vo[3], vo[2])
		face!(n, vi[3], vi[2], vi[1])
		face!(n, vi[2], vi[3], vi[4])

		face!(n, vo[1], vi[1], vi[2])
		face!(n, vi[2], vo[2], vo[1])
		face!(n, vo[4], vi[3], vo[3])
		face!(n, vo[4], vi[4], vi[3])

	end

end

function tube!(net, radius, sides, astart, srcnet, srcvertices, startz, endz)
	for vtxn in srcvertices
		if srcnet.vertices[vtxn[1]].z==startz
			flatDisk!(net, srcnet.vertices[vtxn[1]], (a,z)->radius, astart, sides, "down")
		end
		outerSlice(net, srcnet.vertices[vtxn[1]], srcnet.vertices[vtxn[2]], (a,z)->radius, astart, sides)
		if srcnet.vertices[vtxn[2]].z==endz
			flatDisk!(net, srcnet.vertices[vtxn[2]], (a,z)->radius, astart, sides, "up")
		end
	end

end


const Vertex = SVector{3,Float64}

export Vertex, angleXY, angleYX, angleYZ, angleZY, angleXZ, angleZX, rotate, translate, transformVerts, Edge, Face, Net, vertex!, face!, abc, abci, areaXY, magnitude


import Base.show
function show(io::IO,v::Vertex)
       print(io,"Vertex($(round(v.x, digits=4)), $(round(v.y, digits=4)), $(round(v.z, digits=4)))")
end


#import Base.LinAlg.vecdot
function vecdot(v1::Vertex, v2::Vertex)
	vecdot([v1.x, v1.y, v1.z], [v2.x, v2.y, v2.z])
end

#import Base.LinAlg.cross
function cross(v1::Vertex, v2::Vertex)
	cross([v1.x, v1.y, v1.z], [v2.x, v2.y, v2.z])
end

#import Base.LinAlg.normalize
function normalize(v::Vertex)
      n = normalize([v.x, v.y, v.z])
      Vertex(n[1], n[2], n[3])
end

function angle(x,y) 
	# atan(y,x) - take point (x,y) - angle from x axis to point
	if x == y == 0
		0.0
	else
		atan(y,x)
	end
end

function angleXY(v::Vertex)
	# project onto plane XY, angle from X axis to point		
	angle(v.x, v.y)
end

function angleYX(v::Vertex)
	# project onto plane XY, angle from Y axis to point		
	angle(v.y, v.x)
end

function angleYZ(v::Vertex)
	# project onto plane YZ, angle from Y axis to point
	angle(v.y, v.z)
end

function angleZY(v::Vertex)
	# project onto plane YZ, angle from Z axis to point
	angle(v.z, v.y)
end

function angleXZ(v::Vertex)
	# project onto plane XZ, angle from X axis to point
	angle(v.x, v.z)
end

function angleZX(v::Vertex)
	# project onto plane XZ, angle from Z axis to point
	angle(v.z, v.x)
end

function rotate(v::Vertex, a::Vertex)
	if abs(a.x) > 0
		v = Vertex(v.x, v.y * cos(a.x) - v.z * sin(a.x), v.y * sin(a.x) + v.z * cos(a.x))
	end
	if abs(a.y) > 0
		v = Vertex(v.z * sin(a.y) + v.x * cos(a.y), v.y, v.z * cos(a.y) - v.x * sin(a.y))
	end
	if abs(a.z) > 0
		v = Vertex(v.x * cos(a.z) - v.y * sin(a.z), v.x * sin(a.z) + v.y * cos(a.z), v.z)
	end
	v
end

rotate(x, y, z, xa, ya, za) = rotate(Vertex(x,y,z), Vertex(xa, ya, za))
rotate(v::Vertex, xa, ya, za) = rotate(v, Vertex(xa, ya, za))
rotate(x, y, z, v::Vertex) = rotate(Vertex(x,y,z), v)

translate(v::Vertex, t::Vertex) = v+t

vertex!(n::Net, x, y) = vertex!(n, Vertex(x, y, 0.0))

function vertex!(n::Net, v::Vertex)
	push!(n.vertices, v)
	length(n.vertices)
end

face!(n::Net, a, b, c) = face!(n, Face(Edge(a, b), Edge(b, c), Edge(c, a)))

face!(n::Net, a, b, c, ab, bc, ca) = face!(n, Face(Edge(a, b, ab), Edge(b, c, bc), Edge(c, a, ca)))

function face!(n::Net, f::Face)
	push!(n.faces, f)
	length(n.faces)
end

abci(n::Net, f::Face) = [f.ab.from, f.ab.to, f.ca.from]
abci(n::Net, i::Integer) = abci(n, n.faces[i])

abc(n::Net, f::Face) = n.vertices[abci(n, f)]
abc(n::Net, i::Integer) = abc(n, n.faces[i])

transformVerts(n::Net, f, vs::Integer, ve::Integer) = n.vertices[vs:ve] = map(f, n.vertices[vs:ve])

magnitude(v::Vertex) = sqrt(v.x^2 + v.y^2 + v.z^2)

function areaXY(n::Net, f::Face) 
	a,b,c = abc(n, f)

	ab = b-a
	ca = a-c

	h = abs(magnitude(ab) * sin(angleXY(ab) - angleXY(ca)))
	area = 0.5 * h * magnitude(ca)
end

function quad!(v1, v2, v3, v4)
	a = vertex!(n, v1)
	b = vertex!(n, v2)
	c = vertex!(n, v3)
	d = vertex!(n, v4)

	face!(n, a, b, c)
	face!(n, a, c, d)
end


function cube(w, h, d)
	n = Net()
	# front
	quad(n, Vertex(0,0,0), Vertex(0,h,0), Vertex(w,h,0), Vertex(0,w,0))
	# left
	quad(n, Vertex(0,0,d), Vertex(0,h,d), Vertex(0,h,0), Vertex(0,0,0))
	# back
	quad(n, Vertex(w,0,d), Vertex(w,h,d), Vertex(0,h,d), Vertex(0,0,d))
	# right
	quad(n, Vertex(w,0,0), Vertex(w,h,0), Vertex(w,h,d), Vertex(w,0,d))
	# top
	quad(n, Vertex(0,h,0), Vertex(0,h,d), Vertex(w,h,d), Vertex(w,h,0))
	# bottom
	quad(n, Vertex(0,0,d), Vertex(0,0,0), Vertex(w,0,d), Vertex(w,0,0))
end

###
end
