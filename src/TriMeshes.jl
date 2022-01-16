module Mesh

using Printf
using LinearAlgebra
using Vertices

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
	faces::Vector{Face}
	edges::Vector{Edge}
	Net() = new(Vector{Vertex}(), Vector{Face}(), Vector{Edge}())
end

function vertex!(n::Net, x::Float64, y::Float64, z::Float64)
	vertex!(n, Vertex(x, y, z))
end

function vertex!(v::Vector{Vertex}, x::Float64, y::Float64, z::Float64)
	vertex!(v, Vertex(x, y, z))
end

function vertex!(n::Net, x::Real, y::Real, z::Real)
	vertex!(n, Vertex(Float64(x), Float64(y), Float64(z)))
end

function vertex!(v::Vector{Vertex}, x::Real, y::Real, z::Real)
	vertex!(v, Vertex(Float64(x), Float64(y), Float64(z)))
end

function vertex!(n::Net, v::Vertex)
	vertex!(n.vertices, v)
end

function vertex!(vs::Vector{Vertex}, v::Vertex)
	push!(vs, v)
	length(vs)
end

function scale(n::Net, x, y, z) 
	map!(v->Vertex(v.x*x, v.y*y, v.z*z), n.vertices, n.vertices)
end

function face!(n::Net, v1::Integer, v2::Integer, v3::Integer)
	push!(n.faces, Face(Edge(v1, v2), Edge(v2, v3), Edge(v3, v1)))
end

function face!(f::Vector{Face}, v1::Integer, v2::Integer, v3::Integer)
	push!(f, Face(Edge(v1, v2), Edge(v2, v3), Edge(v3, v1)))
end

function apply(net, vs, fn)
	for v in vs
		net.vertices[v] = fn(net.vertices[v])
	end
end

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

function quadVSlice(net::Net, a::Float64, astep::Float64, zmin::Real, zmax::Real, rs::Vector)

	v1L = vertex!(net, rs[1] * cos(a), rs[1] * sin(a), zmin)
	v2L = vertex!(net, rs[2] * cos(a+astep), rs[2] * sin(a+astep), zmin)
	v1H = vertex!(net, rs[3] * cos(a), rs[3] * sin(a), zmax)
	v2H = vertex!(net, rs[4] * cos(a+astep), rs[4] * sin(a+astep), zmax)

	return [v1L, v2L, v1H, v2H]
end

function outerSlice(net::Net, az2r::Function, asteps::Int64, zmin::Real, zmax::Real)
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

function quadVSlice(net::Net, o::Vertex, a::Float64, astep::Float64, zstep::Real, rs::Vector)

	v1L = vertex!(net, o.x+rs[1] * cos(a), o.y+rs[1] * sin(a), o.z)
	v2L = vertex!(net, o.x+rs[2] * cos(a+astep), o.y+rs[2] * sin(a+astep), o.z)
	v1H = vertex!(net, o.x+rs[3] * cos(a), o.y+rs[3] * sin(a), o.z+zstep)
	v2H = vertex!(net, o.x+rs[4] * cos(a+astep), o.y+rs[4] * sin(a+astep), o.z+zstep)

	return [v1L, v2L, v1H, v2H]
end


function outerSlice(net::Net, o::Vertex, az2r::Function, asteps::Int64, zstep::Real)
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

function quadVSlice(net::Net, o::Vertex, t::Vertex, a::Float64, astep::Float64, rs::Vector)

	v1L = vertex!(net, o.x+rs[1] * cos(a), o.y+rs[1] * sin(a), o.z)
	v2L = vertex!(net, o.x+rs[2] * cos(a+astep), o.y+rs[2] * sin(a+astep), o.z)
	v1H = vertex!(net, t.x+rs[3] * cos(a), t.y+rs[3] * sin(a), t.z)
	v2H = vertex!(net, t.x+rs[4] * cos(a+astep), t.y+rs[4] * sin(a+astep), t.z)

	return [v1L, v2L, v1H, v2H]
end

function outerSlice(net::Net, o::Vertex, t::Vertex, az2r::Function, astart::Float64, asteps::Int64)
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

function stl(n::Net)
	stl(n, STDOUT)
end
 
function stl(n::Net, fn::String)
	fid = open(fn, "w+")
	stl(n, fid)
	close(fid)
end

function stl(n::Net, io::IO)
	println(io, "solid Mesh.jl")
	for f in n.faces
		abf = n.vertices[f.AB.from]
		abt = n.vertices[f.AB.to]
		bct = n.vertices[f.BC.to]
		@printf(io, "facet normal 0 0 0\n\touter loop\n\t\tvertex %0.2f %0.2f %0.2f\n\t\tvertex %0.2f %0.2f %0.2f\n\t\tvertex %0.2f %0.2f %0.2f\nendloop\n", abf.x, abf.y, abf.z, abt.x, abt.y, abt.z, bct.x, bct.y, bct.z)
	end
	println(io, "endsolid Mesh.jl")
end
### STAHP

end
