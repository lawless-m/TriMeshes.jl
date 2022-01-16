module TriMeshes

export Vertex, angleXY, angleYX, angleYZ, angleZY, angleXZ, angleZX, rotate, translate, transformVerts, Edge, Face, Net, vertex!, face!, abc, abci, areaXY, magnitude

struct Vertex
	x
	y
	z
	Vertex(x,y,z) = new(x,y,z)
	Vertex(x,y) = new(x,y,zero(typeof(x)))
end

import Base.show
function show(io::IO,v::Vertex)
       print(io,"Vertex($(round(v.x, digits=4)), $(round(v.y, digits=4)), $(round(v.z, digits=4)))")
end

import Base.+
function +(a::Vertex, b::Vertex)
	Vertex(a.x+b.x, a.y+b.y, a.z+b.z)
end

import Base.-
function -(a::Vertex, b::Vertex)
	Vertex(a.x-b.x, a.y-b.y, a.z-b.z)
end

import Base./
function /(a::Vertex, s::Real)
	Vertex(a.x/s, a.y/s, a.z/s)
end
 
import Base.*
function *(a::Vertex, s::Real)
	Vertex(a.x*s, a.y*s, a.z*s)
end
function *(s::Real, a::Vertex)
	Vertex(a.x*s, a.y*s, a.z*s)
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

function rotate(x::Real, y::Real, z::Real, xa::Real, ya::Real, za::Real)
	rotate(Vertex(x,y,z), Vertex(xa, ya, za))
end

function rotate(v::Vertex, xa::Real, ya::Real, za::Real)
	rotate(v, Vertex(xa, ya, za))
end

function rotate(x::Real, y::Real, z::Real, v::Vertex)
	rotate(Vertex(x,y,z), v)
end

function translate(v::Vertex, t::Vertex)
	v+t
end


struct Edge
	from::Integer
	to::Integer
	style::String
	Edge(f, t) = new(f, t, "")
	Edge(f, t, s) = new(f, t, s)
end

struct Face
	ab::Edge
	bc::Edge
	ca::Edge
	fill::String
	Face(ab, bc, ca) = new(ab, bc, ca, "")
	Face(ab, bc, ca, fill) = new(ab, bc, ca, fill)
end

struct Net
	vertices::Vector{Vertex}
	faces::Vector{Face}
	edges::Vector{Edge}
	Net() = new(Vector{Vertex}(), Vector{Face}(), Vector{Edge}())
end

vertex!(n::Net, x::Float64, y::Float64, z::Float64) = vertex!(n, Vertex(x, y, z))
vertex!(n::Net, x::Float64, y::Float64) = vertex!(n, Vertex(x, y, 0.0))

vertex!(n::Net, x::Real, y::Real, z::Real) = vertex!(n, Vertex(Float64(x), Float64(y), Float64(z)))
vertex!(n::Net, x::Real, y::Real) = vertex!(n, Vertex(Float64(x), Float64(y), 0.0))

function vertex!(n::Net, v::Vertex)
	push!(n.vertices, v)
	length(n.vertices)
end

function face!(n::Net, a::Integer, b::Integer, c::Integer)
	face!(n, Face(Edge(a, b), Edge(b, c), Edge(c, a)))
end

function face!(n::Net, a::Integer, b::Integer, c::Integer, ab, bc, ca)
	face!(n, Face(Edge(a, b, ab), Edge(b, c, bc), Edge(c, a, ca)))
end

function face!(n::Net, f::Face)
	push!(n.faces, f)
	length(n.faces)
end

abci(n::Net, f::Face) = [f.ab.from, f.ab.to, f.ca.from]
abci(n::Net, i::Integer) = abci(n, n.faces[i])

abc(n::Net, f::Face) = n.vertices[abci(n, f)]
abc(n::Net, i::Integer) = abc(n, n.faces[i])

function transformVerts(n::Net, f, vs::Integer, ve::Integer)
	n.vertices[vs:ve] = map(f, n.vertices[vs:ve])
end

magnitude(v::Vertex) = sqrt(v.x^2 + v.y^2 + v.z^2)

function areaXY(n::Net, f::Face) 
	a,b,c = abc(n, f)

	ab = b-a
println(ab)
	ca = a-c

	h = abs(magnitude(ab) * sin(angleXY(ab) - angleXY(ca)))
	area = 0.5 * h * magnitude(ca)
end


### STAHP

end
