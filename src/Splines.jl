module Splines

using Mesh
using Polynomials
using Printf

export Spline, Bearing, bearing!, controlForX, yVal

struct Bearing
	xs::Vector{Float64}
	ys::Vector{Float64}
	angle::Float64
	Bearing(x, y, deg) = Bearing((x, y), deg)
	function Bearing((x, y), deg)
		xs = zeros(Float64, 3)
		ys = zeros(Float64, 3)
		rad = 2pi*deg/360
		spd = 0.01
		xs[1], ys[1] = x-spd * cos(rad), y-spd * sin(rad)
		xs[2], ys[2] = x, y
		xs[3], ys[3] = x+spd * cos(rad), y+spd * sin(rad)

		new(xs, ys, rad)
	end
end

struct Control
	xs::Float64
	xe::Float64
	p::Polynomial{Float64}

	function Control(b1, b2)
		p = fit(vcat(b1.xs, b2.xs), vcat(b1.ys, b2.ys))
		new(b1.xs[2], b2.xs[2], p)
	end
end


struct Spline
	bearings::Vector{Bearing}
	controls::Vector{Control}
	Spline(b::Bearing) = new([b], Vector{Control}())
end


function controlForX(s::Spline, x::Float64)
	if x < s.controls[1].xs
		return 1
	end
	if x >= s.controls[length(s.controls)].xs
		return length(s.controls)
	end
		
	for i in 2:length(s.controls)
		if s.controls[i-1].xs <= x < s.controls[i].xs
			return i-1
		end
	end
	return length(s.controls)
end

yVal(s::Spline, x::Float64) = s.controls[controlForX(s, x)].p(x)

#== not finished, not sure if even sensible
function yVals(s::Spline, xs::Float64, xe::Float64, steps::Int)

	xs, xe = min(xs, xe), max(xs, xe)
	cs = controlForX(s, xs)
	ce = controlForX(s, xe)


	if cs == ce
		xstep = xe-xs / steps
		return xs:xstep:xe, [s.controls[cs].p(x) for x in xs:xstep:xe]
	end

	
	stepsPerC = steps / (1+(ce - cs))
	for c in cs:ce
		
	end
end
==#

function bearing!(s::Spline, b::Bearing)
	o = s.bearings[length(s.bearings)]
	push!(s.bearings, b)
	push!(s.controls, Control(o, b))
end

### STAHP
end
