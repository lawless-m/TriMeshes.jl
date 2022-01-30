
include("Spline_exports.jl")

struct Bearing
	xys::Vector{Point2D}
	angle::Float64
	function Bearing(xy::Point2D, deg)
		rad = deg2rad(deg)
		spd = 0.01
		xys = Point2D[
			Point2D(xy.x-spd * cos(rad), yx.y-spd * sin(rad)),
			xy,
		 	Point2D(xy.x+spd * cos(rad), xy.y+spd * sin(rad))
		]
		new(xys, rad)
	end
end

xs(b::Bearing) = map(xy->xy.x for xy in xys)
ys(b::Bearing) = map(xy->xy.y for xy in xys)

struct Control
	xs::Float64
	xe::Float64
	p::Polynomial{Float64}

	function Control(b1, b2)
		p = fit(vcat(xs(b1), xs(b2)), vcat(ys(b1), ys(b2)))
		new(b1.xys[2].x, b2.xys[2].x, p)
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

