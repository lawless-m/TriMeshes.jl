#!/usr/local/bin/julia 

using TriMeshes

debug = false

n = Net()


function profile(path::Vertex)
	[(1,2), (2,1.6), (3,1.75), (3.5, 2), (4,2.2), (5,2), (6,1), (7,0.2), (8, 0.2), (9, 0.8), (9.3, 1), (10, 1.6), (10.3, 2), (10.7, 3), (10.8, 4), (10.9,5), (10.7, 6), (10.5, 7), (10,7.2), (9,7.3), (8, 7.1), (7.8, 7), (7.3, 6), (7,5.8), (6.4, 5), (6,4.75), (5,4), (4,4.5), (3.5, 5), (3.2, 6), (3,6.3), (2.3,7), (2,7.01), (1,6), (0.7, 5), (0.5,4), (0.6, 3)], Vertex(0,0,1)
end

function arc(t::Real)
	Vertex(30 * cos(2pi*t), 0, 30 * sin(2pi*t))
end

function normal(t::Real, tstep::Real, fpath)
	t0 = fpath(t - tstep, tstep)
	t1 = fpath(t, tstep)
	t2 = fpath(t + tstep, tstep)
	normalize((t1 - t0) + (t2 - t1))
end

function solid(origin::Vertex, n::Net, slices::Real, fpath, slicer)
	path = Vertex(0,0,0) # to return end point
	tstep = 1/(slices-1)
	ve = 0
	for t in 0.0:tstep:1.0

		if debug @printf("\nT: %0.2f\n", t) end

		path = origin + fpath(t, tstep)
		pathv = vertex!(n, path)
		slicepts, slicenorm = slicer(path) 
		steps = length(slicepts)
		pathnorm = normal(t, tstep, fpath)


		if debug @printf("pathnorm %s\n", pathnorm) end
		if debug @printf("slicenorm %s\n", slicenorm) end

		pathay = angleXZ(pathnorm)
		sliceay = angleXZ(slicenorm)
		if sliceay < 0
			sliceay += 2pi
		end

		pathax = angleYZ(pathnorm)
		sliceax = angleYZ(slicenorm)

		rot = Vertex(0., sliceay-pathay, 0.)

		if debug @printf("pathax:%d sliceax:%d pathay:%d sliceay:%d\n", rad2deg(pathax), rad2deg(sliceax), rad2deg(pathay), rad2deg(sliceay)) end
		if debug @printf("rotating x:%d y:%d z:%d\n", rad2deg(rot.x), rad2deg(rot.y), rad2deg(rot.z)) end

		for (x,y) in slicepts
			v = rotate(Vertex(x, y, 0), rot)
			ve = vertex!(n, v + path)
		end


if t < tstep || t >= 1.0
	face!(n, pathv, pathv+steps, pathv+1)
	for v in pathv+2:pathv+steps
		face!(n, pathv, v-1, v)
	end
end

		if t > 0
			face!(n, pathv+1, pathv+steps, pathv-1)
			face!(n, pathv+1, pathv-steps, pathv-1)
			for s in 1:steps-1
				face!(n, pathv+1+s, pathv+s, pathv-1+s-steps)
				face!(n, pathv-1+s-steps, pathv+s-steps, pathv+1+s )
			end
		end
	end
	path, ve
end	

function sqr(path::Vertex)
	[(-1,-1), (1,-1), (1,1), (-1,1), ], Vertex(0,0,1)
end

function circ(path::Vertex)
	pts = Vector{Tuple{Real, Real}}()
	r = 1.5
	steps = 8
	for a = 2pi/steps:2pi/steps:2pi
		push!(pts, (r * sin(a), r * cos(a)))
	end
	pts, Vertex(0,0,1)
end



function stockp(t::Float64, tstep::Float64)
	r = 10
	x = y = z = 0.0
	function q1(u)
		r+r*cos(pi - pi*u), r*sin(pi - pi*u)
	end
	function q2(u)
		x1,z1 = q1(1.)
		x1 - 0.5 * r * u, z1 -r * u
	end
	function q3(u)
		x2,z2 = q2(1.)
		x2 + r + r*cos(pi - pi*u), z2 + r*-sin(pi - pi*u)
	end
	function q4(u)
		x3, z3 = q3(1.)
		x3 - 0.5r*u, z3 + r*u
	end

	q4m = 1.00
	q3m = q4m - tstep
	q2m = 0.50
	q1m = q2m - tstep

	if t < 0
		x,z = q1((1/q1m)*abs(t))
		z = -z
	elseif t > 1
		x,z = q4((1/(q4m-q3m))*((2-t)-q3m))
		z = -z
	elseif t < q1m
		x,z = q1((1/q1m)t)
	elseif t < q2m
		x,z = q2((1/(q2m-q1m))*(t-q1m))
	elseif t < q3m
		x,z = q3((1/(q3m-q2m))*(t-q2m))
	else
		x,z = q4((1/(q4m-q3m))*(t-q3m))
	end
	z = -z

	Vertex(x,y,z)
end

using SVG
function pathsvg(fname, tstep, fpath)
	s = open(fname, "w+")
	pts = Vector{Tuple{Real,Real}}()
	st = fpath(1., tstep)
	for stch = 1:1
		for t=-tstep:tstep:1+tstep
			v = fpath(t, tstep)
			push!(pts, (stch*st.x+v.x, stch*st.z+v.z))
		end
	end
	SVG.open(s, 500.,500.)
	SVG.polyline(s, pts, SVG.blackline(3.))
	SVG.close(s)
end


pathsvg("stk.svg", 1/38, stockp)


function sweep(n::Net, steps, multiples=1)
	e = o = Vertex(0,0,0)
	
	for k in 1:multiples
		o, e = solid(o, n, steps, stockp, circ)
	end
	o, e
end

function knit(n, stitches::Integer, rowpairs::Integer)
	ty = 15
	vs = 1
	for rp = 1:rowpairs
		pe, ve = sweep(n, 16, stitches)
		transformVerts(n, v->rotate(v, Vertex(deg2rad(40), 0, 0)), vs, ve)
		transformVerts(n, v->translate(v, Vertex(0,ty,0)), vs, ve)

		vs = ve + 1
		pe, ve = sweep(n, 16, stitches)
		transformVerts(n, v->rotate(v, Vertex(deg2rad(-40), 0, 0)), vs, ve)
		transformVerts(n, v->translate(v, Vertex(15,ty,0)), vs, ve)
		vs = ve + 1
		ty += 15
	end

	STL_ASCII(n, "stk.stl")
end


knit(n, 3, 2)

