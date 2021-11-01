#!/usr/local/bin/julia 

using Mesh

n = Net()

function slice!(n::Net, r1, r2, z1, z2, rstp, steps)
	r = r1
	for z in z1:z2
		vo = vertex!(n, 0.0, 0.0, z)
		for t = rstp:rstp:2pi
			vertex!(n, r * sin(t), r * cos(t), z)
		end

		face!(n, vo, vo+steps, vo+1)
		for v in vo+2:vo+steps
			face!(n, vo, v-1, v)
		end
	
		if z > z1
			face!(n, vo+1, vo+steps, vo-1)
			face!(n, vo+1, vo-steps, vo-1)
			for s in 1:steps-1
				face!(n, vo+s+1, vo+s, vo+s-steps-1)
				face!(n, vo+s-steps-1, vo+s-steps, vo+s+1 )
			end
		end
		r = r2
	end
end


function solid(n::Net) 	
	r = 30
	steps = 5
	rstp = 2pi / steps
	zstp = 21
	
	for z in 0.0:zstp:360.0
		r = 30 + 15 * sin(pi * 2z/360)
		vo = vertex!(n, 0.0, 0.0, z)
		ve = vo
		for t = rstp:rstp:2pi
			ve = vertex!(n, r * sin(t), r * cos(t), z)
		end
		if z < zstp || z > 300 - zstp
			face!(n, vo, vo+steps, vo+1)
			for v in vo+2:vo+steps
				face!(n, vo, v-1, v)
			end
		end
	
		if z > 0
			face!(n, vo+1, vo+steps, vo-1)
			face!(n, vo+1, vo-steps, vo-1)
			for s in 1:steps-1
				face!(n, vo+s+1, vo+s, vo+s-steps-1)
				face!(n, vo+s-steps-1, vo+s-steps, vo+s+1 )
			end
		end
	end
	
	STL_ASCII(n, "circ.stl")
end

function hollow(n::Net)
	r = 30
	steps = 14
	rstp = 2pi / steps
	zstp = 21
	thk = 5
	vso = veo = vsi = vei = 0
	
	
	slice!(n, 30, 30, -thk, 0, rstp, steps)

	z = 0

	for z in 0:zstp:360.0
		r = 30 + 15 * sin(pi * 2z/360)
		vso = length(n.Vertices)+1
		veo = vso
		t = 0
		while t < 2pi
			veo = vertex!(n, r * sin(t), r * cos(t), z)
			t += rstp
		end
		vsi= veo+1
		vei = vsi
		t = 0
		while t < 2pi
			vei = vertex!(n, (r-thk) * sin(t), (r-thk) * cos(t), z)
			t += rstp
		end

		if z > 0
			# outside
			for v in vso:veo-1
				face!(n, v-2steps, v, v-2steps+1)
				face!(n, v, v+1, v-2steps+1)
			end
			face!(n, veo, vso, vso- 2steps)
			face!(n, veo-2steps, veo, vso- 2steps)
	

			# inside
			for v in vsi:vei-1
				face!(n, v-2steps, v, v-2steps+1)
				face!(n, v, v+1, v-2steps+1)
			end
			face!(n, vei, vsi, vsi- 2steps)
			face!(n, vei-2steps, vei, vsi- 2steps)
		end
	end

	for k in 0:steps-2
		face!(n, vso+k, vsi+k, vsi+k+1)
		face!(n, vsi+k+1, vso+k+1, vso+k)
	end
	face!(n, vso+steps-1, vsi+steps-1, vso)
	face!(n, vso, vsi+steps-1, vsi)

	STL_ASCII(n, "circH.stl")


end



hollow(n)
