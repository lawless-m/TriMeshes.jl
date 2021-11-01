#!/usr/local/bin/julia 

# ./Mesh_test.jl

push!(LOAD_PATH, "/home/matt/Documents/julia/3dprinting/")


using Mesh
using STL

n = Net()

v1 = vertex!(n, 0.0, 0.0, 0.0)
v2 = vertex!(n, 0.0, 1.0, 0.0)
v3 = vertex!(n, 3.0, 0.0, 0.0)
v4 = vertex!(n, 0.0, 0.0, 1.0)
face!(n, v1, v2, v3)
face!(n, v3, v4, v2)
face!(n, v1, v4, v3)
face!(n, v4, v1, v2)

STL_ASCII(n, "tri.stl")

