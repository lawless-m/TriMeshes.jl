using TriMeshes
using Test

@testset "TriMeshes.jl" begin
    @test areaXY(Vertex(0,0,0), Vertex(0,1,0), Vertex(1,0,0)) == 0.5
    @test magnitude(Vertex(0,1,0)) == magnitude(Vertex(1,0,0))
end
