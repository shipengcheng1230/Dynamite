using Test

gmsh.clear()
gmsh.model.occ.addRectangle(0, 0, 0, 1, 1);
@addOption begin
    "Mesh.CharacteristicLengthMax", .3
    "Mesh.CharacteristicLengthMin", .3
end
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)

mi = get_mesh_information(2)
edge2ele = build_edge_element_map(Tri3, mi["eletags"], mi["edgenodes"])
ele2nbr = build_nbr_edge_list(Tri3, mi["eletag2index"], edge2ele)
ele2normal, ele2length = build_normal_length_list(Tri3, mi["nodecoords"], mi["nodetag2index"], mi["elenodes"])

@testset "Geometry relationship for Tri3" begin
    for (iele, nbrs) in enumerate(ele2nbr)
        for (iedge, nbr) in enumerate(nbrs)
            nbrindex, nbriedge = nbr
            if nbrindex == -1 # boundary edge
                continue
            end
            # neighbor normal vector
            @test ele2normal[iele][iedge] ≈ - ele2normal[nbrindex][nbriedge]
            # neighbor edge length
            @test ele2length[iele][iedge] ≈ ele2length[nbrindex][nbriedge]

            rnbrs = ele2nbr[nbrindex]
            meetneighbor = false
            for (jedge, rnbr) in enumerate(rnbrs)
                if rnbr[1] == iele
                    # neighbor edge index reciprocity
                    @test iedge == rnbr[2]
                    @test jedge == nbriedge
                    meetneighbor = true
                end
            end
            @test meetneighbor # must meet with each other
        end
    end
end
