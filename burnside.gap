Read("./utils.gap");

BuildMixedAction := function(G, config)
    local entry, H, mult, cosets, sz, i, offsets, coset_lists, sizes, total, action;
    offsets := [];
    coset_lists := [];
    sizes := [];
    total := 0;
    for entry in config do
        H := entry[1];
        mult := entry[2];
        cosets := RightCosets(G, H);
        sz := Size(cosets);
        Add(offsets, total);
        Add(coset_lists, cosets);
        Add(sizes, sz);
        # Print(offsets, sizes, total, "\n");
        total := total + mult * sz;
    od;
    action := function(g)
        local result, i, mult, sz, cosets, offset_val, copy, pos, new_pos;
        result := [];
        for i in [1..Length(config)] do
            # progress_bar(i, Length(config));
            mult := config[i][2];
            sz := sizes[i];
            cosets := coset_lists[i];
            offset_val := offsets[i];
            for copy in [1..mult] do
                for pos in [1..sz] do
                    new_pos := Position(cosets, cosets[pos] * g);
                    Add(result, offset_val + (copy - 1) * sz + new_pos);
                od;
            od;
        od;
        # Print("\n");
        return PermList(result);
    end;
    return action;
end;



BurnsideCoeffToGroupHom := function(G, coeffs)
    local  gens, cc, reps, orbits_lengths, target_size, config, action, image_perms, hom;
    gens := GeneratorsOfGroup(G);
    cc := SubgroupsCC(G);
    reps := cc[1];
    orbits_lengths := cc[2];

    target_size := Sum(List([1..Length(reps)], i -> orbits_lengths[i]*coeffs[i]));
    # Print("Target = ", target_size, "\n");
    config :=  List([1..Minimum(Length(reps), Length(coeffs))], i -> [reps[i], coeffs[i]]);
    action := BuildMixedAction(G, config);

    image_perms := List([1..Length(gens)] , i -> action(gens[i]));
    #Complete permutation repesentations
    hom := GroupHomomorphismByImages(G, SymmetricGroup(target_size),
                                        gens, image_perms);
    return hom;
end;


GroupHomomorphismConjugate := function(hom, g)
    local gens, perms, perms_conj, conj_hom;
    gens := GeneratorsOfGroup(Source(hom));
    perms := List(gens, g -> hom(g));
    perms_conj := List(perms, p -> g*p*g^-1);
    conj_hom := GroupHomomorphismByImages(Source(hom), Range(hom), gens, perms_conj);
    return conj_hom;
end;



test := function(gens, x)
    local G, reps, burns_coeffs, hom, orbit_x, g, conj_hom;

    G := GroupWithGenerators(gens);

    reps := SubgroupsCC(G)[1];

    burns_coeffs := List([1..Length(reps)], i -> Random([1..10]));
    Print(burns_coeffs, "\n");

    hom := BurnsideCoeffToGroupHom(G, burns_coeffs);

    #Calculate single Orbit from x

    orbit_x := Orbit(hom(G), x, OnPoints);

    Print("Orbit of ", x, " = ", orbit_x , "\n");
    #Apply Conjugation to hom

    g := Random(Range(hom));
    conj_hom := GroupHomomorphismConjugate(hom, g);
    orbit_x := Orbit(conj_hom(G), x, OnPoints);

    Print("Conjugate Orbit of ", x, " = ", orbit_x, "\n" );

end;

# gens := [(1,2), (2,3), (3,4), (4,5), (5,6), (6,7)];
# x := 1008;
# test(gens, x);
