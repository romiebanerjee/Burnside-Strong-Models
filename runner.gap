Read("./burnside.gap");
Read("./diophantine.gap");

CreateRunner := function(gens, target)
    local obj, group, Gr, orbits_lengths, solver, actions_list, 
    mcmc_perms, mcmc_coeffs, compute_hom,
    zero_delta, proposal, next_coeff, delta,
    hom_current, hom, hom_conj, coeff_init, perm_current,
    perm_init, hom_init, target_sum, points, coeff_list,
    no_of_iters, p, reset_episode, V;

    Gr := GroupWithGenerators(gens);
    orbits_lengths := SubgroupsCC(Gr)[2];
    solver := SolutionGenerator(orbits_lengths, target); 

    obj := rec(

    group := Gr,
    target_sum := target,
    denominations := orbits_lengths,
    sol := solver, 
    decay := 0.5,
    V := [[], []],
    returns := [[], []],
    x := Random(1, target),
    
    # Data storage lists
    reset_episode := function(self)
        
        self.coeff_init := self.sol.random(Random(1,1000));
        self.perm_init := Random(SymmetricGroup(self.target_sum));
        self.hom_init := GroupHomomorphismConjugate(BurnsideCoeffToGroupHom(self.group, self.coeff_init), self.perm_init);
        # self.points := [];
        self.coeff_list := [self.coeff_init];
        self.perm_list := [self.perm_init];
        # self.hom_list := [self.hom_init];
        self.actions_list := [];
        self.hom_current := self.hom_init;
        self.perm_current := self.perm_init;
        self.coeff_current := self.coeff_init;
        self.cell_complexities := [self.cell_complexity(self)];
        self.rewards := [];
        
    end,
    
    # MCMC methods
    mcmc_perms := function(self, g)
        local next_perm, zero_delta, hom_conj, cc, reward;
        # Print("\r MCMC perm ");
        # g := Random(SymmetricGroup(self.target_sum));
        next_perm := self.perm_current * g;
        zero_delta := List([1..Length(self.denominations)], i -> 0);
        hom_conj := GroupHomomorphismConjugate(self.hom_current, next_perm);

        Add(self.actions_list, [zero_delta, g]);
        Add(self.perm_list, next_perm);
        Add(self.coeff_list, self.coeff_current);
        # Add(self.hom_list, hom_conj);

        self.hom_current := hom_conj;
        self.perm_current := next_perm;
        cc := self.cell_complexity(self);
        Add(self.cell_complexities, cc);

        reward := self.reward(self);
        Add(self.rewards, reward);
    
    end,
    
    mcmc_coeffs := function(self, proposal)
        local next_coeff, delta, hom, cc, reward;
        # Print("MCMC coeff");
        # proposal := self.sol.propose_move(self.coeff_current);
        next_coeff := proposal[1];
        delta := proposal[2];
        hom := BurnsideCoeffToGroupHom(self.group, next_coeff);

        Add(self.coeff_list, next_coeff);
        Add(self.actions_list, [delta, ()]);
        Add(self.perm_list, self.perm_current);
        # Add(self.hom_list, hom);

        self.hom_current := hom;
        self.coeff_current := next_coeff;
        cc := self.cell_complexity(self);
        Add(self.cell_complexities, cc);

        reward := self.reward(self);
        Add(self.rewards, reward);
    end,


    episode := function(self, no_of_iters, bias)
        local counts, i, proposal, g;
        counts := [0,0];
        Print("Running episode ... \n");
        self.reset_episode(self);
        for i in [1..no_of_iters] do

            if Random(0, 100) < bias then #  biased coin toss
                proposal := self.sol.propose_move(self.coeff_current);
                self.mcmc_coeffs(self, proposal);
                counts[1] := counts[1] + 1;
            else
                g := Random(SymmetricGroup(self.target_sum));
                self.mcmc_perms(self, g);
                counts[2] := counts[2] + 1;
            fi;
            progress_bar_with_types(i, no_of_iters, counts);
        od;
        Print("\n");
    end,

    policy_episode := function(self, samples, no_of_iters, bias)
    local i, j , g, counts, sample_state_values, sampled_states, best_next_state, sampled_actions, sampled_actions_types,
    next_action, next_action_type, next_state, best_next_action, best_next_action_type;
        
        counts := [0,0];
        self.reset_episode(self);
  

        for i in [1..no_of_iters] do 
            sample_state_values := [];
            sampled_states := [];
            sampled_actions := [];
            sampled_actions_types := [];

            for j in [1..samples] do
                proposal := self.sol.propose_move(self.coeff_current);
                next_coeff := proposal[1];
                g := Random(SymmetricGroup(self.target_sum));

                if Random(0, 100) < bias then 
                    next_state := [next_coeff, self.perm_current];
                    next_action := proposal;
                    next_action_type := 'C';
                else
                    next_state :=  [self.coeff_current, self.perm_current * g];
                    next_action := g;
                    next_action_type := 'P';
                fi;

                Add(sampled_states, next_state);
                Add(sampled_actions, next_action);
                Add(sampled_actions_types, next_action_type);


                if next_state in self.V[1] then
                    Print("BING! \n");
                    Add(sample_state_values, self.V[2][Position(self.V[1], next_state)]);
                else
                    Add(sample_state_values, RandomFloat());
                fi;
            od;

            # Print(sampled_actions_types ,"\n");
            # Print(sample_state_values, "\n");

            best_next_action := sampled_actions[Argmax(sample_state_values)];
            best_next_action_type := sampled_actions_types[Argmax(sample_state_values)];
            
            if best_next_action_type = 'C' then
                self.mcmc_coeffs(self, best_next_action);
                counts[1] := counts[1] + 1;
            else 
                self.mcmc_perms(self, best_next_action);
                counts[2] := counts[2] + 1;
            fi;
            progress_bar_with_types(i, no_of_iters, counts);

        od;
        Print("\n");
    end,

    cell_complexity := function(self)
        local K_cell, S_cell, cell, no_of_singelton_cells;
        no_of_singelton_cells := self.coeff_current[Length(self.coeff_current)];
        cell := Orbit(self.hom_current(self.group), self.x , OnPoints);
        K_cell := LogInt(Sum(self.coeff_current) - no_of_singelton_cells + 1 ,2);
        if Size(cell) = 1 then
            S_cell := LogInt(self.coeff_current[Length(self.coeff_current)], 2);
        else
            S_cell := LogInt(Size(cell), 2);
        fi;
        return [K_cell, S_cell, K_cell + S_cell];
    end,

    reward := function(self)
        local episode_id, reward;
        episode_id := Length(self.cell_complexities);
        reward := - self.cell_complexities[episode_id][3] + self.cell_complexities[episode_id - 1][3];
        return reward;
    end,

    get_returns := function(self)
        local t, episode_length, G, state;
        G := 0;
        episode_length := Length(self.cell_complexities);
        for t in [1..episode_length -1] do 
            G := self.decay*G + self.rewards[episode_length - t];
            state := [self.coeff_list[episode_length -t],  self.perm_list[episode_length - t]];
            if  state in self.returns[1] then 
                Add(self.returns[2][Position(self.returns[1], state)], G);
            else    
                Add(self.returns[1], [self.coeff_list[episode_length - t], self.perm_list[episode_length - t]]);
                Add(self.returns[2], [G]);
            fi;
        od;
    end,

    update_V := function(self)
        local values_list, return_values;
        values_list := [];
        for return_values in self.returns[2] do 
            Add(values_list, Sum(return_values)/Length(return_values));
        od;
        self.V[1] := self.returns[1];
        self.V[2] := values_list;
    end,

    monte_carlo_es := function(self, no_of_episodes, mc_samples, episode_length, bias)
        local i;
        for i in [1..no_of_episodes] do 
        self.policy_episode(self, mc_samples, episode_length, bias);
        self.get_returns(self);
        self.update_V(self);
        PrintTo("values.bin", " V := ", self.V, ";");
        PrintTo("point.bin", "point := ", self.x , ";" );

    od;
    end,
);
return obj;
end;







