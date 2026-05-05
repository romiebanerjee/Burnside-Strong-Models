# SolutionGenerator in GAP 
# This implements a Metropolis-Hastings random walk for Diophantine equation solutions

# Create a new SolutionGenerator object
SolutionGenerator := function(denominations, target_sum)
    local obj, compute_null_basis_func, compute_vector_complexity_func, 
          random_basis_coefficients_func, random_delta_func, delta_is_valid_func,
          propose_move_func, compute_acceptance_ratio_func, find_initial_solution_func,
          metropolis_step_func, random_func, k, null_basis ;

    k := Length(denominations);

    
    # Method: Compute null basis
    compute_null_basis_func := function()
        local num_basis, top, bottom, i, i2, j, basis, row;
        num_basis := k -1;
        
        # Top part: n[-1] * identity matrix
        top := [];
        for i in [1..num_basis] do
            top[i] := [];
            for i2 in [1..num_basis] do
                top[i][i2] := 0;
            od;
            top[i][i] := denominations[k];
        od;
        
        # Bottom part: -n[i] for i from 1 to num_basis
        bottom := [];
        for i in [1..num_basis] do
            bottom[i] := -denominations[i];
        od;
        
        # Combine top and bottom to create basis matrix (k x (k-1))
        basis := [];
        for i in [1..k] do
            row := [];
            for j in [1..num_basis] do
                if i <= num_basis then
                    Add(row, top[i][j]);
                else
                    Add(row, bottom[j]);
                fi;
            od;
            basis[i] := row;
        od;
        
        return basis;
    end;
    
    
    # Method: Generate random basis coefficients
    random_basis_coefficients_func := function(a)
        local num_basis, coeffs, lower, upper, i, j, dot_product ;
    
        num_basis := k - 1;
        coeffs := List([1..num_basis], x -> 0);
        
        lower := Int(-a[1] / denominations[k]);
        upper := Int(a[k] / denominations[1]);
        
        coeffs[1] := Random(lower, upper);
        
        for i in [2..num_basis] do
            lower := Int(-a[i] / denominations[k]);
            
            # Compute dot product of coeffs with self.n[1..i-1]
            dot_product := 0;
            for j in [1..i-1] do
                dot_product := dot_product + coeffs[j] * denominations[j];
            od;
            
            upper := Int((a[k] - dot_product) / denominations[1]);
            coeffs[i] := Random(lower, upper);
        od;
        
        return coeffs;
    end;
    
    # Method: Generate random delta vector
    random_delta_func := function(a)
        local coefficients, basis, delta, i, j;
        coefficients := random_basis_coefficients_func(a);
        basis := compute_null_basis_func();
        
        
        # Compute delta = basis * coefficients
        delta := List([1..k], x -> 0);
        for i in [1..k] do
            for j in [1..Length(coefficients)] do
                delta[i] := delta[i] + basis[i][j] * coefficients[j];
            od;
        od;
        
        return delta;
    end;
    
    # Method: Check if delta is valid
    delta_is_valid_func := function(delta, current_sol)
        local new_sol, i;
        new_sol := current_sol + delta;
        for i in [1..k] do
            if new_sol[i] < 0 then
                return false;
            fi;
        od;
        return true;
    end;
    
    # Method: Propose a move
    propose_move_func := function(current_sol)
        local delta, new_sol, acceptance_ratio;
        
        delta := random_delta_func(current_sol);
        
        # Check validity
        if not delta_is_valid_func(delta, current_sol) then
            return [current_sol, delta, 0.0];
        fi;
        
        new_sol := current_sol + delta;
        acceptance_ratio := compute_acceptance_ratio_func(current_sol, delta);
        
        return [new_sol, delta, acceptance_ratio];
    end;
    
    # Method: Compute acceptance ratio
    compute_acceptance_ratio_func := function(current_sol, delta)
        # For uniform sampling with symmetric proposal
        return 1.0;
    end;
    
    # Method: Find initial solution using greedy algorithm
    find_initial_solution_func := function()
        local indices, solution, remaining, idx, count, i, sorted;
        
        # Create list of indices and sort by denomination descending
        indices := [1..k];
        sorted := SortedList(indices, function(a,b) return denominations[a] > denominations[b]; end);
        
        solution := List([1..k], x -> 0);
        remaining := target_sum;
        
        for idx in sorted do
            if remaining >= denominations[idx] then
                count := QuoInt(remaining, denominations[idx]);
                solution[idx] := count;
                remaining := remaining - count * denominations[idx];
            fi;
        od;
        
        if remaining = 0 then
            return solution;
        fi;
        
        # Fallback: use first denomination
        solution := List([1..k], x -> 0);
        solution[1] := QuoInt(target_sum, denominations[1]);
        return solution;
    end;
    
    # Method: Metropolis-Hastings step
    metropolis_step_func := function(current_state)
        local proposal_results, proposed, delta, acceptance_ratio;
        
        proposal_results := propose_move_func(current_state);
        proposed := proposal_results[1];
        delta := proposal_results[2];
        acceptance_ratio := proposal_results[3];
        
        # Accept or reject
        if acceptance_ratio > 0. then
            return proposed;
        else
            return current_state;
        fi;
    end;

    random_func := function(iters)
        local initial, walk, i, proposal; 
        initial := find_initial_solution_func();
        walk := [initial];

        for i in [1..iters] do
            proposal := propose_move_func(walk[Length(walk)]);
            Add(walk, proposal[1]);
        od;
        return walk[Length(walk)];
    end;


    # Create the object as a record
    obj := rec(
        # Method assignments
        propose_move := propose_move_func,
        find_initial_solution := find_initial_solution_func,
        metropolis_step := metropolis_step_func, 
        random := random_func,
    );

    return obj;
end;

# Main function to demonstrate usage
Main := function()
    local denominations, target, i, sol, sample;
    
    denominations := [100, 55, 20, 12, 2, 2, 1, 1];
    target := 19290;
    
    Print("Creating SolutionGenerator...\n");
    sol := SolutionGenerator(denominations, target);
    
    sample := sol.random(12345);
    Print("Random Sample = ", sample, "\n");
    
    
    Print("\nDone!\n");
end;

# Run the main function
# Main();