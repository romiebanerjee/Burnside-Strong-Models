progress_bar := function(i, total)
    local percent, barLength, bar, barString, spaces;
    percent := Int( i / total * 100 );
    # Create the progress bar (20 characters wide)
    barLength := Int( (i / total) * 20 );
    bar := ListWithIdenticalEntries(barLength, '=');
    spaces := ListWithIdenticalEntries(20 - barLength, ' ');
    # Combine the bar and spaces into one string
    barString := Concatenation(bar, spaces);
    # Print the progress bar with \r to overwrite the line
    Print("\rProgress: [", barString, "] ", percent, "%  \c");
end;

progress_bar_with_types := function(i, total, counts)
    local percent, barLength, bar, barString, spaces, c , p ;
    c := counts[1];
    p := counts[2];
    percent := Int( i / total * 100 );
    # Create the progress bar (20 characters wide)
    barLength := Int( (i/total) * 20);
    bar := ListWithIdenticalEntries(barLength, '=');
    spaces := ListWithIdenticalEntries(20 - barLength, ' ');
    # Combine the bar and spaces into one string
    barString := Concatenation(bar, spaces);
    # Print the progress bar with \r to overwrite the line
    Print("\rProgress: [", barString, "] ", percent, "% ", c, "C ,", p, "P \c");
end;


SubgroupsCC := function(G)
    local c, cc, rep, reps, os, orbits_lengths;
    cc := ConjugacyClassesSubgroups(G);
        reps := [];
        orbits_lengths := [];
        for c in cc do
            rep := Representative(c);
            os := Order(G)/Order(rep);
            Add(reps, rep);
            Add(orbits_lengths, os);
        od;
        return List([reps, orbits_lengths]);
end;

Argmax := function(list)
    local max_val, max_pos, i, current;  
    if Length(list) = 0 then
        return fail;
    fi; 
    # Convert first element to float for comparison if needed
    if IsInt(list[1]) then
        max_val := Float(list[1]);
    else
        max_val := list[1];
    fi;
    max_pos := 1;  
    for i in [2..Length(list)] do
        # Convert to float for safe comparison
        if IsInt(list[i]) then
            current := Float(list[i]);
        else
            current := list[i];
        fi;       
        if current > max_val then
            max_val := current;
            max_pos := i;
        fi;
    od;
    return max_pos;
end;

RandomFloat := function()#random float between -10 and 10
    return 2.0 * Random(0, 1000000) / 100000 - 1.0;
end;
# # Works with mixed types
# mixed := [10, 3.14, 42, 2.718, 100];
# Print(ArgmaxSafe(mixed));  # returns 5 (position of 100)