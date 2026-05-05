Read("./runner.gap");

gens := [(1,2), (2,3), (3,4), (4,5), (5,6)];
target_sum := 2^15;

Print("Inititializing runner .. \n");
runner := CreateRunner(gens,target_sum);
Print("Point = ", runner.x, "\n");
no_of_episodes := 1000 ;
mc_samples := 100 ;
episode_length := 1000 ;
bias := 5 ;
runner.monte_carlo_es(runner, no_of_episodes, mc_samples, episode_length, bias);;
