
initial_puzzle{2} = [0,9,0,0,0,0,5,0,0;
           0,0,1,8,9,0,0,2,4;
           0,0,0,0,0,0,7,0,9;
           0,0,4,0,8,2,0,0,0;
           8,0,0,0,6,0,0,0,3;
           0,0,0,3,5,0,2,0,0;
           5,0,9,0,0,0,0,0,0;
           7,4,0,0,2,5,1,0,0;
           0,0,2,0,0,0,0,7,0];
initial_puzzle{3} = [0,6,0,1,0,0,0,0,5;
                     0,0,0,2,0,0,0,0,7;
                     5,3,0,0,0,0,4,1,0;
                     0,0,0,5,0,2,0,3,0;
                     0,2,0,0,0,0,0,4,0;
                     0,4,0,6,0,9,0,0,0;
                     0,9,1,0,0,0,0,7,4;
                     7,0,0,0,0,3,0,0,0;
                     6,0,0,0,0,1,0,9,0];
initial_puzzle{1} = [5,3,0,0,7,0,0,0,0;
                     6,0,0,1,9,5,0,0,0;
                     0,9,8,0,0,0,0,6,0;
                     8,0,0,0,6,0,0,0,3;
                     4,0,0,8,0,3,0,0,1;
                     7,0,0,0,2,0,0,0,6;
                     0,6,0,0,0,0,2,8,0;
                     0,0,0,4,1,9,0,0,5;
                     0,0,0,0,8,0,0,7,9];
                   
initial_puzzle{4} = [0,2,0,1,0,0,0,4,0;
                     0,0,0,0,9,0,1,0,5;
                     0,3,0,0,0,8,0,0,0;
                     0,0,0,6,4,0,0,8,0;
                     0,7,0,0,3,0,0,0,2;
                     0,9,0,7,0,0,0,0,1;
                     0,0,7,0,0,0,0,0,0;
                     0,1,9,0,5,0,0,0,0;
                     0,0,0,0,2,4,0,0,0];
                   
initial_puzzle{5} = [6,0,0,0,0,0,3,0,0;
                     0,0,0,9,0,0,0,0,0;
                     2,0,0,4,0,8,0,1,0;
                     0,5,0,0,3,0,0,7,2;
                     0,1,0,0,0,0,0,0,0;
                     0,0,0,0,0,7,9,0,0;
                     5,0,6,0,0,9,0,0,8;
                     7,0,0,8,5,0,0,0,0;
                     0,0,0,0,2,0,4,0,0];
               
initial_puzzle{7} = [0,0,0,0,0,0,0,0,0;
                     0,0,0,0,0,3,0,8,5;
                     0,0,1,0,2,0,0,0,0;
                     0,0,0,5,0,7,0,0,0;
                     0,0,4,0,0,0,1,0,0;
                     0,9,0,0,0,0,0,0,0;
                     5,0,0,0,0,0,0,7,3;
                     0,0,2,0,1,0,0,0,0;
                     0,0,0,0,4,0,0,0,9];
                   
initial_puzzle{6} = [8,0,0,0,0,0,0,0,0;
                  0,0,3,6,0,0,0,0,0;
                  0,7,0,0,9,0,2,0,0;
                  0,5,0,0,0,7,0,0,0;
                  0,0,0,0,4,5,7,0,0;
                  0,0,0,1,0,0,0,3,0;
                  0,0,1,0,0,0,0,6,8;
                  0,0,8,5,0,0,0,1,0;
                  0,9,0,0,0,0,4,0,0];


generations = 2000;
pop_size = 100;
crossover_fraction = 0.8;
mutation_fraction = 0.01;
replacement = 'gen_elitism'; %'ST' 'RTS' 'RW' 'gen_elitism'
num_elites = 1;
radius = 20;
CF = 75;
period = 50;
amp = 99;
extra = num_elites; %[period, amp];
filename = 'sudoku_out_20180228_GE.mat';

num_trials = 3;
num_puzzles = 2;
sudoku_score = 1000*ones(num_trials,num_puzzles);
num_gens = zeros(num_trials,num_puzzles);
elapsed_time = zeros(num_trials,num_puzzles);
success_rate = 2*ones(1,num_puzzles);
avg_time = zeros(1,num_puzzles);
avg_gens = zeros(1,num_puzzles);

for puzzle_chosen = 1:num_puzzles
fprintf('\n Puzzle %d, Number of givens: %d',puzzle_chosen,nnz(initial_puzzle{puzzle_chosen}));
%print_sudoku_h(initial_puzzle{puzzle_chosen})


solution = dlxsolver(initial_puzzle{puzzle_chosen},1);
%print_sudoku_h(solution)

%{
pop_size = [20,50,100,200];
elite_fraction = [0.1,0.2,0.3];
mutation_fraction = [0.1,0.25,0.5];
diversify_gen = [20,50,100];
elapsed_time = zeros(length(pop_size),length(elite_fraction),length(mutation_fraction),length(diversify_gen));
sudoku_score = zeros(length(pop_size),length(elite_fraction),length(mutation_fraction),length(diversify_gen));
sudoku_puzzle = cell(length(pop_size),length(elite_fraction),length(mutation_fraction),length(diversify_gen));
%block_order = [9,8,7,3,6,5,2,1,4];
%}

for i = 1:num_trials
tic;
%A = solve_sudoku_brute(initial_puzzle,cycle_limit,block_order);
%A = solve_sudoku_brute_cell(initial_puzzle,cycle_limit);
%A = reshape(dlxsolver(initial_puzzle(:),1),[9,9]);
%[sudoku_puzzle,cycles,checks] = solve_sudoku_brute_rc(initial_puzzle{puzzle_chosen},row,cycle_limit);
[~,sudoku_score(i,puzzle_chosen),num_gens(i,puzzle_chosen),~,~,~,~] = solve_sudoku_ga_tests(initial_puzzle{puzzle_chosen},solution,generations,pop_size,replacement,extra);

elapsed_time(i,puzzle_chosen) = toc;
end
%print_sudoku_h(sudoku_puzzle);
%print_sudoku(reshape(A,[9,9]));
%{
figure(puzzle_chosen);
plot(entropy_hist);
title(sprintf('Entropy History for Elitism, Puzzle %d',puzzle_chosen));
ylabel('Entropy');
xlabel('Generations');
%}

success_rate(puzzle_chosen) = nnz(~sudoku_score(:,puzzle_chosen))/num_trials;
avg_time(puzzle_chosen) = sum((~sudoku_score(:,puzzle_chosen)).*elapsed_time(:,puzzle_chosen))/nnz(~sudoku_score(:,puzzle_chosen));
avg_gens(puzzle_chosen) = sum((~sudoku_score(:,puzzle_chosen)).*num_gens(:,puzzle_chosen))/nnz(~sudoku_score(:,puzzle_chosen));
fprintf('\n Success Rate: %5.3f\n', success_rate(puzzle_chosen));
fprintf('\n Average Time on Successes: %5.3f\n', avg_time(puzzle_chosen));
fprintf('\n Average Gens on Successes: %5.3f\n', avg_gens(puzzle_chosen));
save(filename,'success_rate','avg_time','avg_gens','sudoku_score','elapsed_time','num_gens');
end