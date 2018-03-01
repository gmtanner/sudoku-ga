function [final_puzzle,score,current_gen,best_fitness,mean_fitness,median_fitness,entropy_hist] = solve_sudoku_ga_tests(initial_puzzle,solution,generations,pop_size,replacement,varargin)
% Solves 9x9 soduko puzzle by genetic algorithm
% input: 
%   initial_puzzle := 9x9 array with 0 representing empty cells and 1-9 for 
%              filled in cells
%   solution := 9x9 array with solution to sudoku puzzle
%   replacement := replacement strategy
%            gen_elitism: elitism with varargin{1} elites
%            multi_dyn_e: multi_dyn with elitism (always select best nondominated)
%                         varargin{1} = radius of domination
%            multi_dyn_r: multi_dyn with random choice after first choice
%                         varargin{1} = radius of domination
%            RW: Replace-worst
%            RTS: Restricted Tournament Selection Crowding
%                 varargin{1} = radius of comparison
%            ST: Saw-tooth GA, periodic re-initialization
% output:
%   final_puzzle := best solution found to sudoku puzzle

%% Parameters

%num_elites = 0;
today = 20180228;
num_parents = 2;
tourney_size = 2;
crossover_percent = 0.8;
mutation_fraction = 0.01;
elitist = false;
entropy_hist = zeros(generations+1,1);
best_fitness = zeros(generations+1,1);
mean_fitness = zeros(generations+1,1);
median_fitness = zeros(generations+1,1);
if(num_parents*tourney_size>pop_size)
  disp('Error: num_parents*tourney_size>pop_size');
  return;
end
pop_local_search = true;
max_swaps = 5;
if(strcmp(replacement,'multi_dyn_e')||strcmp(replacement,'multi_dyn_r'))
  D_init = varargin{1};
  num_children = crossover_percent*pop_size;
elseif(strcmp(replacement,'gen_elitism'))
  num_elites = varargin{1};
  elitist = true;
  num_children = pop_size - num_elites;
elseif(strcmp(replacement,'RW'))
  num_children = crossover_percent*pop_size;
elseif(strcmp(replacement,'RTS'))
  CF = varargin{1};
  num_children = crossover_percent*pop_size;
elseif(strcmp(replacement,'ST'))
  period = varargin{1}(1);
  amp = varargin{1}(2);
  num_elites = 1;
  elitist = true;
  avg_pop_size = pop_size;
  num_children = avg_pop_size+amp-num_elites;
  max_pop_size = avg_pop_size + amp;
  pop_size = max_pop_size;
end

% replacement = 'multi_dyn_e';
% multi_dyn_e: multi_dyn with elitism (always select best nondominated)
% multi_dyn_r: multi_dyn with random choice after first choice
%num_children = pop_size * crossover_fraction;
%num_mutations = ceil(mutation_fraction * nnz(~initial_puzzle));
%D_init = 30;

% debugging
filename = sprintf(strcat('sudoku_out_d%d_',replacement,'_p%d.txt'),today,pop_size);
%fileID = fopen(filename,'wt');
debug_level = 0;
gen_out = false;
%gen_out_freq = 25;


% convert initial_puzzle matrix to block form: B(:,:,i,j) is the (i,j) block 
initial_B = zeros(3,3,3,3);
solution_B = zeros(3,3,3,3);
for i = 1:3
  for j = 1:3
    initial_B(:,:,i,j) = initial_puzzle((i-1)*3+1:3*i,(j-1)*3+1:3*j);
    solution_B(:,:,i,j) = solution((i-1)*3+1:3*i,(j-1)*3+1:3*j);
  end %j
end %i

% get location of blanks and possible values
[indices,values,perm_len,givens] = initialize_sudoku(initial_B);
%immigration_chart = zeros(3,3,3,3,9);
if(~check_rows_sudoku(initial_B))
  disp('Error: initial_puzzle sudoku does not pass row check\n');
  return;
end
if(~check_cols_sudoku(initial_B))
  disp('Error: initial_puzzle sudoku does not pass column check\n');
  return;
end

% generate initial population
pop = cell(pop_size,1);
pop_score = zeros(pop_size,1);
pop_dist = zeros(pop_size,1);
for i = 1:pop_size
  pop{i} = initial_B;
  for j = 1:3
    for k = 1:3
      block = pop{i}(:,:,j,k);
      block(indices{j,k}) = values{j,k}(randperm(perm_len(j,k)));
      pop{i}(:,:,j,k) = block;
    end
  end
  if(pop_local_search)
    [pop{i},pop_score(i)] = local_search(pop{i},indices,perm_len,initial_B,givens,max_swaps);
  else
    pop_score(i) = score_sudoku(pop{i},initial_B,givens);
  end
  pop_dist(i) = nnz(pop{i} ~= solution_B);
end
best_fitness(1) = min(pop_score);
mean_fitness(1) = mean(pop_score);
median_fitness(1) = median(pop_score);

%Compute Entropy
entropy = zeros(3,3);
entropy_chart = zeros(3,3,3,3,9);
for i = 1:pop_size
    for j = 1:9
        entropy_chart(:,:,:,:,j) = entropy_chart(:,:,:,:,j) + (pop{i}==j);
    end
end
entropy_chart = entropy_chart/pop_size;
for i = 1:3
    for j = 1:3
        log_prob = log(entropy_chart(:,:,i,j,:))/log(perm_len(i,j));
        %log_prop(~isfinite(log_prob)) = 0;
        entropy(i,j) = sum(sum(sum(-1*entropy_chart(:,:,i,j,:).*log_prob,'omitnan')))/perm_len(i,j);
    end
end
entropy_hist(1) = sum(sum(entropy))/9;

% Initial Population output
if(gen_out)
    filename = sprintf(strcat('sudoku_out_',replacement,'%d_g%d_p%d.txt'),today,1,pop_size);
    fileID_gen = fopen(filename,'wt');
    fprintf(fileID_gen,'%5.3f\n',sum(sum(entropy)));
    for i = 1:pop_size
        fprintf(fileID_gen,'%d,%d',pop_score(i),pop_dist(i));
        fprintf(fileID_gen,',%d',pop{i});
        fprintf(fileID_gen,'\n');
    end
    fclose(fileID_gen);
    
    filename = sprintf(strcat('sudoku_out_dist_',replacement,'%d_g%d_p%d.txt'),today,1,pop_size);
    fileID_dist = fopen(filename,'wt');
    for i = 1:pop_size
      fprintf(fileID_dist,'%d',sum(sum(sum(sum(pop{i}==pop{1})))));
      for j = 2:pop_size
        fprintf(fileID_dist,',%d',sum(sum(sum(sum(pop{i}==pop{j})))));
      end
      fprintf(fileID_dist,'\n');
    end
    fclose(fileID_dist);
end

current_gen = 1;
if(debug_level>1)
  fprintf(fileID,'Generation %d: ',current_gen);
  fprintf(fileID,'%d,',pop_score);
  fprintf(fileID,'Entropy: ');
  fprintf(fileID,'%5.3f, ',entropy);
  fprintf(fileID,'\n');
end

offspring = cell(num_children,1);
next_gen = cell(pop_size,1);
next_gen_ind = zeros(pop_size,1);
next_gen_score = zeros(pop_size,1);
parents = zeros(num_parents,1);
child = zeros(3,3,3,3);
offspring_score = zeros(num_children,1);
% Evolution Process
while((current_gen<=generations)&&(min(pop_score)>0))
  
  
  % Pass elites to next generation
  if(elitist)
    [~,sorted_order] = sort(pop_score);
    next_gen(1:num_elites) = pop(sorted_order(1:num_elites));
    next_gen_score(1:num_elites) = pop_score(sorted_order(1:num_elites));
  end
  
  if(strcmp(replacement,'ST'))
    pop_size_new = floor(max_pop_size - 2*amp/(period-1)*(current_gen-1-period*floor((current_gen-1)/period)));
    num_children = pop_size_new - num_elites;
  end
  % Crossover
  for i = 1:num_children
    % Tournament selection of parents
    candidates = reshape(randperm(pop_size,num_parents*tourney_size),tourney_size,num_parents);
    for j = 1:num_parents
      [~,inds] = min(pop_score(candidates(:,j)));
      parents(j) = candidates(inds,j);
    end
    for j = 1:3 % iterate through blocks
      for k = 1:3
        child(:,:,j,k) = pop{parents(randi(num_parents))}(:,:,j,k);
      end 
    end
    % Mutation
    [mj,mk] = find(rand(3)<mutation_fraction*perm_len);
    for l = 1:size(mj,1)
        j = mj(l);
        k = mk(l);
        block = child(:,:,j,k);
        swap = randperm(perm_len(j,k),2);
        block(indices{j,k}(swap)) = block(indices{j,k}([swap(2),swap(1)]));
        child(:,:,j,k) = block;
    end
    % Local Search
    if(pop_local_search)
        [offspring{i},offspring_score(i)] = local_search(child,indices,perm_len,initial_B,givens,max_swaps);
    else
        offspring{i} = child;
        offspring_score(i) = score_sudoku(child,initial_B,givens);
    end
  end
  
  % Replacement Strategy
  if(strcmp(replacement,'multi_dyn_r')) 
      dcn = 1000*ones(pop_size+num_children,1);
      pool = [pop;offspring];
      pool_score = [pop_score;offspring_score];
      [min_score,min_index] = min(pool_score);
      next_gen{1} = pool{min_index};
      next_gen_ind(1) = min_index;
      next_gen_score(1) = min_score;

      dominated = false;
      D = D_init*(1-current_gen/generations);
      for i = 2:pop_size
         if ~dominated
             for k = 1:pop_size+num_children
               dcn(k) = min(dcn(k),nnz(pool{k} ~= next_gen{i-1}));
             end
             ND = find(dcn>D); % non-dominated individuals
             if ~isempty(ND)
                 next_gen_ind(i) = ND(randi(length(ND)));
                 next_gen{i} = pool{next_gen_ind(i)};
                 next_gen_score(i) = pool_score(next_gen_ind(i));
             else % if all individuals are dominated, choose a random individual from remaining pool
                 dominated = true;
                 ND = setdiff(1:pop_size+num_children,next_gen_ind);
                 next_gen_ind(i) = ND(randi(length(ND)));
                 next_gen{i} = pool{next_gen_ind(i)};
                 next_gen_score(i) = pool_score(next_gen_ind(i));
             end
         else % if all individuals are dominated, choose a random individual from remaining pool
             ND = setdiff(1:pop_size+num_children,next_gen_ind);
             next_gen_ind(i) = ND(randi(length(ND)));
             next_gen{i} = pool{next_gen_ind(i)};
             next_gen_score(i) = pool_score(next_gen_ind(i));
         end
      end
  elseif(strcmp(replacement,'multi_dyn_e'))
      dcn = 1000*ones(pop_size+num_children,1);
      pool = [pop;offspring];
      pool_score = [pop_score;offspring_score];
      [sorted_score,sorted_index] = sort(pool_score);
      next_gen{1} = pool{sorted_index(1)};
      next_gen_ind(1) = sorted_index(1);
      next_gen_score(1) = sorted_score(1);

      dominated = false;
      D = D_init*(1-current_gen/generations);
      for i = 2:pop_size
         if ~dominated
             for k = 1:pop_size+num_children
               dcn(k) = min(dcn(k),nnz(pool{sorted_index(k)} ~= next_gen{i-1}));
             end
             ND = find(dcn>D,1); % non-dominated individuals
             if ~isempty(ND)
                 next_gen_ind(i) = sorted_index(ND);
                 next_gen{i} = pool{next_gen_ind(i)};
                 next_gen_score(i) = pool_score(next_gen_ind(i));
             else % if all individuals are dominated, choose a best individual from remaining pool
                 dominated = true;
                 ND = setdiff(1:pop_size+num_children,next_gen_ind);
                 next_gen_ind(i) = ND(randi(length(ND)));
                 next_gen{i} = pool{next_gen_ind(i)};
                 next_gen_score(i) = pool_score(next_gen_ind(i));
             end
         else % if all individuals are dominated, choose a random individual from remaining pool
             ND = setdiff(1:pop_size+num_children,next_gen_ind);
             next_gen_ind(i) = ND(randi(length(ND)));
             next_gen{i} = pool{next_gen_ind(i)};
             next_gen_score(i) = pool_score(next_gen_ind(i));
         end
      end
  elseif(strcmp(replacement,'gen_elitism'))
    next_gen(num_elites+1:num_elites+num_children) = offspring;
    next_gen_score(num_elites+1:num_elites+num_children) = offspring_score;
  elseif(strcmp(replacement,'RW'))
    [~,sorted_order] = sort([pop_score;offspring_score]);
    for i = 1:pop_size
      if (sorted_order(i) > pop_size)
        next_gen{i} = offspring{sorted_order(i)-pop_size};
        next_gen_score(i) = offspring_score(sorted_order(i)-pop_size);
      else
        next_gen{i} = pop{sorted_order(i)};
        next_gen_score(i) = pop_score(sorted_order(i));
      end
    end
  elseif(strcmp(replacement,'ST'))
    if(mod(current_gen,period)==0)
      pop_size = max_pop_size;
      for i = (num_elites+1):pop_size
        next_gen{i} = initial_B;
        for j = 1:3
          for k = 1:3
            block = next_gen{i}(:,:,j,k);
            block(indices{j,k}) = values{j,k}(randperm(perm_len(j,k)));
            next_gen{i}(:,:,j,k) = block;
          end
        end
        if(pop_local_search)
          [next_gen{i},next_gen_score(i)] = local_search(next_gen{i},indices,perm_len,initial_B,givens,max_swaps);
        else
          next_gen_score(i) = score_sudoku(next_gen{i},initial_B,givens);
        end
      end
    else
      pop_size = pop_size_new;
      next_gen(num_elites+1:pop_size) = offspring(1:num_children);
      next_gen_score(num_elites+1:pop_size) = offspring_score(1:num_children);
    end
  elseif(strcmp(replacement,'RTS'))
    next_gen = pop;
    next_gen_score = pop_score;
    CF_dist = zeros(CF,1);
    for i = 1:num_children
      sample = randsample(pop_size,CF);
      for j = 1:CF
        CF_dist(j) = nnz(offspring{i} ~= next_gen{sample(j)});
      end
      [~,min_ind] = min(CF_dist);
      if(offspring_score(i) < next_gen_score(sample(min_ind)))
        next_gen{sample(min_ind)} = offspring{i};
      end
    end
  else
    disp('Warning: replacment strategy invalid');
  end  
  
  pop = next_gen;
  pop_score = next_gen_score;
  current_gen = current_gen+1;
  for k = 1:pop_size
      pop_dist(k) = nnz(pop{k} ~= solution_B);
  end
  

  %Compute Entropy
  entropy_chart = zeros(3,3,3,3,9);
  for i = 1:pop_size
    for j = 1:9
        entropy_chart(:,:,:,:,j) = entropy_chart(:,:,:,:,j) + (pop{i}==j);
    end
  end
  entropy_chart = entropy_chart/pop_size;
  entropy = zeros(3,3);
  for i = 1:3
    for j = 1:3
      log_prob = log(entropy_chart(:,:,i,j,:))/log(perm_len(i,j));
      %log_prop(~isfinite(log_prob)) = 0;
      entropy(i,j) = sum(sum(sum(-1*entropy_chart(:,:,i,j,:).*log_prob,'omitnan')))/perm_len(i,j);
    end
  end
  
  entropy_hist(current_gen) = sum(sum(entropy))/9;
  best_fitness(current_gen) = min(pop_score);
  mean_fitness(current_gen) = mean(pop_score);
  median_fitness(current_gen) = median(pop_score);
  if(debug_level>1)
    fprintf(fileID,'Generation %d: ',current_gen);
    fprintf(fileID,'%d,',pop_score);
    fprintf(fileID,'Entropy: ');
    fprintf(fileID,'%5.3f, ',entropy);
    fprintf(fileID,'\n');
    fprintf(fileID,'%d,',pop_dist);
    fprintf(fileID,'\n');
  end
end


% convert B to 9x9 form
[score,min_ind] = min(pop_score);
%pop{min_ind} = local_search(pop{min_ind},indices,perm_len,initial_B,givens,max_swaps);
for i = 1:3
  for j = 1:3
    final_puzzle((i-1)*3+1:3*i,(j-1)*3+1:3*j) = pop{min_ind}(:,:,i,j);
  end %j
end %i


%fclose(fileID);
end %solve_sudoku_brute

function [indices,values,perm_len,givens] = initialize_sudoku(B)
indices = cell(3,3);
values = cell(3,3);
givens = cell(3,3,2); 
perm_len = zeros(3,3);
for i = 1:3
  for j = 1:3
    block = B(:,:,i,j);
    values{i,j} = setdiff(1:9,nonzeros(block)); %missing values in block(i,j)
    indices{i,j} = find(~block); %indices of missing values in block (i,j)
    perm_len(i,j) = length(values{i,j}); %number of missingvalues in block(i,j)
    givens{i,j,1} = nonzeros(B(i,:,j,:)); %givens{i,j,1} is a vector of givens in row(i,j)
    givens{i,j,2} = nonzeros(B(:,i,:,j)); %givens{i,j,2} is a vector of givens in col(i,j)
    if(length(indices{i,j})~=length(values{i,j})) 
      disp('Error initialize_sudoku: length(ind)~=length(vals)\n');
    end %if
  end %j
end %i
end %initialize_sudoku

function [best_candidate,current_score] = local_search(B,indices,perm_len,initial_B,givens,max_swaps)
current_score = score_comp_sudoku(B,initial_B,givens);
best_candidate = B;
%test_vals = zeros(max_swaps+1,1);
%test_vals(1) = sum(sum(sum(current_score)));
for index = 1:max_swaps
    i = randi(3);
    j = randi(3);
    candidate = best_candidate;
    swap = indices{i,j}(randperm(perm_len(i,j),2));
    [m,n] = ind2sub([3,3],swap);
    spec_score = current_score(m(1),i,1) + current_score(m(2),i,1) + current_score(n(1),j,2) + current_score(n(2),j,2);
    if(spec_score>0)
        block = best_candidate(:,:,i,j);
        block(swap) = block(swap([2,1]));
        candidate(:,:,i,j) = block;
        candidate_score = score_spec_sudoku(candidate,initial_B,givens,m,n,i,j);
        if(sum(candidate_score)<spec_score)
          best_candidate = candidate;
          current_score(m(1),i,1) = candidate_score(1);
          current_score(m(2),i,1) = candidate_score(2);
          current_score(n(1),j,2) = candidate_score(3);
          current_score(n(2),j,2) = candidate_score(4);
        end
    end
    %test_vals(index+1) = sum(sum(sum(current_score)));
end
current_score = sum(sum(sum(current_score)));
end 

%{
function best_candidate = local_search(B,indices)
current_score = score_sudoku(B);
best_candidate = B;
for i = 1:3
  for j = 1:3
    candidate = best_candidate;
    swaps = nchoosek(indices{i,j},2);
    for k = 1:length(swaps)
      block = best_candidate(:,:,i,j);
      block(swaps(k,[1,2])) = block(swaps(k,[2,1]));
      candidate(:,:,i,j) = block;
      candidate_score = score_sudoku(candidate);
      if(candidate_score<current_score)
        current_score = candidate_score;
        best_candidate = candidate;
      end
    end 
  end 
end
end
%}

function test_flag = check_rows_sudoku(B)
test_flag = true;
for i = 1:3
  for j = 1:3
    for k = 1:9
      if(nnz(B(i,:,j,:)==k)>1)
        test_flag = false;
        return
      end %if
    end %k
  end %j
end %i
end %check_rows_sudoku

function test_flag = check_cols_sudoku(B)
test_flag = true;
for i = 1:3
  for j = 1:3
    for k = 1:9
      if(nnz(B(:,i,:,j)==k)>1)
        test_flag = false;
        return
      end %if
    end %k
  end %j
end %i
end %check_cols_sudoku

function test_flag = check_sudoku(B)
test_flag = true;
for i = 1:3
  for j = 1:3
    for k = 1:9
      if(nnz(B(i,:,j,:)==k)>1)
        test_flag = false;
        return
      elseif(nnz(B(:,i,:,j)==k)>1)
        test_flag = false;
        return
      end %if
    end %k
  end %j
end %i
end %check_sudoku

function score = score_sudoku(B,initial_B,givens)
score = 0;
alt_B = B-initial_B;
alt_score = 0;
for i = 1:3
  for j = 1:3
    score = score + score_line(B(i,:,j,:)); %9 - length(unique(B(i,:,j,:)));
    score = score + score_line(B(:,i,:,j)); %9 - length(unique(B(:,i,:,j)));
    alt_score = alt_score + sum(ismember(givens{i,j,1},alt_B(i,:,j,:))); %length(intersect(alt_B(i,:,j,:),initial_B(i,:,j,:))) - 1;
    alt_score = alt_score + sum(ismember(givens{i,j,2},alt_B(:,i,:,j))); %length(intersect(alt_B(:,i,:,j),initial_B(:,i,:,j))) - 1;
  end %j
end %i
score = score + 100*alt_score;
end %score_sudoku



function score_comp = score_comp_sudoku(B,initial_B,givens)
score_comp = zeros(3,3,2);
alt_B = B-initial_B;
for i = 1:3
  for j = 1:3
    alt_score = sum(ismember(givens{i,j,1},alt_B(i,:,j,:))); %length(intersect(alt_B(i,:,j,:),initial_B(i,:,j,:))) - 1;
    score_comp(i,j,1) = score_line(B(i,:,j,:)) + 100*alt_score; %9 - length(unique(B(i,:,j,:))) + 100*alt_score;
    alt_score = sum(ismember(givens{i,j,2},alt_B(:,i,:,j))); %length(intersect(alt_B(:,i,:,j),initial_B(:,i,:,j))) - 1;
    score_comp(i,j,2) = score_line(B(:,i,:,j)) + 100*alt_score; %9 - length(unique(B(:,i,:,j))) + 100*alt_score;
  end %j
end %i
end %score_comp_sudoku

function score = score_spec_sudoku(B,initial_B,givens,i,j,k,l)
alt_B = B - initial_B;
score = zeros(4,1);
alt_score = sum(ismember(givens{i(1),k,1},alt_B(i(1),:,k,:))); %length(intersect(alt_B(i(1),:,k,:),initial_B(i(1),:,k,:))) - 1;
score(1) = score_line(B(i(1),:,k,:)) + 100*alt_score; %9 - length(unique(B(i(1),:,k,:))) + 100*alt_score;
alt_score = sum(ismember(givens{i(2),k,1},alt_B(i(2),:,k,:))); %length(intersect(alt_B(i(2),:,k,:),initial_B(i(2),:,k,:))) - 1;
score(2) = score_line(B(i(2),:,k,:)) + 100*alt_score; %9 - length(unique(B(i(2),:,k,:))) + 100*alt_score;
alt_score = sum(ismember(givens{j(1),l,2},alt_B(:,j(1),:,l))); %length(intersect(alt_B(:,j(1),:,l),initial_B(:,j(1),:,l))) - 1;
score(3) = score_line(B(:,j(1),:,l)) + 100*alt_score; %9 - length(unique(B(:,j(1),:,l))) + 100*alt_score;
alt_score = sum(ismember(givens{j(2),l,2},alt_B(:,j(2),:,l))); %length(intersect(alt_B(:,j(2),:,l),initial_B(:,j(2),:,l))) - 1;
score(4) = score_line(B(:,j(2),:,l)) + 100*alt_score; %9 - length(unique(B(:,j(2),:,l))) + 100*alt_score;

end %score_spec_sudoku

function [counters,current_perm] = reset_perm(counters,current_perm,perm_len,i,j)
  counters{i,j} = [2,ones(1,perm_len(i,j)-1)];
  current_perm{i,j} = 1:perm_len(i,j);
end %reset_perm

function line_score = score_line(v)
line_score = sum(~diff(sort(v(:))));
end
