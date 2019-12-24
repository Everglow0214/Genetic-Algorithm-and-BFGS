clear all

pop_size = 400; % Rows
chrom_len = 40; % Columns
chrom_num = 2; % Variable(s) in the objective function
chrom_range = [-1 -1;1 1]; % Range of variable(s)
iteration = 4000;
pc = 0.6; % Probability of crossover
pm = 0.4; % Probability of mutation

pop = initpop(pop_size, chrom_len); % Initialization
obj_value = cal_obj_value(pop, chrom_num, chrom_range);
[best_individual_binary, best_value] = best(pop, obj_value);
best_individual = binary2decimal(best_individual_binary, chrom_num, chrom_range);

result = zeros(1,chrom_num+1); % The first several elements are the best variable, and the
%                                last one is the optimal value
result(1:chrom_num) = best_individual;
result(chrom_num+1) = best_value;

for i = 1:iteration
    obj_value = cal_obj_value(pop, chrom_num, chrom_range);
    fit_value = obj_value;
    new_pop = select(pop, fit_value);
    new_pop = crossover(new_pop, pc);
    new_pop = mutation(new_pop, pm);
    [best_individual_binary, best_fit] = best(pop, fit_value);
    
    % Determine which one is optimal
    if result(chrom_num+1) < best_fit
        best_individual = binary2decimal(best_individual_binary, chrom_num, chrom_range);
        result(1:chrom_num) = best_individual;
        result(chrom_num+1) = best_fit;
    end
end
result(chrom_num+1) = -result(chrom_num+1)

% Initialization
function pop = initpop(pop_size, chrom_len)
    pop = round(rand(pop_size, chrom_len));
end

% Convert the binary to a decimal number
function pop2 = binary2decimal(pop, chrom_num, chrom_range)
    [len_x,len_y] = size(pop);
    gene_binary = len_y/chrom_num;
    for i = 1:chrom_num
        for j = 1:(len_y/chrom_num)
            % 1(or 0) * 2^(k-1)
            pop1(:,(i-1)*gene_binary+j) = 2.^(gene_binary-j).*pop(:,(i-1)*gene_binary+j);
        end
    end
    temp = zeros(len_x,chrom_num);
    pop2 = zeros(len_x,chrom_num);
    for i = 1:chrom_num
        % Summation
        temp(:,i) = sum(pop1(:,(1+gene_binary*(i-1)):(gene_binary*i)),2);
        % Make the variables stay in the range
        % For example, [0,1]*4+2 -> [2,6]
        pop2(:,i) = temp(:,i)*(chrom_range(2,i)-chrom_range(1,i))/(2^gene_binary-1)+chrom_range(1,i);
    end
end        

% Calculate values of the objective function, given binary arguments
function [obj_value] = cal_obj_value(pop, chrom_num, chrom_range)
    x = binary2decimal(pop, chrom_num, chrom_range);
    [len_x,len_y] = size(pop);
    obj_value = zeros(len_x,1);
    if chrom_num == 1
        obj_value = -(-10.*cos(3.*x).^2-(x-5).^2+250);
    elseif chrom_num == 2
        x1 = x(:,1);
        x2 = x(:,2);
        % obj_value = -(3.*(1-x1).^2.*exp(-x1.^2-(x2+1).^2) - 10.*(x1./5-x1.^3-x2.^5).*exp(-x1.^2-x2.^2)  - 1/3.*exp(-(x1+1).^2-x2.^2));
        % obj_value = -(x1.^2 + x2.^2 - 0.5.*cos(pi.*x1) - 0.5*cos(2.*pi.*x2) + 1);
        obj_value = -(x1.^2 + x2.^2 - 0.7.*cos(2.*pi.*x1).*cos(3.*pi.*x2) + 0.7);
        % obj_value = -(x1.^2 + 2.*x2.^2 - 0.3.*cos(4.*pi.*x1) - 0.3.*cos(5.*pi.*x2) + 0.6);
    else
        disp('Wrong variable numbers!!');
    end
end  

function [new_pop] = select(pop, fit_value)
    [len_x,len_y] = size(pop);
    total_fit = sum(fit_value);
    p_fit_value = fit_value/total_fit;
    p_fit_value = cumsum(p_fit_value);
    ms = sort(rand(len_x,1));
    fitin = 1;
    newin = 1;
    while newin <= len_x
        if (ms(newin))<p_fit_value(fitin)
            new_pop(newin,:) = pop(fitin,:);
            newin = newin + 1;
        else
            fitin = fitin + 1;
        end
    end
end    

% Crossover
function [new_pop] = crossover(pop, pc)
    [len_x,len_y] = size(pop);
    new_pop = ones(len_x,len_y);
    for i = 1:2:len_x-1
        if rand < pc
            cpoint = round(rand*len_y);
            new_pop(i,:) = [pop(i,1:cpoint), pop(i+1,cpoint+1:len_y)];
            new_pop(i+1,:) = [pop(i+1,1:cpoint), pop(i,cpoint+1:len_y)];
        else
            new_pop(i,:) = pop(i,:);
            new_pop(i+1,:) = pop(i+1,:);
        end
    end
end

% Mutation
function [new_pop] = mutation(pop, pm)
    [len_x,len_y] = size(pop);
    new_pop = ones(len_x,len_y);
    for i = 1:len_x
        if rand < pm
            mpoint = round(rand*len_y);
            if mpoint <= 0
                mpoint = 1;
            end
            new_pop(i,:) = pop(i,:);
            if new_pop(i,mpoint) == 0
                new_pop(i,mpoint) = 1;
            elseif new_pop(i,mpoint) == 1
                new_pop(i,mpoint) = 0;
            else
                disp('Wrong!!')
            end
        else
            new_pop(i,:) = pop(i,:);
        end
    end
end

% Obtain the optimal solution (with the greatest value of the objective function)
function [best_individual,best_fit] = best(pop,fit_value)
	[len_x,len_y] = size(pop);
    best_individual = pop(1,:);
    best_fit = fit_value(1);
    for i = 2:len_x
        if fit_value(i) > best_fit
            best_individual = pop(i,:);
            best_fit = fit_value(i);
        end
    end
end