function [lines1_mean, lines2_mean, Dkl] = bootstraplines(lines1,lines2, settings)

% This function computes mean traces after bootstrapping the data

%input:
%lines1 and lines2, traces based on which bootstrapping can be performed
%Settings =[x,y]. x determines how many bootstraps are performed, y sets
%how many lines should be included in each bootstrap sample. when y is not
%specified the bootstrap sample has the same size as the original data set. 

%output
%lines1_mean - average lines from the first setting obtained by each bootstrap run
%lines2_mean - average lines from the second setting obtained by each bootstrap run
%Dkl - kullback leibler divergence for the differce between the first and
%second set of line profiles
%%
% lines1 = Bdel_v_BRP_scaled;
% lines2 = CaM_v_BRP_scaled;
% settings = [10,15];
if length(settings)==1
r1 = randi ([1 size(lines1, 1)], settings(1), size(lines1, 1));
r2 = randi ([1 size(lines2, 1)], settings(1), size(lines2, 1));
else
    r1 = randi ([1 size(lines1, 1)], settings(1), settings(2));
r2 = randi ([1 size(lines2, 1)], settings(1), settings(2));
end
lines1_mean = zeros(settings(1), size(lines1,2));
lines2_mean = zeros(settings(1), size(lines2,2));
Dkl = zeros(settings(1),1);
for i =1:settings(1)
    lines1_set = lines1(r1(i,:),:);
    lines2_set = lines2(r2(i,:),:);
    lines1_mean(i,:) = (mean(lines1_set));
    lines2_mean(i,:) = (mean(lines2_set));
    if round(sum(lines1_mean(i,:)),2)~=1 | round(sum(lines2_mean(i,:)),2)~=1
        disp('distribution not summing to 1')
    end
    Dkl(i) = compute_Dkl(lines1_mean(i,:),lines2_mean(i,:));
end

figure()
plot(lines1_mean', 'k')
hold on 
plot(lines2_mean', 'r')
end

function Dkl = compute_Dkl(Q, P)
    for i = 1:length(P)
      Dkl_i(i) = P(i)*log(P(i)/Q(i));
    end
    Dkl = sum(Dkl_i);
end