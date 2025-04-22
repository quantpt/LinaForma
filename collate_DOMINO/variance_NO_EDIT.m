% This script calculates the number of phases
function [v,ass] = variance_NO_EDIT(data)

data = data(:,2:end);
for i = 1:size(data,1)
    t = data(i,:); t = table2cell(t);

% Find number of rows
empty = cellfun(@isempty, t);
phases = length(t) - sum(empty);

% Find all instances of (2)
tmp2 = cellfun(@(x) contains(x, '(2)'), t);
tmp = sum(tmp2);
add = tmp;
phases = phases + add;

% Find all instances of (3)
tmp3 = cellfun(@(x) contains(x, '(3)'), t);
tmp = sum(tmp3);
add = tmp*2;
phases = phases + add;

% Find all instances of (4)
tmp4 = cellfun(@(x) contains(x, '(4)'), t);
tmp = sum(tmp4);
add = tmp*3;
v(i) = phases + add;
ass{i} = strjoin({t{:}});
end
v= v';
ass = ass';
end