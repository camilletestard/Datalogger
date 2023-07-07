function [draws] = randomDraw(distribution,n)
if (sum(distribution) ~= 1)
    error('Input distribution is not valid. sum(dist) ~= 1')
else
    cumDist = cumsum(distribution);
    for i=1:n
        r = rand;
        inds=find(cumDist > r);
        draws(i)=inds(1);
    end
end
end

