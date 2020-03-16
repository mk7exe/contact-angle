function w = weighted_distribution(Distribution,Weight)

D=find(Distribution>-1);

% calculating weighted histogram acording to the area of triangles
edges = [0:180];
vals = Distribution(D);
weights = Weight(D);

Nedge = length(edges);
w = zeros(size(edges));
    
for n = 1:Nedge-1
    ind = find(vals >= edges(n) & vals < edges(n+1));
    if ~isempty(ind)
        w(n) = sum(weights(ind));
    end
end

ind = find(vals == edges(end));
if ~isempty(ind)
    w(Nedge) = sum(weights(ind));
end
end
