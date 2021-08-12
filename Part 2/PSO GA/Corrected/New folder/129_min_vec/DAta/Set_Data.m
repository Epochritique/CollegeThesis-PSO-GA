for i=1:131
    for j=1:131
        if dist_ij(i,j) == -1
            dist_ij(i,j) = inf;
        end
    end
end

for k=1:131
    for i=1:131
        for j=1:131
            if dist_ij(i,j) > dist_ij(i,k) + dist_ij(k,j)
                dist_ij(i,j) = dist_ij(i,k) + dist_ij(k,j);
            end
        end
    end
end