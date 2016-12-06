function G = calculateHedgesG(input)

%input is an n x 3 matrix containing:
% column 1 = measurements
% column 2 = population identifiers
% column 3 = trait identifiers
% all of the above coded as numbers

npops = length(unique(input(:,2)));
ntraits = length(unique(input(:,3)));

meaningfulComps = npops * (npops - 1) / 2;
rows = meaningfulComps * ntraits;

G = zeros(rows,1);

wr = 1;
for t = 1:ntraits
    wt = find(input(:,3)==t);
    for i = 1:(npops-1)
        wi = find(input(:,2)==i);
        for j = (i+1):npops
            wj = find(input(:,2)==j);
            is = intersect(wt,wi);
            js = intersect(wt,wj);
            p1 = input(is,1);
            p2 = input(js,1);
            m1 = mean(p1);
            m2 = mean(p2);
            n1 = length(p1);
            n2 = length(p2);
            v1 = var(p1);
            v2 = var(p2);
            sstar = sqrt( (((n1-1)*v1) + ((n2-1)*v2))/(n1+n2-2) );
            G(wr) = (m1-m2)/sstar;
            wr = wr + 1;
        end
    end
end