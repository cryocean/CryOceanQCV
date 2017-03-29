function kcross = crossOver(x,y,orbit)
% USAGE : kcross = crossOver(x,y,orbit)
%
% Find the indices of pair of points corresponding to crossovers
%
% INPUTS :    x : vector containing the longitudes of all records
%             y : vector containing the latitudes of all records
%             orbit : vector the same size as x and y containing the orbit number for
%                     each record
%
% OUTPUTS :    kcross : 2 x n matrix containing the indices of the pair of points
%                       corresponding to crossovers
%
% author : Francisco Mir Calafat (francisco.calafat@noc.ac.uk)
% NOTE: this code assumes that each point is crossed once at most.

x = x(:);
y = y(:);
orbit = orbit(:);
orbitu = unique(orbit);
norbit = length(orbitu);

kcross = cell(norbit-1,1);
for j=1:norbit-1
    k1 = orbit == orbitu(j);
    k2 = find(~k1 & cumsum(k1)); % test 1 to 2:norbit, 2 to 3:norbit, etc.
    k1 = find(k1);
    x1 = x(k1);
    y1 = y(k1);
    x2 = x(k2);
    y2 = y(k2);
    
    if length(k2) > 1 % added on 4 Oct 2016
    [k,d] = knnsearch([x2 y2],[x1 y1]); % closest point (x2,y2) to (x1,y1)
    K = find(d' < 3*0.0585); % select close points only
    if length(K) < 2
        continue
    end
    Kd = diff(K);
    Kd = find(Kd ~= 1);
    nx = length(Kd)+1; % number of pair segments that may intersect
    Kd = [0 Kd length(K)]; %#ok
    % store the coordinates and indices of all pair of segments
    X1 = cell(nx,1);
    Y1 = cell(nx,1);
    X2 = cell(nx,1);
    Y2 = cell(nx,1);
    K1 = cell(nx,1);
    K2 = cell(nx,1);
    for i=1:nx
        K1{i} = unique(K(Kd(i)+1):K(Kd(i+1)));
        K2{i} = unique(k(K1{i}));
        % if there is only one point we add two more points to test for
        % possible intersection, but only if the prior and next points are
        % close enough to form a straight line
        if length(K1{i}) == 1
            if K1{i} ~= 1 && K1{i} ~= length(x1)
                ktmp = [K1{i}-1 K1{i} K1{i}+1];
            elseif K1{i} == 1
                ktmp = [K1{i} K1{i}+1 K1{i}+2];
            elseif K1{i} == length(x1)
                ktmp = [K1{i}-2 K1{i}-1 K1{i}];
            end
            d1 = sqrt((x1(ktmp(1))-x1(ktmp(2)))^2 + (y1(ktmp(1))-y1(ktmp(2)))^2);
            d2 = sqrt((x1(ktmp(3))-x1(ktmp(2)))^2 + (y1(ktmp(3))-y1(ktmp(2)))^2);
            if d1 < 10*0.0585 && d2 < 10*0.0585;
                K1{i} = ktmp';
            else
                continue
            end
        end
        if length(K2{i}) == 1
            if K2{i} ~= 1 && K2{i} ~= length(x2)
                ktmp = [K2{i}-1 K2{i} K2{i}+1];
            elseif K2{i} == 1
                ktmp = [K2{i} K2{i}+1 K2{i}+2];
            elseif K2{i} == length(x2)
                ktmp = [K2{i}-2 K2{i}-1 K2{i}];
            end
            d1 = sqrt((x2(ktmp(1))-x2(ktmp(2)))^2 + (y2(ktmp(1))-y2(ktmp(2)))^2);
            d2 = sqrt((x2(ktmp(3))-x2(ktmp(2)))^2 + (y2(ktmp(3))-y2(ktmp(2)))^2);
            if d1 < 10*0.0585 && d2 < 10*0.0585;
                K2{i} = ktmp';
            else
                continue
            end
        end
        X1{i} = x1(K1{i});
        Y1{i} = y1(K1{i});
        X2{i} = x2(K2{i});
        Y2{i} = y2(K2{i});
    end
    
    % test if segments intersect and calculate intersection point
    ind1 = zeros(nx,1);
    ind2 = zeros(nx,1);
    for i=1:nx
        if ~isempty(X1{i}) && ~isempty(X2{i})
            p0 = [X1{i}(1);Y1{i}(1)];
            p1 = [X1{i}(end);Y1{i}(end)];
            q0 = [X2{i}(1);Y2{i}(1)];
            q1 = [X2{i}(end);Y2{i}(end)];
            rhs = q0-p0;
            lhs = [p1-p0 q0-q1];
            xy = lhs\rhs;
            if sum(xy >= 0 & xy <=1) == 2
                xpoint = p0 + xy(1).*(p1-p0);
                kx1 = dsearchn([X1{i} Y1{i}],xpoint');
                ind1(i) = K1{i}(kx1);
                kx2 = dsearchn([X2{i} Y2{i}],[X1{i}(kx1) Y1{i}(kx1)]);
                ind2(i) = K2{i}(kx2);
            end
        end
    end
    ind1(ind1 == 0) = [];
    ind2(ind2 == 0) = [];
    kcross{j} = [k1(ind1) k2(ind2)];
    end % close "if length(k2) > 1". Added on 4 Oct 2016
end
kempty = cellfun(@isempty,kcross);
kcross(kempty) = [];
kcross = cell2mat(kcross)';
    
    
