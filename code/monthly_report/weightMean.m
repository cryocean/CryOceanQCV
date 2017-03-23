function [mu,s] = weightMean(z,y,lat,dy)

% weighted mean
wi = cosd(lat);
Ki = true(size(lat));
mu = 0;
for j=1:length(lat)
    k = y >= lat(j)-dy & y < lat(j)+dy;
    if sum(k) == 0
        Ki(j) = false;
    else
        mu = mu + mean(z(k))*wi(j);
    end
end
mu = mu/sum(wi(Ki));

%weighted standard deviation
ai = cosd(lat);
Ki = true(size(lat));
s = 0;
for j=1:length(lat)
    k = y >= lat(j)-dy & y < lat(j)+dy;
    if sum(k) == 0
        Ki(j) = false;
    else
        nk = sum(k);
        wi = ai(j)/nk;
        s = s + sum(wi.*(z(k)-mu).^2);
    end
end
s = s/sum(ai(Ki));
s = sqrt(s);

