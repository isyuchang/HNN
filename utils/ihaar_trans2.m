function out = iHaart2(LL,HL,LH,HH,M,N)
[m,n,c] = size(LL); m = 2*m; n = 2*n;
out = zeros(m,n,c);
out(1:2:m, 1:2:n, :) = LL - HL - LH + HH;
out(2:2:m, 1:2:n, :) = LL - HL + LH - HH;
out(1:2:m, 2:2:n, :) = LL + HL - LH - HH;
out(2:2:m, 2:2:n, :) = LL + HL + LH + HH;
out = out/2;
out = out(1:M, 1:N, :);
end

