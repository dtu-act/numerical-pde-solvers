function C = calcNormConstant(P1,P2)
    Cxy = P2./P1;
    Cxy(isinf(Cxy)|isnan(Cxy)) = 0;
    
    N = length(find(abs(Cxy) > 0));
    C = sum(Cxy, 'all')/N;
end