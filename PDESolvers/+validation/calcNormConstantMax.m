function C = calcNormConstantMax(P1,P2)
    %[~,i] = max(abs(P1));
    [~,i] = max(P1);
    C = P2(i)/P1(i);    
end