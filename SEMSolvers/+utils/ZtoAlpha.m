function alpha = ZtoAlpha(Z,c,rho)
    alpha = 1 - ((Z - rho*c)/(Z + rho*c))^2;
end