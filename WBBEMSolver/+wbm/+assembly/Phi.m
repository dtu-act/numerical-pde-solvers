function out = Phi(phi,w_vec,nw)
    W = repmat(w_vec, [1,nw]);
    out = (W.*phi);
end