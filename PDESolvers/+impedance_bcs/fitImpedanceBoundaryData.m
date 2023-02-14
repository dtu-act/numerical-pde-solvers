function [impedance_data, Y, fit, f_range] = fitImpedanceBoundaryData(rho, c, f_range, sigma, N, d)        
    Zc = rho*c*(1+0.07*(f_range./sigma).^(-0.632) - 1i*0.107*(f_range./sigma).^(-0.632));
    kt = 2*pi*f_range/c.*(1+0.109*(f_range./sigma).^(-0.618) - 1i*0.16*(f_range./sigma).^(-0.618));
    
    Z = -1i*Zc.*cot(kt*d);
    Y = 1./Z;
    
    % Fit impedance to a rational function
    w = 2*pi*f_range;
    s = 1i*w;
    opts.cmplx_ss = 1; % imaginary output: 0 false, 1 true
    opts.spy2 = 0;
    Ns = length(s);
    weight = ones(1,Ns);
    
    % Complex starting poles
    bet=linspace(w(1),w(Ns),N/2);
    poles=[];
    for n=1:length(bet)
      alf=-bet(n)*1e-2;
      poles=[poles (alf-1i*bet(n)) (alf+1i*bet(n)) ]; 
    end
    
    Niter=3;
    for iter=1:Niter
      if iter==Niter, opts.skip_res=0; end
      disp(['   Iter ' num2str(iter)])
      [SER,poles,rmserr,fit]=vectfit3(Y,s,poles,weight,opts);
      rms(iter,1)=rmserr;
    end

    % We know it returns 2 real poles and one complex pole
    lambda(1) = -full(SER.A(1,1));
    lambda(2) = -full(SER.A(2,2));
    assert(isreal(lambda(1))); assert(isreal(lambda(2)))
    
    al = -full(real(SER.A(3,3))); 
    be = -full(imag(SER.A(3,3)));

    A(1) = SER.C(1);
    A(2) = SER.C(2); 
    assert(isreal(A(1))); assert(isreal(A(2)))
    B = real(SER.C(3));
    C = imag(SER.C(3));
    Yinf = SER.D;
    
    impedance_data.lambdas = lambda;
    impedance_data.alpha = al;
    impedance_data.beta = be;
    impedance_data.Yinf = Yinf;
    impedance_data.A = A;
    impedance_data.B = B;
    impedance_data.C = C;