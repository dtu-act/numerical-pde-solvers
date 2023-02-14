function plotImpedanceFitting(f, c_norm, Y, fit)
    f_range = [min(f), max(f)]*c_norm;
    figure()
    subplot(221)
    semilogx(f*c_norm,imag(Y))
    hold on
    semilogx(f*c_norm,imag(fit))
    xlabel('Frequency [Hz]');
    ylabel('Im(Y)')
    legend('Initial admittance curve','Rational function fit');
    xlim(f_range)

    subplot(222)
    semilogx(f*c_norm,real(Y))
    hold on
    semilogx(f*c_norm,real(fit))
    xlabel('Frequency [Hz]');
    ylabel('Re(Y)')
    legend('Initial admittance curve','Rational function fit');
    xlim(f_range)

    subplot(223)
    semilogx(f*c_norm,angle(Y)*180/pi)
    hold on
    semilogx(f*c_norm,angle(fit)*180/pi)
    xlabel('Frequency [Hz]');
    ylabel('\angle Y')
    legend('Initial admittance curve','Rational function fit');
    xlim(f_range)


    subplot(224)
    semilogx(f*c_norm,abs(Y))
    hold on
    semilogx(f*c_norm,abs(fit))
    xlabel('Frequency [Hz]');
    ylabel('|Y|')
    legend('Initial admittance curve','Rational function fit');
    xlim(f_range)
    
    hold off
end