clear all
close all

import utilsDD.*
import models.types.*
import models.*
import analysis.*

iter = 2000;

% "FDTD_FDTD" | "FOURIER_FOURIER" | "FOURIER_FDTD" | "SEM_SEM" | "SEM_FDTD";
compare_methods = "SEM_FDTD";

scheme1_order = 6;
scheme2_order = 2;
interface_order = 2;

fn = 1000;
ppw = 10;

xminmax = [0,2];

c=343;
rho=1.2;

[dx, nmodes] = calcSpatial(xminmax(2),fn,ppw,c);

iloc = InterfaceLocation1D.LEFT;

switch compare_methods
    case "FDTD_FDTD"
        % "spatialAdjustment" | "temporalAdjustment" | "none"
        adjust = "none";
        
        %compare_FDTD(iter, xminmax, dx, fn, c, rho, scheme1_order, interface_order)
        analysis.compare_FDTD_FDTD(adjust, iter, xminmax, dx, fn, c, rho, scheme1_order, scheme2_order, iloc)
    case "FOURIER_FOURIER"
        % "spatialAdjustment" | "temporalAdjustment" | "none"
        adjust = "none";
        
        analysis.compare_FOURIER_FOURIER(adjust, iter, xminmax, dx, nmodes, fn, c, rho, scheme1_order, scheme2_order, iloc)
    case "FOURIER_FDTD"
        % "spatialAdjustment" | "temporalAdjustment" | "none"
        adjust = "none";
        analysis.compare_FOURIER_FDTD(adjust, iter, xminmax, dx, nmodes, fn, c, rho, scheme1_order, scheme2_order, iloc)
    case "SEM_SEM"
        % TODO: assertion even number of elements required
        iface1 = SEMUniformInterface(dx/2, 4, 1, false, InterfaceLocation1D.LEFT);
        iface2 = SEMUniformInterface(dx/2, 4, 1, false, InterfaceLocation1D.RIGHT);
        adjust = "none";

        %compare_SEM(iter, xminmax, dx, fn, c, rho, scheme1_order, interface_order)
        analysis.compare_SEM_SEM(adjust, iter, xminmax, dx, fn, c, rho, ...\
            scheme1_order, scheme2_order, interface_order, interface_order, iface1, iface2, iloc)
    case "SEM_FDTD"
        uniform_distr = true;
        iface = SEMUniformInterface(dx, 2, scheme1_order, uniform_distr, iloc);
        
        analysis.compare_SEM_FDTD(iter, xminmax, dx, fn, c, rho, ...\
            scheme1_order, scheme2_order, interface_order, iface)
    otherwise
        error('method not supported')
end