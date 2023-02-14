clear all
close all

% NOTE: passivity not enforced - see vectorfitDriverLR.m inside
% libParanumal git code.

base_path = '/Users/nikolasborrel/data/libP/setups';
filename = sprintf('%s/freq_dep_lr.dat',base_path);

c_phys = 343.0; % m/s speed of sound
rho = 1.2;
c = c_phys;
c_norm = c/c_phys;
f_range = [50,5000]; % Range where boundary conditions are defined
f = (f_range(1):1:f_range(end))*c_norm; % Range where boundary conditions are defined

sigma = 8000*c_norm; % Flow resistivity
dmat = 0.10; % Thickness of porous material inside cavity
Npoles = 14; % Order of approximation on rational function for BC's

[impedance_data, Y, fit, fspan] = impedance_bcs.fitImpedanceBoundaryData(rho, c, f, sigma, Npoles, dmat);
impedance_data.type = "freq_dep";

impedance_bcs.plotImpedanceFitting(fspan, c_norm, Y, fit)

fileID = fopen(filename,'w');
fprintf(fileID,'%d %d %d\n',[Npoles,length(impedance_data.A),length(impedance_data.B)]);
fprintf(fileID,'%.12f\n',impedance_data.A);
fprintf(fileID,'%.12f\n',impedance_data.B);
fprintf(fileID,'%.12f\n',impedance_data.C);
fprintf(fileID,'%.12f\n',impedance_data.lambdas);
fprintf(fileID,'%.12f\n',impedance_data.alpha);
fprintf(fileID,'%.12f\n',impedance_data.beta);
fprintf(fileID,'%.12f\n',impedance_data.Yinf);
fprintf(fileID,'-----\n');

fprintf(fileID,'sigma = %f\n',sigma);
fprintf(fileID,'dmat = %f\n',dmat);
fprintf(fileID,'freqRange = [%d,%d]\n',f_range);

fclose(fileID);