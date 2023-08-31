clear all
close all
addpath('../shared')
addpath('../SEMSolvers')

[path_ref_dir] = paths.setupPaths("LOCAL");

compare_with_SEM = true;

% === PARAMETERS === %
Nw_greens = 500;

f = 300;
c = 343; rho = 1.225;
acoustic_params = models.AcousticParameters(f,c,rho);
boundary_type = models.BoundaryCondition.Velocity;
lx = 2;
ly = 2;
bbox = [0,0;lx,ly];

Porder = 4;
%ppw = 20;
%[dx, NeX, NeY] = utils.calcGridResolution(f,acoustic_params.c,Porder,ppw,lx,ly);
dx = 0.1;
NeX = lx/dx;
NeY = ly/dx;
assert(NeX == round(NeX) && NeY == round(NeY));

xy0_rel = [3/5,3/5];
xy0 = meshing.calculateSourcePosition(NeX,NeY,bbox,xy0_rel);

source = models.PointSource(xy0,1);
% sigma = 0.2; % -> 1000Hz % source = setup.gauss2DSetup(c,sigma,xy0);

uniform_mesh = true;

r0 = 0.4;

% === SOLVE === %
XY = meshing.mesh2D(bbox, NeX, NeY);
P = greens.greensNeumannUnstructuredSpatial2D(Nw_greens,Nw_greens,XY,xy0,acoustic_params.k);

if compare_with_SEM
    path_data = paths.refDataPath(path_ref_dir,'SEM',Porder,NeX,NeY,acoustic_params,xy0,bbox,boundary_type);
    data = load(path_data);
    Pref = data.P;
    conn = data.conn;
    Pref = Pref(unique(conn(:,1:3)));
    fd = @(p) min(-drectangle(p,0,lx,0,ly), dcircle(p,xy0(1),xy0(2),r0));
    Proi = validation.extractROI(P(:),XY,fd);
    Proi_ref = validation.extractROI(Pref(:),XY,fd);
    %Proi = validation.calcNormConstantMax(Proi,Proi_ref)*Proi;

    l2_err = norm(Proi-Proi_ref, 2);
    fprintf('L2 err: %e\n', l2_err)

    plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'SEM', 'GREENS', 90)
    plotting.plotHelmholtzError(Proi, Proi_ref, XY(:,1), XY(:,2), 'SEM', 'GREENS', 45)
else
    fd = @(p) min(-drectangle(p,0,lx,0,ly), dcircle(p,xy0(1),xy0(2),r0));
    Proi = validation.extractROI(P(:),XY,fd);
    plotting.plotHelmholtz(Proi, XY(:,1), XY(:,2), 'SEM')
end