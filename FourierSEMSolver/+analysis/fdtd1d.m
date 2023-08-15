%parameters
c=343;
spatialWidth=100;
len=2; %waveguide length (in m)

% num of nodal points in the space direction
dx=len/spatialWidth;

%temporal resolution. Courant s=c*dt/dx = 1 <-> dt = dx/c?
s=1.0;
dt=s*dx/c;

%medium density and medium bulk modulus
rho=1.21;
K=rho*c^2;

%num iterations
iter=200;

%excitation point
i0=spatialWidth/2;

%init
p=zeros(spatialWidth,1);
v=zeros(spatialWidth,1);

sigma=0.0005;
t0=3*sigma;
s = @(t) exp(-((t-t0)./0.0005).^2 ); %gaussian excitation pulse

for n=2:iter
  t=n*dt;
  for i=2:spatialWidth-1    
      p(i) = p(i)-K*dt/dx*(v(i+1)-v(i));
      if i==i0 
          p(i) = p(i) + s(t);
      end
      v(i)=v(i)-(1/rho)*dt/dx*(p(i)-p(i-1));
  end
  figure(1)
  plot(p)
  axis([0 spatialWidth -1.5 1.5]);
  drawnow
end