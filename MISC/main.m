clear all

profile on

%% add mixer to the path
path(path,'MIXER');
path(path,'PSEUDO');
path(path,'SPLINES');

[s,w] = unix('python < parse_inp.py'); %% generate parsec.dat
parsec_dat;                            %% execute to load data.

fd_order = 4;                          %% order of finite difference scheme 
maxits = 40;                           %% max SCF iterations 
tol    = 1.e-04;                       %% tolerance for SCF iteration.

fid = fopen('./rsdft.out', 'w');       %% create a output file

h  = grid_spacing;                     %% define physical domain and grid spacing

nx = fix(2*sph_rad/h) + 1; 
h  = 2*sph_rad/(nx-1);                 %% re-adjust h 
ny = nx;                               %% [-radis, radius] in each direction
nz = nx;
Domain = struct('radius', sph_rad,'nx',nx,'ny',ny,'nz',nz,'h', h);

fprintf(fid, ' Number of states: \t%d\n\n', num_states);
fprintf(fid,  'Atom data:\n -------------\n');
fprintf(fid, ' Total # of atom types is %d\n', at);
atom_count = 0;
for atm_typ = 1:at
  xyz = Atoms(atm_typ).coord;
  fprintf(fid, ' There are %d %s atoms\n', length(xyz(:,1)), Atoms(atm_typ).typ);
  fprintf(fid, ' and their initial corrdinates are:\n\n');
  fprintf(fid, '\tx [a.u.]\t\ty [a.u.]\t\tz [a.u.]\n');
  atom_count = atom_count + length(xyz(:,1));
  for i=1:length(xyz(:,1)) 
    fprintf(fid, '\t%.6f\t\t%.6f\t\t%.6f\n', xyz(i,1), xyz(i,2), xyz(i,3));
  end
  fprintf(fid, '\n');
end
fprintf(fid, ' Total number of atoms = %d\n\n', atom_count);

fprintf(fid, ' Number of states: \t%d\n\n', num_states);
fprintf(fid,' h grid spacing %10.5f  \n', h);
fprintf(fid,' Radius %10.5f  \n', sph_rad);
fprintf(fid,' Number of grid points %10.5f  \n', nx);

%% construct Laplacian operator
A  = (1/(h*h))*fd3d(nx, ny, nz, fd_order);

%% call real space dft code
[rh vals Wave] = rsdft(A, num_states, Domain, Atoms, tol, maxits, fid);

fclose(fid); %% Close output file

%% Display Graphics
rad = sph_rad;
figure(1);
title('Atomic Positions');
box = [-rad rad -rad rad -rad rad];
axis(box) 
hold
for atm_typ = 1:at
  xyz = Atoms(atm_typ).coord;
  for i=1:length(xyz(:,1)) 
    plot3(xyz(i,1),xyz(i,2),xyz(i,3),'o','linewidth',5)
  end
end



nplane=nx/2+0.5
rhoplot = reshape(rh,nx,ny,nz);

figure(2) 
title('Charge density at z = middle vertical plane');
contour(rhoplot(:,:,nplane));


figure(3) 
title('Charge density at z = middle vertical plane');
surface(rhoplot(:,:,nplane));

profile viewer
p = profile('info');
profsave(p,'profile_results')



