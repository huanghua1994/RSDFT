 function [rho, lam, W] = rsdft(nev, Domain, Atoms, tol, maxits, fid) 
%% function [rho, lam, W] = rsdft(A, nev, Domain, Atoms, tol, maxits, fid) 
%% IN:
%% A      = sparse matrix representing the discretization of the
%%          Laplacean -- 
%% nev    = number of eigenvalues - this is the number of occupied states
%% Domain = struct containing information on the physical domain.
%% Atoms  = struct containing information on the atoms.
%% tol    = tolerance parameter for stopping scf iteration.
%% maxtis = maximum number of SCF iterations allowed. 
%% fid    = output file id
%%
%% OUT   : 
%% rho    = final charge density found
%% lam    = eigenvalues computed - their number may be larger than nev
%% W      = set of wave functions. 
%%============================================================ 
%%-------------------- include file --  

%%%
%%%  Read atoms types and positions
%%%  Derive everything else!
%%%
%%%%Defaut Parameters
%%%%
fd_order = 8;                         %% order of finite difference scheme 
maxits   = 40;                         %% max SCF iterations 
tol      = 1.e-05;                     %% tolerance for SCF iteration.
Fermi_temp  =  1000.0 ;                %%  Smear out fermi level

%%%%%%%%%
%%%%%%%%%  Grid spacing from table (takes smallest h)
%%%%%%%%%

%% Boundary_Sphere_Radius 
 sph_rad = 6. ; 
%% grid-spacing 
 grid_spacing  = 0.5 ; 
%% number of states  
 num_states  =  16 ; 






fid = fopen('./rsdft.out', 'w');       %% create a output file

h  = grid_spacing;                     %% define phys. domain and grid spacing

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
  fprintf(fid, ' and their coordinates are:\n\n');
  fprintf(fid, '\tx [a.u.]\t\ty [a.u.]\t\tz [a.u.]\n');
  atom_count = atom_count + length(xyz(:,1));
  for i=1:length(xyz(:,1)) 
    fprintf(fid, '\t%.6f\t\t%.6f\t\t%.6f\n', xyz(i,1), xyz(i,2), xyz(i,3));
  end
  fprintf(fid, '\n');
end
fprintf(fid,' --------------------------------------------------\n');
fprintf(fid,' Total number of atoms :         %d\n\n', atom_count);
fprintf(fid,' Number of states:               %10d \n', num_states);
fprintf(fid,' h grid spacing :                %10.5f  \n', h);
fprintf(fid,' Hamiltonian size :              %10d  \n', nx*ny*nz);
fprintf(fid,' Sphere Radius :                 %10.5f   \n', sph_rad);
fprintf(fid,' # grid points in each direction %10d  \n', nx);
fprintf(fid,' --------------------------------------------------\n');
%% construct Laplacian operator
A  = (1/(h*h))*fd3d(nx, ny, nz, fd_order);










include;
%parsec_dat;         %% duplicated read... but OK for now
n = size(A,1);
Hpot = zeros(n,1);
pot = Hpot;
index0=length(pot);
poldeg;
err = 1 + tol;
its = 0;
Ry=13.605698066;
%%
%%  Set up ion-core potentials/initial screening
%%  Import screening from the atom, i.e., hartree potential and charge.
%%
%%  Ion-Ion repulsion
[E_nuc0]=nuclear(Domain, Atoms);
disp('Nuclear repulsion in eV')
disp(E_nuc0*Ry)
%%
%%--------------------  Local part of the pseudopotential
%%
[rho0, hpot0, Ppot]  = pseudoDiag(Domain, Atoms);
%%%%%%%
hpsum0=sum(rho0.*hpot0)*Ry;

fprintf(fid,' Initial Hartree energy (eV) = %10.5f  \n', hpsum0) ;
%%-------------------- count # atoms for stats
n_atoms  = 0;
for ii = 1:length(Atoms) 
  n_atoms = n_atoms + length(Atoms(ii).coord(:,1));
end
%%
%%--------------------  Non-local part of the pseudopotential
%%
      [vnl] = pseudoNL(Domain, Atoms) ; 
%%
%%--------------------  Find a initial guess for the screening potential
%%--------------------  Use screening from Gaussian density
%%
  indx1=length(rho0);
  h=Domain.h;
%%
  rhoxc = rho0' ./ h^3;
  [XCpot,exc] = exc_nspn(Domain, rhoxc, fid);
  xcpot=XCpot';
  Nelec = nelectrons(Atoms);
  Fermi_temp;
%%
%%------------------- open output file (wfn.dat)
%%
  wfnid = fopen('./wfn.dat', 'wb');
  indx2=length(Ppot);
  indx3=length(hpot0);
  indx4=length(XCpot);
%% indx5=length(vnl)
   pot = Ppot + hpot0 + 0.5*xcpot  ;
%%-------------------- SCF LOOP
%%-------------------- when 'preconditioning' is used fall ilu0
 if (CG_prec) 
      disp('calling ilu0 ...') 
      [L, U] = luinc(A,'0');
      disp(' done.') 
      PRE = struct('L',L, 'U',U);
 end

%%-------------------- clear persistent variables in mixer.
 clear mixer;
%%-------------------- SCF LOOP starts here
 fprintf(fid, '\n----------------------------------\n\n');
 while (err > tol & its <= maxits) 
     its = its+1;
     fprintf(1,'  SCF iter # %d  ... ',its)
%%-------------------- redefine Hamiltonian 
     B = 0.5* A + spdiags(pot, 0, n, n) + vnl ;
%%--------------------diagmeth defined in include
     tic;
     if (diagmeth ==1 | (its == 1 & diagmeth == 0)) 
        disp('calling lanczos..') 
        v = randn(n,1); 
        [W, lam] = lan(B, nev+15, v, nev+400, 1.e-06);
     elseif (its == 1 & diagmeth == 2)
        disp('calling chsubsp..') 
        v = randn(n,1); 
        [W, lam] = chsubsp(poldeg, nev+15, B) ;
     else 
        disp('calling chebsf..') 
        [W, lam] = chefsi1(W, lam, poldeg, nev, B) ;
     end
%%
     diag_time = toc;
%%---------------------print results
     fprintf(fid,' \n \n SCF iter # %d  ... \n',its);
     fprintf(fid,' Diagonalization time [sec] :\t%f\n\n', diag_time);
%%-------------------- get occupation factors and fermi-level 
 [Fermi_level, occup] = occupations(lam(1:nev), Fermi_temp, Nelec, 0.000001);
%%
 fprintf(fid, '   State  Eigenvalue [Ry]     Eigenvalue [eV]\n\n');
     for i = 1:nev
       eig = lam(i) * 2*Ry;
       ry = eig / Ry; 
       occ=occup(i);
       fprintf(fid, '%5d   %15.10f   %18.10f  %5.2f\n', i, ry, eig, occ); 
     end
%%-------------------- get charge density
     rho = (W(:,1:nev) .* W(:,1:nev)) *2* occup ;
     hrhs = (4*pi/h^3)*(rho-rho0);
     indx=length(rho)  ;
     for j=1:indx    
         rho(j)=rho(j)/h^3;
     end
%% trigger timer
     tic;
     if (CG_prec) 
        Hpot = pcg (A, hrhs, Hpot, 200, 1.e-07, PRE,'precLU');
     else 
        Hpot = pcg (A, hrhs, Hpot, 400, 1.e-07);
     end
     hart_time = toc;
     fprintf(fid, '\nHartree potential time [sec]: \t%f\n\n', hart_time);
     
     [XCpot,exc] = exc_nspn(Domain, rho, fid);
     HHpot=Hpot; 
     potNew = Ppot+0.5*XCpot+Hpot+hpot0;
     err = norm(potNew - pot) / norm(potNew);
     fprintf(fid,'   ...  err = %10.2e  \n', err) ;
     fprintf(1,'   ...  err = %10.2e  \n', err) ;
%-------------------- call mixer 
     pot = mixer(pot, potNew-pot);
%%   pot= 0.4*potNew+0.6*pot; % Simple mixing.
%%-------------------- total energy 
%%  Sum over eigenvalues--Put everything in Ryd
    Esum= sum(lam(1:nev).*occup(1:nev)); 
    Esum0=4*Esum;
%%
%%--------------------   Hartree potential
%%  Factor of two for double counting--converts to Ryd
    Hsum0=sum(rho.*(Hpot+hpot0))*h^3;
%%
%%-------------------- Exchange correlaion sum
%% No factor of two because energy is in Ry
%%
Vxcsum0=sum(rho.*XCpot)*h^3;
Excsum0=exc;
%%--------------------  Total electronic energy
%%
    E_elec0=  Esum0-Hsum0+Excsum0-Vxcsum0;
%%   Add in nuclear-nuclear repulsion term 
    E_total0 = E_elec0 + E_nuc0;
%%
%%--------------------  Convert to eV
    Esum=Ry*Esum0;
    Hsum=Ry*Hsum0;
    Excsum=Ry*Excsum0;
    E_nuc=Ry*E_nuc0;
    E_total=Ry*E_total0;
    E_total2=E_total0;
%%-------------------- actual printing of info on energies..
  fprintf(fid,'  \n\n');
  fprintf(fid,' Total Energies \n\n');
  fprintf(fid,' Sum of eigenvalues      = %10.5f  eV   = ',Esum) ;
  fprintf(fid,'  %10.5f  Ry  \n',Esum/Ry) ;
  fprintf(fid,' Hartree energy          = %10.5f  eV   = ',Hsum) ;
  fprintf(fid,'  %10.5f  Ry  \n',Hsum/Ry) ;
  fprintf(fid,' Exchange-corr. energy   = %10.5f  eV   = ',Excsum) ;
  fprintf(fid,'  %10.5f  Ry  \n',Excsum/Ry) ;
  fprintf(fid,' Ion-ion repulsion       = %10.5f  eV   = ',E_nuc) ;
  fprintf(fid,'  %10.5f  Ry  \n',E_nuc/Ry) ;
  fprintf(fid,' Tot. electronic energy  = %10.5f  eV   = ',E_total) ;
  fprintf(fid,'  %10.5f  Ry  \n',E_total0) ;
  fprintf(fid,' Electronic energy/atom  = %10.5f  eV   = ',E_total/n_atoms) ;
  fprintf(fid,'  %10.5f  Ry  \n',E_total0/n_atoms) ;
%%
end
%%-------------------- free memory (persistent variables) in mixer.
clear mixer;
%%------------------- Output results
    fprintf(fid, '\n Finished\n');
   
    little_big_test = 26;
    fwrite(wfnid, little_big_test, 'uint32');
     
    fwrite(wfnid, Domain.radius, 'uint32');
    fwrite(wfnid, Domain.h, 'double');
     
    pot_length = length(pot(:,1));
    fwrite(wfnid, pot_length, 'uint32');
    fwrite(wfnid, pot, 'double');
     
     rho_length = length(rho(:,1));
     fwrite(wfnid, rho_length, 'uint32');
     fwrite(wfnid, rho, 'double');
     
     w_length = length(W(:,1));
     w_col_length = length(W(1,:));
     fwrite(wfnid, w_length, 'uint32' );
     fwrite(wfnid, nev, 'uint32');
     for i = 1:nev
       fwrite(wfnid, W(:, i), 'double');
     end

fclose(wfnid);
