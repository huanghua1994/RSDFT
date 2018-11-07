clear all

%%% Variables and definitions
%% Laplacian = sparse matrix representing the discretization of the Laplacian
%% n_occ  = number of eigenvalues - this is the number of occupied states
%% Domain = struct containing information on the physical domain.
%% Atoms  = struct containing information on the atoms.
%% tol    = tolerance parameter for stopping scf iteration.
%% maxtis = maximum number of SCF iterations allowed. 
%% fid    = output file id
%%
%% rho = final charge density found
%% ev  = eigenvalues computed - their number may be larger than n_occ
%% psi = set of wave functions. 
%%============================================================ 
%%-------------------- include file --  


%% add mixer to the path
path(path,'MIXER');
%path(path,'PSEUDO');
path(path,'SPLINES');
%%%
%%%  Read atoms types and positions
%%%  
%%%  Derive everything else!
%%%
%%%   Units
include;

% Use the same random seed to get the same result for the same input
rng(2018);  

[fid,message] = fopen('parsec.dat','r');
N_atoms=0;
if fid == -1
    disp(message)
end
disp(' ********************** ')
disp(' DATA INPUT FOR PARSEC')
disp(' *********************  ')
disp('------------------------')
%in_data=input('Input data mode: 1 for manual, 2 for file:  ' );
%disp('   ')
in_data = 2;
if in_data==1
    at=input('Input number of different atomic species: '  );
    disp('   ')
    for i=1:at
        typ = input('Element of species, e.g., Mg, C, O, only first 18 elements supported: ', 's' );
        disp('  Coordinates should be in atomic units ')
        disp('  Example atoms at (0,0,0) should be entered as 0 0 0 on each line ')
        disp('  Terminate with /, i.e., 0 0 0 / for the last entry ')
        xyz=[];
        readxyz=[' '];
        while isempty(strfind(readxyz,'/'))
            readxyz=[readxyz,' ',input('  Input coordinates ','s')];
        end
        xyz=reshape(str2num(readxyz(1:strfind(readxyz,'/')-1)),3,[])';
        n_atom(i)=size(xyz,1);
        % xyz =reshape(xyz,3,[])';
        str=struct('typ',typ,'coord',xyz);
        Atoms(i) = str;
    end
end

if in_data==2
    at = fscanf(fid,'%g');
    for i=1:at
          typ = fscanf(fid,'%s',1);
          readxyz = fscanf(fid,'%g');
          n_atom(i)=length(readxyz)/3;
          xyz = reshape(readxyz,3,[])';
          str=struct('typ',typ,'coord',xyz);
          Atoms(i) = str;
    end
    status = fclose(fid);
end
%%%
%%%%Default Technical Parameters
%%%%
fd_order = 8;                         %% order of finite difference scheme 
maxits   = 40;                         %% max SCF iterations 
tol      = 1.e-04;                     %% tolerance for SCF iteration.
Fermi_temp  =  1000.0 ;                %%  Smear out fermi level
%%%%%%%%%
%%%%%%%%%  Grid spacing from table (takes smallest h)
%%%%%%%%%
% Information about each element
elem    = importdata('elements_new.csv');
N_elements=size(elem.data);
%%%%
%%%%
% This is the number of species that need to be looked at
N_types = length(Atoms);

zelec=0.;
for at_typ=1:N_types
    % This gets the first atom and assigns the atomic symbol to typ
    typ = Atoms(at_typ).typ;
    % This iterates through elem looking for the matching element data

    %%% Loop over number of atoms
    hmin=100.;
    for i = 1: N_elements
        if strcmp(typ,elem.textdata{i})
            index=i;
            Z=elem.data(index,2)*n_atom(at_typ);
            h=elem.data(index,4);
            if h<hmin
                hmin=h;
            end
            zelec=zelec+Z;
        end
    end
end
%%%%
%%%%  hmin is the smallest recommend grid for a particular atom
%%%%  n_occ is the number of eigenvalues to be fond
%%%%
h=hmin;
num_states=zelec+4;
if num_states < 6
    num_states=6;
end
n_occ=num_states;


%%%%
%%% Estimate the radius
%%%
% Find ion-ion core repulsion
rmax=0.;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_types = length(Atoms);
for at_typ=1:N_types
typ = Atoms(at_typ).typ;
%
for i=1:length(elem.data)
    if strcmp(typ,elem.textdata{i})
        index=i;
        break
    end
end
xyz=Atoms(at_typ).coord;
%retrieve corresponding information from elements.csv
rsize=elem.data(index,5);
natoms = size(xyz,1);

    
%%%%
    
%%%%
%% scan all points for atom most removed from the domain center
%%%%
    for at1 = 1: natoms
        xx=xyz(at1,3);
        yy=xyz(at1,2);
        zz=xyz(at1,1);
        rdis=sqrt(xx^2+yy^2+zz^2);
        rs=rdis+rsize;
        if rs>rmax
            rmax=rs;
        end	 
    end
end
rsize;
sph_rad=rmax;


fid = fopen('./rsdft.out', 'w');       %% create a output file
grid_vol = h^3;
nx = fix(2*sph_rad/h) + 1;              %%  Make sure nx is even
nx = 2*round(nx/2+0.01);
sph_rad=0.5*h*(nx-1);                  %%  Re-adjust R
ny = nx;                               %% [-radis, radius] in each direction
nz = nx;
Domain = struct('radius', sph_rad,'nx',nx,'ny',ny,'nz',nz,'h', h);
fprintf('  Grid points distance = %d\n', h);
fprintf('  Grid points on each dimension = %d\n', nx);


fprintf(fid, ' Number of states: \t%d\n\n', num_states);
fprintf(fid,  'Atom data:\n -------------\n');
fprintf(fid, ' Total # of atom types is %d\n', at);
atom_count = 0;
Atoms(1).coord ;
for atm_typ = 1:at
    xyz = Atoms(atm_typ).coord ;
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
%
disp('      ')
disp('      ')
disp('******************')
disp('     OUTPUT       ')
disp('******************')
disp('      ')
disp(' Working.....constructing Laplacian matrix...')

%% construct Laplacian operator
Laplacian = (1/(h*h))*fd3d(nx, ny, nz, fd_order);

nsizeA = size(Laplacian);
n      = size(Laplacian,1);
Hpot   = zeros(n,1);
pot    = Hpot;
index0 = length(pot);
%% 
%%
err = 1 + tol;
its = 0;
Ry  = 13.605698066;
%%
%%  Set up ion-core potentials/initial screening
%%  Import screening from the atom, i.e., hartree potential and charge.
%%
%%  Ion-Ion repulsion
[E_nuc0]=nuclear(Domain, Atoms);
%disp('Nuclear repulsion in eV')
%disp(E_nuc0*Ry)
%%
%%--------------------  Local part of the pseudopotential
%%
disp(' Working.....setting up ionic potential...')
[rho0, Hpot0, Ppot]  = pseudoDiag(Domain, Atoms);
%%%%%%%
hpsum0=sum(rho0.*Hpot0)*Ry;
fprintf(fid,' Initial Hartree energy (eV) = %10.5f  \n', hpsum0) ;
%%-------------------- count # atoms for stats
n_atoms  = 0;
for ii = 1:length(Atoms) 
    n_atoms = n_atoms + length(Atoms(ii).coord(:,1));
end
%%
%%--------------------  Non-local part of the pseudopotential
%%
disp(' Working.....setting up nonlocal part of ionic potential...')
[vnl] = pseudoNL(Domain, Atoms) ;  
%%
%%--------------------  Find a initial guess for the screening potential
%%--------------------  Use screening from Gaussian density
%%
indx1=length(rho0);
h=Domain.h;
%%
rhoxc = rho0' ./ grid_vol;
% Report max and min values of the charge density
dmax = max(rhoxc(:));
dmin = min(rhoxc(:));
fprintf(fid,' max and min values of charge density [e/bohr^3]');
fprintf(fid,'   %10.5e   %10.5e  \n', dmax, dmin);
if (dmin < 0) 
    fprintf(1,'warning in excorr.f: \n');
    fprintf(1,'NEGATIVE CHARGE DENSITY FOUND .\n');
end

zero_rhoxc = (rhoxc == 0);
[xcpot, exc] = XC_Update(rhoxc, grid_vol);
xcpot(zero_rhoxc) = 0.0;
Nelec = nelectrons(Atoms);
%%
%%------------------- open output file (wfn.dat)
%%
wfnid = fopen('./wfn.dat', 'wb');
pot = Ppot + Hpot0 + 0.5*xcpot;
%%-------------------- SCF LOOP
%%-------------------- when 'preconditioning' is used call ilu0
if (CG_prec) 
    disp('calling ilu0 ...') 
    %[L, U] = luinc(Laplacian,'0');
    [L, U] = ilu0(Laplacian);
    disp(' done.') 
    PRE = struct('L',L, 'U',U);
end

%%-------------------- clear persistent variables in mixer.
clear mixer;
%%-------------------- SCF LOOP starts here
scf_time_start = tic;
fprintf(fid, '\n----------------------------------\n\n');
while (err > tol && its <= maxits) 
    its = its+1;
    fprintf(1,'  Working ... SCF iter # %d  ... ',its);
    
    % Update Hamiltonian    
    B = 0.5 * Laplacian + spdiags(pot, 0, n, n) + vnl;
    
    % % Diagonalize Hamiltonian using the option defined in "include.m"
    tic;
    if (diagmeth == 1 || (its == 1 && diagmeth == 0)) 
        disp('calling lanczos..') 
        v = randn(n,1); 
        [psi, ev] = lan(B, n_occ+15, v, n_occ+500, 1.e-05);
    elseif (its == 1 && diagmeth == 2)
        disp('Calling SCF_Step1_CheFSI..')
        [psi, ev] = SCF_Step1_CheFSI(B, n_occ, poldeg);
    else 
        disp('Calling CheFSI..')
        [psi, ev] = CheFSI(B, psi, poldeg, ev);
    end
    diag_time = toc;
    fprintf(fid,' \n \n SCF iter # %d  ... \n',its);
    fprintf(fid,' Diagonalization time [sec] :\t%f\n\n', diag_time);
    
    % Get occupation factors and Fermi-level
    % Increase Fermi temp if does not converge
    [Fermi_level, occup] = occupations(ev(1:n_occ), Fermi_temp, Nelec, 0.000001);
    fprintf(fid, '   State  Eigenvalue [Ry]     Eigenvalue [eV]\n\n');
    for i = 1:n_occ
        eig = ev(i) * 2 * Ry;
        ry  = eig / Ry;
        fprintf(fid, '%5d   %15.10f   %18.10f  %5.2f\n', i, ry, eig, occup(i)); 
    end
    
    % Get charge density
    rho  = (psi(:,1:n_occ) .* psi(:,1:n_occ)) * 2 * occup;
    hrhs = (4*pi/grid_vol) .* (rho - rho0);
    rho  = rho ./ grid_vol;
    
    tic;
    if (CG_prec) 
        Hpot = pcg(Laplacian, hrhs, Hpot, 200, 1.e-04, PRE, 'precLU');
    else 
        Hpot = pcg(Laplacian, hrhs, Hpot, 200, 1.e-04);
    end
    hart_time = toc;
    fprintf(fid, '\nHartree potential time [sec]: \t%f\n\n', hart_time);
    
    % Report max and min values of the charge density
    dmax = max(rho(:));
    dmin = min(rho(:));
    fprintf(fid,' max and min values of charge density [e/bohr^3]');
    fprintf(fid,'   %10.5e   %10.5e  \n', dmax, dmin);
    if (dmin < 0) 
        fprintf(1,'warning in excorr.f: \n');
        fprintf(1,'NEGATIVE CHARGE DENSITY FOUND .\n');
    end
    
    [XCpot, exc] = XC_Update(rho, grid_vol);
    potNew = Ppot + 0.5 * XCpot + Hpot + Hpot0;
    err = norm(potNew - pot) / norm(potNew);
    fprintf(fid,'   ... SCF error = %10.2e  \n', err) ;
    fprintf(1,'   ... SCF error = %10.2e  \n', err) ;

    % Mixer
    pot = mixer(pot, potNew-pot);
end
scf_time = toc(scf_time_start);
fprintf('##### Total SCF time = %.2f (s)\n', scf_time);

disp('          ')
disp('**************************')
disp(' CONVERGED SOLUTION!! ')
disp('**************************')
disp('         ')
fprintf(1, '   State  Eigenvalue [Ry]     Eigenvalue [eV]  Occupation \n');
for i = 1:n_occ
    eig = ev(i) * 2*Ry;
    ry = eig / Ry; 
    occ=2*occup(i);
    fprintf(1, '%5d   %15.4f   %18.3f  %10.2f\n', i, ry, eig, occ); 
end
%%  
%%-------------------- total energy 
%%  Sum over eigenvalues--Put everything in Ryd
Esum= sum(ev(1:n_occ).*occup(1:n_occ)); 
Esum0=4*Esum;
%%
%%--------------------   Hartree potential
%%  Factor of two for double counting--converts to Ryd
Hsum0=sum(rho.*(Hpot+Hpot0))*grid_vol;
%%
%%-------------------- Exchange correlaion sum
%% No factor of two because energy is in Ry
%%
Vxcsum0=sum(rho.*XCpot)*grid_vol;
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
fprintf(fid,' Total electronic energy = %10.5f  eV   = ',E_total) ;
fprintf(fid,'  %10.5f  Ry  \n',E_total0) ;
fprintf(fid,' Electronic energy/atom  = %10.5f  eV   = ',E_total/n_atoms) ;
fprintf(fid,'  %10.5f  Ry  \n',E_total0/n_atoms) ;
  
%%-------------------- display energies..
fprintf(1,'  \n');
fprintf(1,' Total Energies \n\n');
fprintf(1,' Sum of eigenvalues      = %10.3f  eV   = ',Esum) ;
fprintf(1,'  %10.4f  Ry  \n',Esum/Ry) ;
fprintf(1,' Hartree energy          = %10.3f  eV   = ',Hsum) ;
fprintf(1,'  %10.4f  Ry  \n',Hsum/Ry) ;
fprintf(1,' Exchange-corr. energy   = %10.3f  eV   = ',Excsum) ;
fprintf(1,'  %10.4f  Ry  \n',Excsum/Ry) ;
fprintf(1,' Ion-ion repulsion       = %10.3f  eV   = ',E_nuc) ;
fprintf(1,'  %10.4f  Ry  \n',E_nuc/Ry) ;
fprintf(1,' Total electronic energy = %10.3f  eV   = ',E_total) ;
fprintf(1,'  %10.4f  Ry  \n',E_total0) ;
fprintf(1,' Electronic energy/atom  = %10.3f  eV   = ',E_total/n_atoms) ;
fprintf(1,'  %10.4f  Ry  \n',E_total0/n_atoms) ; 
  
  
  
%%

%%-------------------- free memory (persistent variables) in mixer.
clear mixer;
%%------------------- Output results
fprintf(fid, '\n Finished\n');

little_big_test = 26;
fwrite(wfnid, little_big_test, 'uint32');

fwrite(wfnid, Domain.radius, 'double');
fwrite(wfnid, Domain.h, 'double');

pot_length = length(pot(:,1));
fwrite(wfnid, pot_length, 'uint32');
fwrite(wfnid, pot, 'double');

rho_length = length(rho(:,1));
fwrite(wfnid, rho_length, 'uint32');
fwrite(wfnid, rho, 'double');

w_length = length(psi(:,1));
w_col_length = length(psi(1,:));
fwrite(wfnid, w_length, 'uint32' );
fwrite(wfnid, n_occ, 'uint32');
for i = 1:n_occ
    fwrite(wfnid, psi(:, i), 'double');
end

fclose(wfnid);
fclose(fid); %% Close output file


