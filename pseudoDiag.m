function [rho0, hpot0, pot]=pseudoDiag(Domain, Atoms);
%
% Set up initial screening and local ionic pseudopotential
%  
%%%  Initial charge density is from the atom.  
%%%  Atomic Hartree potential used as reference.  It is 
%%%  added and subtracted from the total energy.
%%%  Exchange correlation is determined from a superposition 
%%%  of atomic charge
%%%

% parsec_dat;

%%% Creating the splines:
%atoms_data;
splineData;

%%% Opening file -- Right now, this is not used
%fid=fopen('oxy_s_pot.txt','r'); 
%pot_oxy=fscanf(fid,'%g %g', [2 inf]);
%rr_oxy_s=pot_oxy(1,:);
%pot_oxy_s=pot_oxy(2,:);
%n_oxy=length(rr_oxy_s)
%rc_oxy_s=rr_oxy_s(n_oxy)   %% last (largest) value or r 

%%% Localizing variables
nx   = Domain.nx;
ny   = Domain.ny;
nz   = Domain.nz;
h    = Domain.h;
rad  = Domain.radius;

%%% Pre-allocating memory
ndim = nx*ny*nz;
pot  = zeros(ndim,1);
rho0 = zeros(ndim,1);
hpot0 =zeros(ndim,1);

%%%%
N_types = length(Atoms);
Z_sum=0.;

%%---------------- read elements.csv -- info on atoms
elem=importdata('elements_new.csv') ;

%-------------------- Atom-type loop
for at_typ=1:N_types
    typ = Atoms(at_typ).typ;
%	disp('Element')
%    disp(typ)
    for i=1:length(elem.data)
        if strcmp(typ,elem.textdata{i})
            index=i;
            break
        end
    end
    xyz=Atoms(at_typ).coord;
%-- retrieving corresponding information from elements.csv
    Z=elem.data(index,2);
    natoms = size(xyz,1);
    Z_sum=Z_sum+Z*natoms;
    Rzero=elem.data(index,6);
    trad2=Rzero*Rzero;
                      
%%--Searching for atom's data
    for i=1:length(AtomFuncData)
        if strcmp(typ,AtomFuncData(i).atom)
            index=i;
            break;
        end
    end
    
%%--Localizing variables and initializing arrays
% Find out what columns charge,pot_S and hartree are stored in
% note that we assume that the radius is in column 1
    for k = 1:length(data_list)
    	if strcmp('charge',data_list(k))
    		i_charge = k;
  	end
  	if strcmp('pot_S',data_list(k))
    		i_pot_S = k;
  	end
  	if strcmp('hartree',data_list(k))
    		i_hartree = k;
  	end
    end

    atom_data=AtomFuncData(i).data;
    x_charg=atom_data(:,1);
    y_charg=atom_data(:,i_charge);
    x_pot_s=atom_data(:,1);
    y_pot_s=atom_data(:,i_pot_S);
    x_vhart=atom_data(:,1);
    y_vhart=atom_data(:,i_hartree);
    
%%-------------------- pre-processing the data
    I = preProcess(y_charg);
    x_charg = x_charg(I);
    y_charg = y_charg(I);
        
    I = preProcess(y_pot_s);
    x_pot_s = x_pot_s(I);
    y_pot_s = y_pot_s(I);
        
    I = preProcess(y_vhart);
    x_vhart = x_vhart(I);
    y_vhart = y_vhart(I);
    
%%--Calculating the splines
    [z_chg,c_chg,d_chg]=fspline(x_charg, y_charg);
    [z_p_s,c_p_s,d_p_s]=fspline(x_pot_s, y_pot_s);
    [z_vht,c_vht,d_vht]=fspline(x_vhart, y_vhart);
    
%-- Scan all points for each atom (not optimal but OK)--copied from
%-- ppot.m with a small change of formula to find potential
    for at=1:natoms
        indx = 0;
        for k=1:nz
%-- atom coordinates are given wrt sphere center.. need to adjust
%-- x_point - x_A = (x_point - x_center) - (x_A -x_center) 
            dz = ((k-1)*h - rad - xyz(at,3))^2;
            for j=1:ny
                dy = dz+((j-1)*h - rad  - xyz(at,2))^2;
		
		%%--------initialization of intervals j_ch;j_p_s and j_vht
		j_ch  = 1;
		j_p_s = 1;
		j_vht = 1;
%%		
                for i=1:nx
                    r2 = dy+((i-1)*h - rad  - xyz(at,1))^2;
                    r1=sqrt(r2);
%%----------------- next point in space				 
                    indx = indx+1;

%%----------------- Evaluating the splines
                    [ppot j_p_s]   = fsplevalIO(z_p_s,c_p_s,d_p_s,x_pot_s,r1,j_p_s);
                    [rrho j_ch]    = fsplevalIO(z_chg,c_chg,d_chg,x_charg,r1,j_ch);
                    %rrho = rrho*h^3/(4.*pi);
		    rrho = max(0,rrho*h^3/(4.*pi));
		    [hpot00 j_vht] = fsplevalIO(z_vht,c_vht,d_vht,x_vhart,r1,j_vht);
		    %hpot00 = hpot00*h^3/(4.*pi);                    
                
%%----------------- done atom-specific calculations - now compute
%%----------------- potentials, charge.
                    pot(indx)   = pot(indx)+ppot;
		    rho0(indx)  = rho0(indx)+rrho;
                    hpot0(indx) = hpot0(indx)+hpot00;
                    
                end  %% end x
            end      %% end y
        end          %% end z
    end              %% end atom
end                  %% end atom_typ

%%---- Check charge
%disp('charge total')
Rsum0=sum(rho0);
sum_charg=Rsum0;
%disp(sum_charg)
rho0=Z_sum*rho0/Rsum0;
sum_charg=sum(rho0);
%%disp(' renormalized charge')
%%disp(sum_charg)

