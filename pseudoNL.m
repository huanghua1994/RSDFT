function [vnl]=pseudoNL(Domain, Atoms)
% nlpseudo3 computes the nonlocal pseudopotential contribution to potential energy
% This function is called by rsdft.m
%
% Input Parameters
% Domain - Domain is created in main.m and basically describes the spherical area in 
%          which the nonlocal pseudopotential as well as the local pseudopotential
%          should be analyzed.
% Atoms  - Atoms is created in parsec_dat.m and contains information related to the
%          different atoms of interest.
%
% Output Parameters
% vnl    - a sparse matrix of size Domain.nx*Domain.ny*Domain.nz, where it is assumed
%          that Domain.nx  == Domain.ny == Domain.nz and theoretically has at most
%          Domain.nx^4 non-zero elements

%% Creation of local variables for the variables contained within the Domain struct
% Assume cubic domain. This assumption is valid only because of how Domain is created
% in main.m, since the values for ny and nz are set to be nx.
nx      = Domain.nx;
ny      = Domain.ny;
nz      = Domain.nz;
h       = Domain.h;
rad     = Domain.radius;

% load the spline data
splineData;

% ndim is the dimension of the sparse matrix for the nonlocal pseudopotential (vnl)
ndim    = nx * ny * nz;
% imax is the default maximum nonzero values that should be found in vnl
%%%%----FIX-ME
imax    = optimalSize(Atoms, h);

%%%%%-imax    = ndim*50;

vnl     = spalloc(ndim,ndim,imax);

% This is the number of atoms that need to be looked at
N_types = length(Atoms);
% This is just importing the information about each element
elem    = importdata('elements_new.csv');

% Find out what columns pot_P, pot_S, and wav_P are stored in
pot_P = 0; 
pot_S = 0; 
wav_P = 0;

for i = 1:length(data_list)
  if strcmp('pot_P',data_list(i))
    pot_P = i;
  end
  if strcmp('pot_S',data_list(i))
    pot_S = i;
  end
  if strcmp('wfn_P',data_list(i))
    wav_P = i;
  end
end

%%-------------------- Loop over atom types
typ=Atoms(1).typ;

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
for at_typ=1:N_types
  % This gets the first atom and assigns the atomic symbol to typ
  typ = Atoms(at_typ).typ;
  % This iterates through elem looking for the matching element data

%%% Start of generating the coefficients for splines
  for i = 1: length(AtomFuncData)
    if strcmp(typ,AtomFuncData(i).atom)
      index=i;
      break
    end
  end
%%%% We assume that the first column is for the radius
  xi = AtomFuncData(index).data(:,1);
%%%% Now we need to find the index for the P wave function
%%%% Now we do the same thing for finding the P and S potentials
  [rows cols] = size(AtomFuncData(index).data);

  if cols < pot_P
    potPS = zeros(rows,1);
  else
    potPS = AtomFuncData(index).data(:,pot_P);
  end

  if cols >= pot_S
    potPS = potPS - AtomFuncData(index).data(:,pot_S);
  end
%
  if cols < wav_P
    wfn_P = zeros(rows,1);
  else
    wfn_P = AtomFuncData(index).data(:,wav_P);
  end
  
%%-------------------- pre-process the data
   I = preProcess(wfn_P);
   wfn_P = wfn_P(I);
   xi_wfn_P = xi(I);
   
   I = preProcess(potPS);
   potPS = potPS(I);
   xi_potPS = xi(I);

  [zWav, cWav, dWav] = fspline(xi_wfn_P,wfn_P);
  [zPotPS, cPotPS, dPotPS] = fspline(xi_potPS,potPS);
%%--------------------  End generating the coefficients for splines

  for i = 1:length(elem.data)
    if strcmp(typ,elem.textdata{i})
      index=i;
      break
    end
  end
% This grabs the Z value from elements.csv
% Z=elem.data(index,3);
% xint is equal to A divided by the cube of the grid spacing
  xint=elem.data(index,3)/h^3;

  %%%%%%% DEBUGGING INFORMATION %%%%%%%
%  disp('Element')
%  disp(typ)
%  disp(' Z ')
%  disp(Z)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % xyz corresponds to the atom's position within the Domain
  xyz=Atoms(at_typ).coord;
  % natoms is just the number of typ atoms that are in the system
  natoms = size(xyz,1);
  % Rzero is the core size of the atom

%%%% FIX ME -- THIS SHOULD COMPUTED -- RADIUS WHERE
%%%% R
  Rzero=elem.data(index,6);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This for loop iterates through all the atoms of at_typ
  for at = 1 : natoms
%-------------------- xxa, yya, zza, store x,y,z positions of atom
    xxa=xyz(at,1);
    yya=xyz(at,2);
    zza=xyz(at,3);
%-------------------- Loop over grid points near atom
    i0=round((xxa+rad)/h+1);
    j0=round((yya+rad)/h+1);
    k0=round((zza+rad)/h+1);
%-------------------- Determine the maximum distance from initial (x, y, z)
    span=round(Rzero/h);   %%%-

    indx=0;

    for k=k0-span:k0+span
      zzz = (k-1)*h - rad - zza;
      for j=j0-span:j0+span
        yyy = (j-1)*h - rad - yya;
        j_p_ps = 1; j_wfn = 1;
        for i=i0-span:i0+span
          xxx = (i-1)*h - rad - xxa;
%-------------------- calculate distance from (xxx,yyy,zzz) to atom
          dd2=xxx^2+yyy^2+zzz^2;
          dd1=sqrt(dd2);
% The checking order was changed to take advantage of short circuiting.
% check if the distance is within the core size.
          if (dd1 < Rzero && dd1 > 0)
            indx=indx+1;

            nn(indx)=i + nx * ((j-1) + (k-1)*ny);
% keep track of the x, y, z coordinates and the distance
            xx(indx)=xxx;
            yy(indx)=yyy;
            zz(indx)=zzz;
            dd(indx)=dd1;
            [ vspp(indx)  j_p_ps ] = fsplevalIO(zPotPS,cPotPS,dPotPS,xi_potPS,dd1, j_p_ps);
            [ wavpp(indx) j_wfn  ] = fsplevalIO(zWav,cWav,dWav,xi_wfn_P,dd1, j_wfn);
          end % end of core size if statement
        end % end of for loop concerning xxx
      end % end of for loop concerning yyy
    end % end of for loop concerning zzz
%
    %%%%%%% DEBUGGING INFORMATION %%%%%%%
%   disp('Z xint')
%   disp(Z)
%   disp(xint)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Optimize by pre-allocation
% Here we're creating a imx x 3 matrix of zeroes
% The natoms * portion was removed because the code currently never runs into a
% situation where more than indx^2 + 1 space is needed.
% ??? indx+1 ??
    imx=indx^2 + 1;
    vnll=zeros(imx,3);
% Moved the work from inside the loops to outside the loops. This should improve run
% time since MATLAB is built specifically to perform Matrix operations, which was 
% being done cell by cell within the for loops.
% This is an explanation of what is being done within the following six lines of code.
% { xx, yy, zz } ./ dd calculates px, py, and pz respectively element by element.
% wavpp .* vspp calculates the fac for ulmspx2, ulmspy2, and ulmspz2.
% The operation of { px, py, pz } .* fac results in the values for ulmsp?2.
% The difference between ulmsp?1 and ulmsp?2 is a factor of ./ xint.
% Because of this, ulmsp?2 is calculated first and then ulmsp?1 is calculated by taking
% ulmsp?2 and dividing by ./ xint.
% This should greatly reduce the number of actual calculations that are performed and
% have them be calculated more optimally.
%
    ulmspx = xx ./ dd .* wavpp .* vspp;
    ulmspy = yy ./ dd .* wavpp .* vspp;
    ulmspz = zz ./ dd .* wavpp .* vspp;
    xmatrix = (ulmspx' * ulmspx) ./ xint;
    ymatrix = (ulmspy' * ulmspy) ./ xint;
    zmatrix = (ulmspz' * ulmspz) ./ xint;
    total = xmatrix + ymatrix + zmatrix;
%
    n_index = 0;
%-------------------- Double loop over r and r'
    for i = 1 : indx
      %   index for r
      for j = 1 : indx
        %% Put the two together
        % The following code has been altered from the original to work with the new
        % method of calculating the ulmsp?? outside of the for loops.
        n_index = n_index + 1;
        vnll(n_index, 1) = nn(i);
        vnll(n_index, 2) = nn(j);
        vnll(n_index, 3) = total(i,j);
      end % end of index for r for loop
    end % end of double for loop over r and r'

%   disp(' End atom loop')

% This adds the size information so that spconvert can be used
    vnll(n_index+1,1)=ndim ;
    vnll(n_index+1,2)=ndim ;
    vnll(n_index+1,3)=0 ;
%%######
    vnl=vnl + spconvert(vnll);
  end % end of for loop of a single at_typ
end % end of at_typs for loop

end % end of function
