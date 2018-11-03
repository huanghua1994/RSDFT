close all

%% This file allows to plot the atomic positions in the rsdft call.

%% Load data 

%parsec_dat_bis;
parsec_dat;
%benzene_parsec_dat;
%Si_parsec_dat;
%Si_bfgs_parsec_dat;
  
%% Import the information about each element.
%% Particularly, element.csv allows to obtain the atomic mass and the 
%% atomic radius of each atom.
% The name of this file has to be changed for elements.csv if all is ok
elem    = importdata('elements.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%			   COLORMAP			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This creates the colormap for the plot
% Each main atom (C,O,H,Na,Si,Mg)will be ploted in a particular color
% Each other atoms will be ploted in a particular gray color which depends  
% on its atomic mass
nn  = size(elem.data,1);
map = zeros(nn,3);
mass_max = max(elem.data(:,10));
	
for i = 1:nn
	if strcmp('H',elem.textdata{i})
	   map(i,:)  =[205/255 205/255 193/255];
	else if strcmp('C',elem.textdata{i})
	     	map(i,:)  =[0 0 205/255];
	     else if strcmp('O',elem.textdata{i})
	     	     map(i,:)  =[205/255 0 0];
	          else if strcmp('Na',elem.textdata{i})
		          map(i,:)  =[138/255 43/255 226/255];
		       else if strcmp('Si',elem.textdata{i})
		       	       map(i,:)  =[210/255 180/255 140/255];
			    else if strcmp('Mg',elem.textdata{i})
				     map(i,:)  =[102/255 205/255 0];
				 else	
				 	%cc = (elem.data(i,10)/mass_max);
					cc = 1-((i-1)/nn);
					map(i,:)  =[cc cc cc];
				 end
			    end
		       end
		  end
	     end
	end
end  
%colormap(map)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%		      DISPLAY GRAPHICS			     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
rad = sph_rad;
fig = figure('Name','Atomic Positions');
%title('Atomic Positions','FontAngle','italic','FontWeight','bold','FontSize',16);					
box = [-rad rad -rad rad -rad rad];
axis(box) 
hold
caption0 = [];
grid on
H = zeros(1,at);
M = [];

%% these variables allow to manage the writing of the atoms name
% sol_1 is used with less than 4 atoms --- use the barycentre
% sol_2 (= not sol_1) is used with more than 4 atoms --- simple method
n_atoms = 0;	% Number of atoms
for k = 1:at
	n_atoms = n_atoms + size(Atoms(k).coord,1);
end

sol_1 = (n_atoms <= 4);

if (sol_1)
	xx = 0;
	yy = 0;
	zz = 0;

	for k = 1:at
		for j = 1:size(Atoms(k).coord,1)
			xx = xx+Atoms(k).coord(j,1);
			yy = yy+Atoms(k).coord(j,2);
			zz = zz+Atoms(k).coord(j,3);
		end 
	end
	bary = [xx/n_atoms yy/n_atoms zz/n_atoms];
end

for atm_typ = 1:at

	% find the type of the considered atom
	typ = Atoms(atm_typ).typ;
	% This iterates through elem looking for the matching element data
	for i = 1:nn
		if strcmp(typ,elem.textdata{i})
			index = i;
			break
		end
	end

	%load the coordinates, the atomic mass and the atomic radius of the atom
  	xyz    = Atoms(atm_typ).coord;
	mass   = elem.data(index,10);
	radius = elem.data(index,11);
	
	%This creates handles and matching names for the legend
	H(atm_typ) = plot3(xyz(1,1),xyz(1,2),xyz(1,3),'o','color',map(index,:),'linewidth',radius/12,'Visible','off');
	M = [M {elem.textdata{index}}];
	
	%plot the atom, the linewidth is proportionnal to its radius,
	% and the color is taken from the colormap previously created,
	for k = 1:size(xyz,1)
  		plot3(xyz(k,1),xyz(k,2),xyz(k,3),'o','color',map(index,:),'linewidth',radius/4); 
		hold on
	
	% this adds above of the atom its symbol
	if (~ sol_1)
		Xr = xyz(k,1);
		Yr = xyz(k,2);
		Zr = xyz(k,3)+0.3;		
	else
		c = 0.5;
		Xr = (xyz(k,1)-bary(1))*c + xyz(k,1);
		Yr = (xyz(k,2)-bary(2))*c + xyz(k,2);
		Zr = (xyz(k,3)-bary(3))*c + xyz(k,3);	
	end
	text(Xr,Yr,Zr,elem.textdata{index},'FontSize',12)
		
	%connect this atom to the other
	
	%% this sets the linewidth function of the number of atoms -- usually less than 12 atoms...
	if (n_atoms<=3)
		lw = 5;
		else if (n_atoms<=10)
			lw = 3;
			else
				lw = 2;
		     end
	end
	
  	for j = atm_typ:at
		for i = 1:size(Atoms(j).coord,1)
			x_coord = [xyz(k,1) Atoms(j).coord(i,1)];
			y_coord = [xyz(k,2) Atoms(j).coord(i,2)];
			z_coord = [xyz(k,3) Atoms(j).coord(i,3)];
			plot3(x_coord,y_coord,z_coord,'color',[0.5 0.5 0.5],'linewidth',lw)
  		end
	end
	end
end

axis off 
axis equal
%axis vis3d
set(gcf, 'color', [245/255 245/255 245/255]);
set(gcf, 'InvertHardCopy', 'off');

legend(H,M)

%cursor
dcm_obj = datacursormode(fig);
set(dcm_obj,'DisplayStyle','window',...
'SnapToDataVertex','off','Enable','on')
info_struct = getCursorInfo(dcm_obj);

% Rotation of the 3D-graph
for i=1:300
camorbit(1,1,'camera')
drawnow
end

		
