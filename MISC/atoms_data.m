%% This is the configuration file for generating a shorter 'data_atoms' table 
%% than the currently .txt files, in order to be used directly by the splines.

fid0 = fopen('./splineData.m', 'w');
fprintf(fid0, '%% Spline data for atoms\n\n');
fprintf(fid0, 'nAtoms = 0;\n');

%% Number of points which are going to be kept
nPts = 70;

%% load and read the atoms_list.txt file, allows to add easily other atoms 
fid = fopen('./atoms_list.txt','r');
l = fgetl(fid);
i = 1;

%% "data_list" lists the names of the data gathered in the generated table
%% The data will appear in the order of the data_list
data_list = [{'radius'}];

while l~=-1
      atom = l
% loading the data files from the directory of the atom
%      files = dir(['./PSEUDO/' atom '_atom/*.dat']);
      files = dir(['./' atom '_atom/*.dat']);
      nn = size(files,1);
%  This creates a list with all the name of data files "files_name"
%  and updates data_list
%  To do it working we have to remove 3 data files: atom.dat,
%  plot.dat and pseudo.dat
	
	files_name = [];
	for i = 1:nn
		[pathstr, name, ext, versn] = fileparts(files(i,1).name);
   if (~(strcmp(name,'atom') || strcmp(name,'plot') || strcmp(name,'pseudo')))
			files_name = [files_name;{name}];
			idx = strfind(data_list,name);
			if (sum([idx{:,:}]) == 0)
				data_list = [data_list;{name}];
			end
       	end
    end
	
%% This creates the matching paths and load the data files
%% All is initialized considering the first case as a particular case
	n_files = size(files_name,1);
%	p_path = ['./PSEUDO/' atom '_atom/' files_name{1} '.dat'];
	p_path = ['./' atom '_atom/' files_name{1} '.dat'];
	A = load(p_path);
	n = length(A);
	
	for i = 2:n_files
%%		p_path = ['./PSEUDO/' atom '_atom/' files_name{i} '.dat'];
		p_path = ['./' atom '_atom/' files_name{i} '.dat'];
		B = load(p_path);
		if (length(B)<n)
			n = length(B);
		end
	end
%%-------------------- RADIUS
	h = 1:fix(n/(nPts-1)):n;
	l_h = length(h);
	radius = A(h,1);
	
	M = zeros(l_h,n_files);
	
	for i = 1:n_files
%%		p_path = ['./PSEUDO/' atom '_atom/' files_name{i} '.dat'];
		p_path = ['./' atom '_atom/' files_name{i} '.dat'];
		B = load(p_path);
		M(:,i) = B(h,2);
		
		%% scales the atoms data
%% VARIOUS SCALINGS
	  if (strcmp(files_name(i),'charge'))
			M(:,i) = M(:,i)./(radius.^2);
               elseif (strcmp(files_name(i),'wfn_S') ) 
                        M(:,i) = 0.5* (M(:,i) ./ radius);
	       elseif (strcmp(files_name(i),'wfn_P')) 
                        cons = sqrt(0.75/pi);
		        M(:,i) = cons*(M(:,i) ./ radius);
	       else
    			M(:,i) = 0.5*M(:,i);
               end
        end
%% write data in splineData.m file
fprintf(fid0,'%% %s data\n', atom);
fprintf(fid0,'nAtoms = nAtoms + 1;\n');
fprintf(fid0,'atom = ''%s'';\n',atom);
fprintf(fid0,'data = [...\n ');
%find the data name in the data_list
	n_data = size(data_list,1);
	index = zeros(1,n_data);
	k = 1;
for j = 2:n_data
	for i = 1:n_files
		if (strcmp(data_list(j),files_name(i)))
			index(j) = i;
			break;
		end
		% else the data_list(j) is not a data for this atom... 
		% we have to put an empty row for this data
		index(j) = n_files+k;
		k = k+1;
	end
end
  for j = 1:l_h
  	fprintf(fid0,'  %.4e',radius(j));
  	for k=2:n_data
  		if (index(k)<=n_files)
    			fprintf(fid0,'  %.4e',M(j,index(k)));
    		else
			fprintf(fid0,'   %.4e',NaN);
		end    
    	end
	fprintf(fid0,'  ;...\n '); 
  end
fprintf(fid0,'  ];\n');
fprintf(fid0,'str = struct(''atom'', atom, ''data'', data);\n');
fprintf(fid0,'AtomFuncData(nAtoms) = str;\n\n');
	
%% go to the next atom
	i = i+1;
	l = fgetl(fid);
end
fprintf(fid0, '%% List of data name:\n');
fprintf(fid0, 'data_list = [ ');
for i=1:n_data-1
	fprintf(fid0,'{''%s''}, ',data_list{i});
end
fprintf(fid0, '{''%s''} ];',data_list{n_data});


fclose(fid);
fclose(fid0);


