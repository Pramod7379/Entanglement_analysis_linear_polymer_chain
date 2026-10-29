% At present, this script works for monodisperse polymer chains containing n atoms each. 
% Please change the value of n in 58th line for your sample. 
% In this script, we take input from lammps dump file ('revcon.*' in 
% this case; you have to replace the file nam in line 24) 
% at different points of time in the course of deformation. 
% We use the dump command
% dump       id all custom 10000 rev_final.* id x y z vx vy vz fx fy fz
% 
% This script sorts out different chains and their distinct images and saves
% output data in 'output_mat_file' directory. Some of these data are coordinates of distinct
% images (coord_m_full), number of distinct images for each
% chain(diffseg_chain) their cumulative sum (diffcumseg_chain) and dimension
% of simulation box(slx, sly, slz).



main_dir = pwd;  % get working directory
for num = 1:1 %num is time steps

    clearvars -except num main_dir
    
    cd(strcat(main_dir,'/../simulation_data'))
%     name = strcat('revcon.',num2str(num)); % MD dump file from LAMMPS
    name = 'rev_final.1000000'; % change filename here
    fid = fopen(name);
    for i = 1:5
        fgetl(fid);
    end
    
    l=[];
    for i = 1:3
        l = [l;sscanf(fgetl(fid),'%f')];
    end
    
    slx = l(2) - l(1);
    sly = l(4) - l(3);
    slz = l(6) - l(5);
    fclose(fid);
    
    fid = fopen(name);
    for i = 1:9
        fgetl(fid);
    end
    buffer = fread(fid,Inf);
    fclose(fid);
    
    mkdir('data_only')
    cd('data_only')
    fid = fopen(strcat('revcon_',num2str(num)),'w'); % this file contains only the coordinate of atoms by skipping all other information about simulation box
    fwrite(fid,buffer);
    fclose(fid);
    
    a=load(strcat('revcon_',num2str(num)));
    a=sortrows(a,1);
    ch=0;
    
    num_atom = size(a,1);
    n = 200; %number of atoms in one polymer chain; change according to your requirement
    N = num_atom/n;
   
    
    for ijk=1:n:(num_atom-n+1)
        ch=ch+1;chn(ch)=n;
        for j=1:n
            coord_array(j,1,ch)=a(ijk-1+j,2);
            coord_array(j,2,ch)=a(ijk-1+j,3);
            coord_array(j,3,ch)=a(ijk-1+j,4);
        end
    end
    
    % seg_chain(i)
    % i : 1:N
    seg_chain = zeros(N,1);
    seg_tot = 0;
    
    for i =1 :N
        % extract coordinates of all atoms in a chain
        x = coord_array(:,1,i);
        y = coord_array(:,2,i);
        z = coord_array(:,3,i);
        
        % find the vector distance between coordinates
        diff_x = diff(x);
        diff_y = diff(y);
        diff_z = diff(z);
        
        % identifying indices where the chain crosses the boundaries of the simulation
        % box
        ind_x_neg = find(diff_x < -slx/2);
        ind_x_pos = find(diff_x > slx/2);
        
        ind_y_neg = find(diff_y < -sly/2);
        ind_y_pos = find(diff_y > sly/2);
        
        ind_z_neg = find(diff_z < -slz/2);
        ind_z_pos = find(diff_z > slz/2);
        
        ind = unique(sort([0;ind_x_neg;ind_x_pos;ind_y_neg;ind_y_pos;ind_z_neg;ind_z_pos;n]));
        % number of segments of the chain. the number of segments is equal
        % to the no of box crossings + 1
        seg_chain(i) = length(ind)-1; % total number of discontinuous segments in a chain
        
        if i == 1
            segchain_ind = {ind};
        else
            segchain_ind = [segchain_ind;{ind}];
        end
        segchain_ind_o{i} = segchain_ind{i};
        
        % coordinates of every segment of a chain after projecting it
        % according to its order, i.e. increasing order of the indices from
        % 1 to n. every entry in ind marks the discontinuity of chain at
        % the boundary of simulation box. The for loop below does the same.
        
        for j = 1:length(ind)-1
        
            
            if i == 1 && j == 1
                
                coord_array_seg = {[coord_array(ind(j)+1:ind(j+1),1,i) coord_array(ind(j)+1:ind(j+1),2,i) coord_array(ind(j)+1:ind(j+1),3,i)]};
                
                temp=[i j]
                temp1=coord_array_seg;
                temp2=ind;
                if isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                    
                elseif isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{i} = [coord_array_seg{i};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                end
                temp4=coord_array_seg;
                
            else
                coord_array_seg = [coord_array_seg;{[coord_array(ind(j)+1:ind(j+1),1,i) coord_array(ind(j)+1:ind(j+1),2,i) coord_array(ind(j)+1:ind(j+1),3,i)]}];
                
                if isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)+slz]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0 && isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)-slz]];
                    
                elseif isempty(find(ind_x_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)+slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_x_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i)-slx coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_y_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)+sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_y_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i)-sly coord_array(ind(j+1)+1,3,i)]];
                    
                elseif isempty(find(ind_z_neg == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)+slz]];
                    
                    
                elseif isempty(find(ind_z_pos == ind(j+1), 1)) == 0
                    coord_array_seg{end} = [coord_array_seg{end};[coord_array(ind(j+1)+1,1,i) coord_array(ind(j+1)+1,2,i) coord_array(ind(j+1)+1,3,i)-slz]];
                    
                end
                
            end
        end
        
        seg_tot = seg_tot + length(ind)-1; % total number of segments in the ensemble
        
    end
    
    segchain_ind_o = segchain_ind_o';
    cumseg_chain = cumsum(seg_chain); % cumulative sum of no of segments in chain order. 
    len_seg = zeros(seg_tot,1); % stores length of every segment in the following loop.
    
    for i = 1:seg_tot
        
        len_seg(i) = size(coord_array_seg{i},1);
        
    end
    
    % Following if-else and for-loop skips all the segments that consist of 
    % either 0 or 1 atom as they will be accounted for by other segments.  
    
    if min(len_seg) == 0 || min(len_seg) == 1
        
        coord_array_seg_o = coord_array_seg;
        len_seg_o = len_seg;
        seg_chain_o = seg_chain;
        cumseg_chain_o = cumseg_chain;
        coord_array_seg = {};
        seg_skip = [];
        
        for i = 1:seg_tot 
            if len_seg(i) == 1 || len_seg(i) ==0
                ind_skip = find(cumseg_chain >= i, 1 );
                seg_chain(ind_skip) = seg_chain(ind_skip) - 1;
                seg_skip = [seg_skip;i];
                continue
            else
                if size(coord_array_seg,1) == 0
                    coord_array_seg = {coord_array_seg_o{i}};
                else
                    coord_array_seg = [coord_array_seg;{coord_array_seg_o{i}}];
                end
            end
        end
        
        cumseg_chain = cumsum(seg_chain);
        seg_tot = size(coord_array_seg,1);
        len_seg = zeros(seg_tot,1);
        
        for i = 1:seg_tot
            len_seg(i) = size(coord_array_seg{i},1);
        end
    end
    
    diff_seg_full = [];
    seg_same_full = [];
    
    % In next step, we sort out segments that are the part of different chain 
    % images and discard all segments that are the part of the same chain image      
    for i = 1:N
        
        temp_seg = cumseg_chain(i)-seg_chain(i)+1:cumseg_chain(i); 
        seg_same = [];
        diff_seg = [];
        
        for jk = 1:seg_chain(i)
            
            if isempty(find(seg_same == jk, 1)) == 0
                continue
            end
            
            diff_seg_full = [diff_seg_full;cumseg_chain(i)-seg_chain(i)+jk];
            diff_seg = [diff_seg;cumseg_chain(i)-seg_chain(i)+jk];
            
            coord_m = coord_array(:,:,i);
            
            xx_diff = diff(coord_array(segchain_ind{i}(jk)+1:n,1,i));
            yy_diff = diff(coord_array(segchain_ind{i}(jk)+1:n,2,i));
            zz_diff = diff(coord_array(segchain_ind{i}(jk)+1:n,3,i));
            
            ind_x_sub = find(xx_diff > slx/2) + segchain_ind{i}(jk);
            if isempty(ind_x_sub) == 0
                for j = 1:length(ind_x_sub)
                    coord_m(ind_x_sub(j) + 1:end,1) = coord_m(ind_x_sub(j) + 1:end,1) - slx;
                end
            end
            
            ind_x_add = find(xx_diff < -slx/2) + segchain_ind{i}(jk);
            if isempty(ind_x_add) == 0
                for j = 1:length(ind_x_add)
                    coord_m(ind_x_add(j) + 1:end,1) = coord_m(ind_x_add(j) + 1:end,1) + slx;
                end
            end
            
            ind_y_sub = find(yy_diff > sly/2) + segchain_ind{i}(jk);
            if isempty(ind_y_sub) == 0
                for j = 1:length(ind_y_sub)
                    coord_m(ind_y_sub(j) + 1:end,2) = coord_m(ind_y_sub(j) + 1:end,2) - sly;
                end
            end
            
            ind_y_add = find(yy_diff < -sly/2) + segchain_ind{i}(jk);
            if isempty(ind_y_add) == 0
                for j = 1:length(ind_y_add)
                    coord_m(ind_y_add(j) + 1:end,2) = coord_m(ind_y_add(j) + 1:end,2) + sly;
                end
            end
            
            ind_z_sub = find(zz_diff > slz/2) + segchain_ind{i}(jk);
            if isempty(ind_z_sub) == 0
                for j = 1:length(ind_z_sub)
                    coord_m(ind_z_sub(j) + 1:end,3) = coord_m(ind_z_sub(j) + 1:end,3) - slz;
                end
            end
            
            ind_z_add = find(zz_diff < -slz/2) + segchain_ind{i}(jk);
            if isempty(ind_z_add) == 0
                for j = 1:length(ind_z_add)
                    coord_m(ind_z_add(j) + 1:end,3) = coord_m(ind_z_add(j) + 1:end,3) + slz;
                end
            end
            
            xx_diff = diff(coord_array(segchain_ind{i}(jk+1):-1:1,1,i));
            yy_diff = diff(coord_array(segchain_ind{i}(jk+1):-1:1,2,i));
            zz_diff = diff(coord_array(segchain_ind{i}(jk+1):-1:1,3,i));
            
            ind_x_sub = segchain_ind{i}(jk+1) - find(xx_diff > slx/2);
            if isempty(ind_x_sub) == 0
                for j = 1:length(ind_x_sub)
                    coord_m(ind_x_sub(j):-1:1,1) = coord_m(ind_x_sub(j):-1:1,1) - slx;
                end
            end
            
            ind_x_add = segchain_ind{i}(jk+1) - find(xx_diff < -slx/2);
            if isempty(ind_x_add) == 0
                for j = 1:length(ind_x_add)
                    coord_m(ind_x_add(j):-1:1,1) = coord_m(ind_x_add(j):-1:1,1) + slx;
                end
            end
            
            ind_y_sub = segchain_ind{i}(jk+1) - find(yy_diff > sly/2);
            if isempty(ind_y_sub) == 0
                for j = 1:length(ind_y_sub)
                    coord_m(ind_y_sub(j):-1:1,2) = coord_m(ind_y_sub(j):-1:1,2) - sly;
                end
            end
            
            ind_y_add = segchain_ind{i}(jk+1) - find(yy_diff < -sly/2);
            if isempty(ind_y_add) == 0
                for j = 1:length(ind_y_add)
                    coord_m(ind_y_add(j):-1:1,2) = coord_m(ind_y_add(j):-1:1,2) + sly;
                end
            end
            
            ind_z_sub = segchain_ind{i}(jk+1) - find(zz_diff > slz/2);
            if isempty(ind_z_sub) == 0
                for j = 1:length(ind_z_sub)
                    coord_m(ind_z_sub(j):-1:1,3) = coord_m(ind_z_sub(j):-1:1,3) - slz;
                end
            end
            
            ind_z_add = segchain_ind{i}(jk+1) - find(zz_diff < -slz/2);
            if isempty(ind_z_add) == 0
                for j = 1:length(ind_z_add)
                    coord_m(ind_z_add(j):-1:1,3) = coord_m(ind_z_add(j):-1:1,3) + slz;
                end
            end
            
            if i ==1 && jk == 1
                coord_m_full = {coord_m};
            else
                coord_m_full = [coord_m_full;coord_m];
            end
            
            for kt = jk+1:length(temp_seg)
                
                if segchain_ind{i}(kt+1) ==n
                    temp_seg_array = coord_m(segchain_ind{i}(kt)+1:segchain_ind{i}(kt+1),:);
                else
                    temp_seg_array = coord_m(segchain_ind{i}(kt)+1:segchain_ind{i}(kt+1)+1,:);
                end
                
                diff_mat = abs(temp_seg_array - coord_array_seg{temp_seg(kt)});
                
                if isempty(find(diff_mat > 1, 1))
                    seg_same = [seg_same;kt];
                end
            end  
        end
        
        seg_same_full = [seg_same_full;sort(seg_same)];
        diff_chainseg{i,1} = diff_seg;
        
    end
    
    diffseg_chain = zeros(N,1);
    
    for m = 1:N
        diffseg_chain(m) = size(diff_chainseg{m},1);
    end
    
    diffcumseg_chain = cumsum(diffseg_chain);
 
 
 %  Assigns each of different images of the same chain a unique identity
 %  At this point of the time, we have N chains and their different images(translated in 
 %  x, y and z direction by simulation box dimension but have some part in the simulation box).
 %  In this configuration, we choose an arbitrary image  of a chain and assign [0 0 0],
 %  and other images of the same chain are assigned [d_x, d_y, d_z], where d_i's are  translation 
 % of image in the corresponding directions. 
 idnty = zeros(diffcumseg_chain(end),3);
 
 for i = 1:N
            idnty(diffcumseg_chain(i)-diffseg_chain(i)+1,:) = [0 0 0];
            for j = diffcumseg_chain(i)-diffseg_chain(i)+1:diffcumseg_chain(i)
                idnty(j,1) = idnty(diffcumseg_chain(i)-diffseg_chain(i)+1,1) + round((coord_m_full{j}(1,1) - coord_m_full{diffcumseg_chain(i)-diffseg_chain(i)+1}(1,1))/slx);
                idnty(j,2) = idnty(diffcumseg_chain(i)-diffseg_chain(i)+1,2) + round((coord_m_full{j}(1,2) - coord_m_full{diffcumseg_chain(i)-diffseg_chain(i)+1}(1,2))/sly);
                idnty(j,3) = idnty(diffcumseg_chain(i)-diffseg_chain(i)+1,3) + round((coord_m_full{j}(1,3) - coord_m_full{diffcumseg_chain(i)-diffseg_chain(i)+1}(1,3))/slz);
            end
        end

    

    cd(strcat(main_dir,'/../output_mat_file'));

    save(strcat('result_linknu0'));
    
    display(num);
end

cd(main_dir)
