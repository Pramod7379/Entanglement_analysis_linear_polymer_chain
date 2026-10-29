function [] = Rheological_entanglement()

% In this script 1D clustering algorithm is performed to cluster the entanglement points on each chain

load(strcat(pwd,'/../output_mat_file/result_linknu0.mat'));
ent_points = ent_clustering();
ent_points(cellfun('length',ent_points)<3) = [];

num_atoms_all_chains = [];
clusters_num = [];
ke_mat=[];
for i = 1:length(ent_points)
    silh_coeff_score={};
    % k-means clustering for various k's
    sum_point_cluster_dist = [];
    sil_stand_dev = [];
    for k = 1:length(ent_points{i})

        [IDX, C, SUMD] = kmeans(ent_points{i}, k, 'MaxIter', 1000, 'Distance', 'cityblock');
   
        s = silhouette(ent_points{i}, IDX, 'cityblock');
        
        % silhouette_plot
        %       [s,h] = silhouette(ent_points{i}, IDX, 'cityblock');
        %       figure;
        sum_point_cluster_dist = [sum_point_cluster_dist; sum(SUMD)];
        sil_stand_dev = [sil_stand_dev; sqrt(var(s))];
        silh_coeff_score(k,1)={s};
        
        if k==1
            AA={IDX};
        else
            AA=[AA;{IDX}];
        end


    end
    
    %  elbow_point_ke_detected using knee algorithm
    ke=knee_pt(sum_point_cluster_dist);
    ke_mat=[ke_mat;ke];



    %  average silhouette score for each value of k's

    for N=1:length(silh_coeff_score)

        AAAA = silh_coeff_score{N,1};
        BBBB = AA{N,1};

        % Find unique indices in B
        uniqueIndices = unique(BBBB);

        % Initialize result array
        result = zeros(1, length(uniqueIndices));

        % Iterate over unique indices
        for I = 1:length(uniqueIndices)
            idx = uniqueIndices(I);

            % Calculate sum of elements of A corresponding to current index in B
            result(I) = sum(AAAA(BBBB == idx)) / sum(BBBB == idx);
        end
        summ{N,1}=result;
        avg_silh_score(N,1)=mean(summ{N,1});
    end


    % optimum number of clusters kc

    max_avg_silh_score= max(avg_silh_score(ke:ke+ceil(length(avg_silh_score)/9)));
    max_indices=find(avg_silh_score==max_avg_silh_score);
    kc=max_indices;

    % Final number of clusters is (Ki=kc-number of clusters containing only
    % one topological entanglement)

    new_point=AA{kc,1};
    unique_values = unique(new_point)';
    counts = histcounts(new_point, [unique_values, max(unique_values)+1]);
    Ki=sum(counts > 1);
    if Ki >1
        Ki_mod=Ki;
    else
        Ki_mod=1;
    end

    clusters_num = [clusters_num; Ki_mod];
    num_atoms_all_chains = [num_atoms_all_chains; (n/Ki_mod)];


end
% The average clustered entanglement length "Ne"
num_atoms = mean(num_atoms_all_chains);

disp(['The average clustered entanglement length is ' num2str(num_atoms) 'atoms between two entanglements.']);

  main_dir =pwd;
  cd(strcat(main_dir,'/../output_mat_file'));
  save('rheo.mat', 'clusters_num', 'num_atoms');
  cd(main_dir)



