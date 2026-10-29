 function []= Topological_entanglement()
 
 % In this script topological entanglement length is calculated 
 
 load(strcat(pwd,'/../output_mat_file/result_linknu0.mat'));
 ent_points = ent_clustering();

 ent_points(cellfun('length',ent_points)<1) =[];
 num_atoms_unclustered = [];
 for i = 1:length(ent_points)

  num_atoms_unclustered = [num_atoms_unclustered; (n/length(ent_points{i}))];

 end
% Final average number of atoms between two entanglements "N_t"

num_atoms_unc = mean(num_atoms_unclustered);
disp(['The average unclustered entanglement length is ' num2str(num_atoms_unc) 'atoms between two entanglements.']);

  main_dir =pwd;
  cd(strcat(main_dir,'/../output_mat_file'));
  save('topo.mat', 'ent_points', 'num_atoms_unc');
  cd(main_dir)


 end

