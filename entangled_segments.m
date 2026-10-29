function [] = entangled_segments(r_cut, abs_lnk_lim,n_c)

load(strcat(pwd,'/../output_mat_file/result_linknu0.mat'));
seg_indices = [];

count = 1;

for i = 1:length(diff_seg_full)-1

    
    for j = i+1:length(diff_seg_full)
        
        seg_indices = [seg_indices; segmenting(i, j, coord_m_full{i}, coord_m_full{j}, r_cut, abs_lnk_lim,n_c)];

    end
         disp(num2str(count));
         count = count + 1; 
end
  main_dir =pwd;
  cd(strcat(main_dir,'/../output_mat_file'));
  save('seg_data.mat', 'seg_indices', 'coord_m_full');
  cd(main_dir)
end

