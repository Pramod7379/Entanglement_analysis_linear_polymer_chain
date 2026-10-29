function org_chain_ent_points = ent_clustering()

load(strcat(pwd,'/../output_mat_file/result_linknu0.mat'));
Lx = [l(1) l(2)];
Ly = [l(3) l(4)];
Lz = [l(5) l(6)];
iden = idnty;

load(strcat(pwd,'/../output_mat_file/seg_data.mat'));
ent_points = cell(length(coord_m_full), 1);
for i = 1:length(coord_m_full)
    % Find entanglement points on a chain
    A = find(seg_indices(:, 1) == i);
    ent_A = seg_indices(A, 3);
    B = find(seg_indices(:, 2) == i);
    ent_B = seg_indices(B, 4);
    ent_chain = sort([ent_A; ent_B]);
    % Entanglement points within the box
    for j = 1:length(ent_chain)
        if (coord_m_full{i}(ent_chain(j), 1) <= Lx(2)) && ((coord_m_full{i}(ent_chain(j), 1) >= Lx(1)))
            if (coord_m_full{i}(ent_chain(j), 2) <= Ly(2)) && ((coord_m_full{i}(ent_chain(j), 2) >= Ly(1)))
                if (coord_m_full{i}(ent_chain(j), 3) <= Lz(2)) && ((coord_m_full{i}(ent_chain(j), 3) >= Lz(1)))
                    ent_points{i} = [ent_points{i}; ent_chain(j)];
                end
            end
        end
    end
end

% Original chain nos.
org_chains = [];
for i = 1:length(iden)
    if iden(i, :) == [0 0 0]
        org_chains = [org_chains; i];
    end
end
org_chains = [org_chains; (length(iden) + 1)];

org_chain_ent_points = cell((length(org_chains) - 1), 1);
for i = 1:(length(org_chains) - 1)
    ent_dum = [];
    for j = org_chains(i):(org_chains(i+1) - 1)
        ent_dum = [ent_dum; ent_points{j}];
    end
    org_chain_ent_points{i} = sort(ent_dum);
end

end

