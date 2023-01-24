function [P_mat, R_vec, w, i] = m3d_randpack(P_mat, R_vec, R_gen, i, n_total, box_length, o_f, counter, maxtrials, w)

%% Choose positioning algorithm:
if (max(R_gen(:,2))-min(R_gen(:,2)) == 0)
    
    %% Generate random particle positions for mono-sized particles:
    while (i <= n_total && counter <= maxtrials)
        R_vec(i,:) = R_gen(i,:);
        P_mat(i,:) = [i, box_length+box_length.*rand(1,3)];
        
        %% Compute inter-particle distances:
        N_mat = P_mat;
        dxdydz = bsxfun(@minus, N_mat(:,2:4), P_mat(i,2:4));
        dist_trans = sqrt(sum(dxdydz.*dxdydz, 2));
        dist_nearest = min(dist_trans(dist_trans > 0));
        ind_nearest = find(dist_trans == dist_nearest);
        ind_nearest_dxdydz = ind_nearest(1);
        ind_nearest = N_mat(ind_nearest_dxdydz,1);
        R_trans = o_f*(R_vec(i,2)+R_vec(ind_nearest,2));
        dist_eff_nearest = dist_nearest-R_trans;
        
        %% Check position of particle i:
        if (dist_eff_nearest < 0)
            ++counter;
            continue;
        elseif (dist_eff_nearest >= 0)
            exeyez = dxdydz(ind_nearest_dxdydz,:);
            exeyez = bsxfun(@rdivide, exeyez, dist_nearest);
            
            %% Translate particle i:
            P_mat(i,2:4) = P_mat(i,2:4)+dist_eff_nearest.*exeyez;
            
            %% Update particle positions:
            [P_mat, R_vec] = m3d_particlemap(P_mat, R_vec, i, n_total, box_length);
            waitbar((single(i)/single(n_total)), w, "Generating particle positions...");
            i = i+1;
        end
    end
    
elseif (max(R_gen(:,2))-min(R_gen(:,2)) > 0)
    
    %% Generate random particle positions for distributed particle diameters:
    while (i <= n_total && counter <= maxtrials)
        R_vec(i,:) = R_gen(i,:);
        P_mat(i,:) = [i, box_length+box_length.*rand(1,3)];
        
        %% Compute inter-particle distances:
        N_mat = P_mat;
        dxdydz = bsxfun(@minus, N_mat(:,2:4), P_mat(i,2:4));
        dist_trans = sqrt(sum(dxdydz.*dxdydz, 2));
        R_trans = o_f.*bsxfun(@plus, R_vec(i,2), R_vec(:,2));
        dist_eff = dist_trans-R_trans;
        dist_eff_nearest = min(dist_eff(dist_eff > 0));
        ind_nearest = find(dist_eff == dist_eff_nearest);
        ind_nearest_dxdydz = ind_nearest(1);
        
        %% Check position of particle i:
        if (sum(dist_eff < 0) > 1)
            ++counter;
            continue;
        elseif (sum(dist_eff < 0) <= 1)
            exeyez = dxdydz(ind_nearest_dxdydz,:);
            exeyez = bsxfun(@rdivide, exeyez, dist_trans(ind_nearest_dxdydz));
            
            %% Translate particle i:
            P_mat(i,2:4) = P_mat(i,2:4)+dist_eff_nearest.*exeyez;
            
            %% Update particle positions:
            [P_mat, R_vec] = m3d_particlemap(P_mat, R_vec, i, n_total, box_length);
            waitbar((single(i)/single(n_total)), w, "Generating particle positions...");
            i = i+1;
        end
    end
    
end

end
