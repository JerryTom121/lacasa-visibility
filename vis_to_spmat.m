function K = vis_to_spmat(T)

T = T + 1;
T = T(sum(T,2) ~= 0,:);
K = sparse(T(:,1),T(:,2),1);