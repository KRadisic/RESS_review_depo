function [alpha, c_alpha] = get_multiindex(myPCE_LARS,n_inputs)
% Transform in a good form for the multiindices
for idx = 1:n_inputs
    poly_indices(idx,:) = find(myPCE_LARS.PCE.Basis.Indices(:,idx)>0);
end
%poly_indices3 = find(myPCE_LARS.PCE.Basis.Indices(:,3)>0);

alpha = zeros(size(myPCE_LARS.PCE.Basis.Indices,1), ...
    size(myPCE_LARS.PCE.Basis.Indices,2));
for idx = 1:n_inputs
    alpha(poly_indices(idx,:),idx) = myPCE_LARS.PCE.Basis.Indices(poly_indices(idx,:),idx);
end
%alpha(poly_indices3,3) = myPCE_LARS.PCE.Basis.Indices(poly_indices3,3);

c_alpha = zeros(size(myPCE_LARS.PCE.Basis.Indices,1),1);
for idx = 1:n_inputs
    c_alpha(poly_indices(idx,:)) = myPCE_LARS.PCE.Coefficients(poly_indices(idx,:));
end
% c_alpha(poly_indices3) = myPCE_LARS.PCE.Coefficients(poly_indices3);
end