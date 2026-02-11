function [a_reduced,c_reduced] = discard_small_regressors(a,c,eps,n_inputs)
    a_reduced = a;
    c_reduced = c;
    a_reduced(c.^2 < eps,:) = [];
    c_reduced(c.^2 < eps) = [];
    % add the constant term manually, if it is needed
    % if the function removes it
    if sum(a_reduced(1,:)) ~= 0
        a_reduced = [zeros(1,n_inputs);a_reduced];
        c_reduced = [0;c_reduced];
    end
end
