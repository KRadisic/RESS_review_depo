function [a, c] = union_multiindices(a_old, c_old, a_new, c_new)
if (size(a_old,1) ~= length(c_old) | size(a_new,1) ~= length(c_new))
    msg = 'The number of coefficients must be compatible with the basis dimension.';
    error(msg)
elseif (size(a_old,1) == 0)
    a = a_new;
    c = c_new;
elseif (size(a_new,1) == 0)
    a = a_old;
    c = c_old;
else
a = a_old;
c = c_old;
i = 1;
j = 1;
    while  (j <= size(a_new,1)) % & i <= size(a_old,1)
            % if the fixed new line already exists in the old
            if ((a_old(i,1) == a_new(j,1)) & ...
                (a_old(i,2) == a_new(j,2))) 
                % just change coeff, and go to the next new line
                %disp('1st case');
                c(i) = max(abs(c_old(i)),abs(c_new(j)));
                j = j+1;
                i = 1;
            % if we went through all the old lines and the new line does not exist in the old 
            elseif (i == size(a_old,1))
                % add in the new line, with the good coeff
                a(length(a)+1,:) = a_new(j,:);
                c(length(c)+1) = c_new(j);
                % go to the next new line
                %disp('2nd case');
                j = j+1;
                i = 1;
            % while we're not reaching the end of old line
            % if the new line not equal to the current old one, go to next old line
            else
                %disp('3rd case');                
                i = i+1;
            end
        % we either got till the end of old lines and it wasn't equal to
        % anyone, or we foung an old line along the way to which it was
        % equal to. Thus, no other incrementation is needed for index j.
    end
end
end
