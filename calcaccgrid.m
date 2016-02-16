function A = calcaccgrid(fd, area)
   
    %% Compute accumulation area for each grid cell by solving sparse system of
    %% equations in cell areas

    % preallocate sparse matrix arrays
    [m,n] = size(fd.grid);
    len = m*n;
    ii = zeros(len, 1);
    jj = zeros(len, 1);
    val = zeros(len, 1);
    k = 0;

    % populate coefficient matrix
    for(i=1:m)
        for(j=1:n)

            [i_out, j_out, flag] = getfdneighb(i, j, [m,n], fd.grid(i, j));

            source = sub2ind([m, n], i, j);

            % resize as necessary
            if(k > len)
                len = 2*len;
                ii(len) = 0;
                jj(len) = 0;
                val(len) = 0;
            end
            
            % populate diagonal
            k = k+1;
            ii(k) = source;
            jj(k) = source;
            val(k) = 1;
            
            % off diagonal entries
            if(flag)
                k = k+1;
                dest = sub2ind([m, n], i_out, j_out);
                ii(k) = dest;
                jj(k) = source;
                val(k) = -1;
            end
            
        end
    end
   
    % solve system of equations
    M = sparse(ii(1:k), jj(1:k), val(1:k));
    a = reshape(area, [m*n, 1]); % vec out area grid for right-hand side
    x = M\a;
    A = reshape(x, [m, n]);

end
