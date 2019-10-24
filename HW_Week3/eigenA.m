function [minEigen, maxEigen] = eigenA(n)
    A = Amatrix(n);
    [V,D] = eig(A);
    d = diag(D);
    minEigen = min(d);
    maxEigen = max(d);
end

function A = Amatrix(n)
    % A is a (2,2)-tensor(rank 4), with Following properties:
    % A_{i,j}^{i+1,j  } = 1/h^2;
    % A_{i,j}^{i-1,j  } = 1/h^2;
    % A_{i,j}^{i  ,j+1} = 1/h^2;
    % A_{i,j}^{i  ,j-1} = 1/h^2;
    % A_{i,j}^{i  ,j  } = -4/h^2;
    % With the formalizm of:
    % \Delta u = u_xx + u_yy = f, Numerically derivated as (mersh aera to n*n)
    % \Delta u = \frac{ u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j} }{h^2}
    % Corresponding to the vectorized order:
    % u'=[u_11, u_21, u_31,...,u_n1, u_12, u_22, ..., u_nn]
    % f'=[f_11, f_21, f_31,...,f_n1, f_12, f_22, ..., f_nn]
    % Au = f
    % The A(n^2 by n^2) Matrix column and rows correspond to tensor script (n by n by n by n) by:
    % A_{ii,jj}^{kk,ll} = A( kk+ll*(n-1), ii+jj*(n-1) )
    try
        assert(n>0 && round(n)==n )
    catch
        error('MyComponent:incorrectInput',...
            'Error.\n n should be a positive interger instead of %f', n);
    end
    A = zeros(n*n,n*n);
    for ii = 1:n
        for jj = 1:n
            if    ii+1 + n* (jj   -1)<n*n+1 && ii + n* (jj -1) < n*n+1 && ii+1 + n* (jj   -1)>0
                A(ii+1 + n* (jj   -1),         ii + n* (jj -1) ) = 1;
            end
            if    ii-1 + n* (jj   -1)<n*n+1 && ii + n* (jj -1) < n*n+1 && ii-1 + n* (jj   -1)>0
                A(ii-1 + n* (jj   -1),         ii + n* (jj -1) ) = 1;
            end
            if    ii   + n* (jj+1 -1)<n*n+1 && ii + n* (jj -1) < n*n+1 && ii   + n* (jj+1 -1)>0
                A(ii   + n* (jj+1 -1),         ii + n* (jj -1) ) = 1;
            end
            if    ii   + n* (jj-1 -1)<n*n+1 && ii + n* (jj -1) < n*n+1 && ii   + n* (jj-1 -1)>0
                A(ii   + n* (jj-1 -1),         ii + n* (jj -1) ) = 1;
            end
            A(ii   + n* (jj   -1), ii + n* (jj -1) ) = -4;
        end
    end
%     A = A/h;
end