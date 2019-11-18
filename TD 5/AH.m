% take the volatily matrix, the step j, dT,dK, the strike list K, Cobs
% matrix, x the volatility
function output = AH(tf_matrix,j,x,C,dT,dK,K)
%complete the matrix between the observation we have
    m = round((tf_matrix(1:end-1,1) + tf_matrix(2:end,1)) / 2); % list of 1 or 0
    NK = length(K);
    %ksigma vector that we need to complete between the observation 
    kvol = NaN(NK,1);
    k = 1;
    for p = 1:length(m)
        while k < m(p)
            kvol(k) = x(p)
            k = k+1;
        end
        kvol(k:end) = x(end);
    end

    %We construct the matrix A
    % using formula form slide 16 lecture 4 with r=q=0
    % Ii ci ri have a common part 0.5dt/dK^2 kvol^2
    % Ii = -0.5dt/dK^2 sigma(i,j)2k(i,j)2
    % ci = 1+dt/dK^2sigma(i,j)2k(i,j)2 
    % ri = -0.5dt/dK^2 sigma(i,j)2k(i,j)2
    %A(0,0) = x0 = 1
    z = 0.5 * dT/dK^2 * kvol(2:end).^2;
    z(end) = [];  
    % returns an array containing a copy of z arange in a matrix 1 line and 3 col
    diag = repmat(z,1,3).*[-1 2 -1]; % to have the correct Ii Ci ri 
    % we create A line by line
    % 1st line 
    A = zeros(1,NK);
    % for line 2 to line NK-1 we create a line only with 
    % 3 elements in the dia Ii,Ci,ri (that we have using the vector z
    for m = 1:NK-2
        % A =  O diag 0 
        A = [A; zeros(1,m-1) diag(m,:) zeros(1,NK-2-m)];
    end
    % add last line
    A = [A; zeros(1,NK)];
    %using the fact that
    % ci = 1+dt/dK^2sigma(i,j)2k(i,j)2 and that we fill only
    % dt/dK^2sigma(i,j)2k(i,j)2 in the diag we add Identity matrix
    A = A + eye(NK);
    %Cmodel(tj+1) = A(tj+1)^-1 Cmodel(tj)
    output = A\C(:,j);
end