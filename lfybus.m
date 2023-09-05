% Define the imaginary unit as 'j'
j = sqrt(-1);  
i = sqrt(-1);  

% Starting bus numbers for each line
nl = linedata(:, 1); 

 % Ending bus numbers for each line
nr = linedata(:, 2); 

% Resistance of each line
R = linedata(:, 3); 

% Reactance of each line
X = linedata(:, 4);   

% Line charging susceptance (imaginary part)
Bc = j * linedata(:, 5);  

% Line tap ratios
a = linedata(:, 6);   

% Number of branches (lines)
nbr = length(linedata(:, 1));  

% Total number of buses in the system
nbus = max(max(nl), max(nr));  

% Line impedance (real + imaginary part)
Z = R + j * X;  

% Branch admittance (inverse of impedance)
y = ones(nbr, 1) ./ Z; 

% Initialize Ybus to zero
Ybus = zeros(nbus, nbus);  

% Formation of the off-diagonal elements
for k = 1:nbr
    if a(k) <= 0
        a(k) = 1;  % If tap ratio is non-positive, set it to 1
    else
        % No action needed if tap ratio is positive
    end
    
    Ybus(nl(k), nr(k)) = Ybus(nl(k), nr(k)) - y(k) / a(k);
    Ybus(nr(k), nl(k)) = Ybus(nl(k), nr(k));  % Ybus is symmetric, so the element is copied
end

% Formation of the diagonal elements
for n = 1:nbus
    for k = 1:nbr
        if nl(k) == n
            Ybus(n, n) = Ybus(n, n) + y(k) / (a(k)^2) + Bc(k);
        elseif nr(k) == n
            Ybus(n, n) = Ybus(n, n) + y(k) + Bc(k);
        else
            % No action needed for other elements
        end
    end
end


% Separate real and imaginary parts of variables
R_real = real(R);
R_imaginary = imag(R);
X_real = real(X);
X_imaginary = imag(X);
Bc_real = real(Bc);
Bc_imaginary = imag(Bc);
Z_real = real(Z);
Z_imaginary = imag(Z);
y_real = real(y);
y_imaginary = imag(y);
Ybus_real = real(Ybus);
Ybus_imaginary = imag(Ybus);