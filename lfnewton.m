ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;
nbus = length(busdata(:,1));  % Get the number of buses in the system
for k=1:nbus
    n=busdata(k,1);
     % Assign values from busdata to variables
    kb(n)=busdata(k,2);
    Vm(n)=busdata(k,3);
    delta(n)=busdata(k, 4); 

    % Assign values from busdata to variables
    Pd(n)=busdata(k,5);
    Qd(n)=busdata(k,6); 
    Pg(n)=busdata(k,7);
    Qg(n) = busdata(k,8);  
    Qmin(n)=busdata(k, 9);
    
    % Assign values from busdata to variables
    Qmax(n)=busdata(k, 10);  

    % Assign values from busdata to variables
    Qsh(n)=busdata(k, 11);  
    
    if Vm(n) <= 0  
        Vm(n) = 1.0; 

    % If voltage magnitude is less than or equal to 0, set it to 1 and assign V(n) as 1+j*0
        V(n) = 1 + j*0;  
    else 
        % Convert delta from degrees to radians
        delta(n) = pi/180*delta(n);  

        % Calculate the complex voltage
        V(n) = Vm(n)*(cos(delta(n)) + j*sin(delta(n)));  

        % Calculate the active power injection
        P(n)=(Pg(n)-Pd(n));  

        % Calculate the reactive power injection
        Q(n)=(Qg(n)-Qd(n)+ Qsh(n)); 

        % Calculate the complex power injection
        S(n) = P(n) + j*Q(n);  
    end
end

for k=1:nbus
    if kb(k) == 1
        % Count the number of slack buses
        ns = ns+1;  
    else
    end
    if kb(k) == 2

    % Count the number of PV buses
        ng = ng+1;  
    else
    end

    % Store the count of PV buses for each bus
    ngs(k) = ng;  

    % Store the count of slack buses for each bus
    nss(k) = ns;  
end

% Get the magnitude of the bus admittance matrix
Ym=abs(Ybus);  

% Get the phase angle of the bus admittance matrix
t = angle(Ybus);  

% Calculate the number of unknowns
m=2*nbus-ng-2*ns;  
maxerror = 1; 
converge=1;
iter = 0;

% Start of iterations
clear A  DC   J  DX

% Test for max. power mismatch
while maxerror >= accuracy && iter <= maxiter 
    for i=1:m
        for k=1:m

           % Initializing Jacobian matrix
            A(i,k)=0;      
        end
    end
    
    % Incrementing the iteration counter
    iter = iter+1;
    
    
    for n=1:nbus
        % Adjusted bus index for matrices
        nn=n-nss(n);

        % Adjusted bus index for matrices
        lm=nbus+n-ngs(n)-nss(n)-ns;

        % Initialize element J11,J22,J33,J44 of Jacobian matrix
        J11=0; J22=0; J33=0; J44=0;
        
        for i=1:nbr
            if nl(i) == n || nr(i) == n
                if nl(i) == n
                    l = nr(i); % Neighboring bus index
                end
                if nr(i) == n
                    l = nl(i); % Neighboring bus index
                end
                
                 % Element J11 and J33 calculation
                J11=J11+ Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                J33=J33+ Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
                
                if kb(n)~=1
                    % Element J22 and J44 calculation
                    J22=J22+ Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
                    J44=J44+ Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                else
                end
                                
             
                if kb(n) ~= 1  && kb(l) ~=1
                    % Adjusted bus index for matrices
                    lk = nbus+l-ngs(l)-nss(l)-ns; 
                    ll = l -nss(l);
                    
                    % off-diagonal elements of J1
                    A(nn, ll) =-Vm(n)*Vm(l)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                    
                    if kb(l) == 0  % off-diagonal elements of J2
                        A(nn, lk) =Vm(n)*Ym(n,l)*cos(t(n,l)- delta(n) + delta(l));
                    end
                    
                    if kb(n) == 0  % off-diagonal elements of J3
                        A(lm, ll) =-Vm(n)*Vm(l)*Ym(n,l)*cos(t(n,l)- delta(n)+delta(l));
                    end
                    
                    if kb(n) == 0 && kb(l) == 0  % off-diagonal elements of J4
                        A(lm, lk) =-Vm(n)*Ym(n,l)*sin(t(n,l)- delta(n) + delta(l));
                    end
                else
                end
            else
            end
        end



        % Active power injection at bus n
        Pk = Vm(n)^2*Ym(n,n)*cos(t(n,n))+J33;

        % Reactive power injection at bus n
        Qk = -Vm(n)^2*Ym(n,n)*sin(t(n,n))-J11;
        
        if kb(n) == 1
          % Swing bus P and Q

            P(n)=Pk; % Set active power injection at Swing bus n
            Q(n) = Qk; % Set reactive power injection at Swing bus n 
        end   
        
        if kb(n) == 2
            Q(n)=Qk;% Set reactive power injection at PV bus n
            
            if Qmax(n) ~= 0

                % Calculate generator Mvar injection
                Qgc = Q(n)*basemva + Qd(n) - Qsh(n);

        % Check if the calculated generator Mvar is within limits
        % If not within limits, adjust the voltage magnitude to bring it within the specified limits
        % Between the 2nd and 6th iterations
                if iter <= 7
                    if iter > 2
                        if Qgc  < Qmin(n)
    
          % Adjust voltage magnitude to bring generator Mvar within limits
                            Vm(n) = Vm(n) + 0.01; % Increase voltage magnitude by 0.01 pu
                        elseif Qgc  > Qmax(n)  
                            Vm(n) = Vm(n) - 0.01;  % Decrease voltage magnitude by 0.01 pu 
                        end
                    else
                    end
                else
                end
            else
            end
        end
        
        if kb(n) ~= 1
            % diagonal elements of J1
            A(nn,nn) = J11; 
            DC(nn) = P(n)-Pk;
        end
        
        if kb(n) == 0
            % diagonal elements of J2,J3,J4
            A(nn,lm) = 2*Vm(n)*Ym(n,n)*cos(t(n,n))+J22;  
            A(lm,nn)= J33;   
            A(lm,lm) =-2*Vm(n)*Ym(n,n)*sin(t(n,n))-J44; 

            %updates the reactive power mismatch for the bus n in the DC vector.
            DC(lm) = Q(n)-Qk;
        end
    end
if iter == 1
                 
    % Extract the J11 matrix from the Jacobian matrix A
    J11 = A(1:nbus-ng-2*ns, 1:nbus-ng-2*ns);

    % Display the J11 matrix at iteration 1
    disp('J11 at iteration 1:');
    disp(J11);

    % Extract the J22 matrix from the Jacobian matrix A
    J22 = A(nbus-ng-2*ns+1:nbus-ng-ns, nbus-ng-2*ns+1:nbus-ng-ns);

    % Display the J22 matrix at iteration 1
    disp('J22 at iteration 1:');
    disp(J22);

    % Extract the J33 matrix from the Jacobian matrix A
    J33 = A(nbus-ng-ns+1:nbus-ns, nbus-ng-ns+1:nbus-ns);

    % Display the J33 matrix at iteration 1
    disp('J33 at iteration 1:');
    disp(J33);

    % Extract the J44 matrix from the Jacobian matrix A
    J44 = A(nbus-ns+1:nbus, nbus-ns+1:nbus);

    % Display the J44 matrix at iteration 1
    disp('J44 at iteration 1:');
    disp(J44);  
    
end


    
   % Solve for the change in state variables (voltages and angles)
    DX=A\DC';
    
    for n=1:nbus
        nn=n-nss(n); % Adjust bus index for state variables (angles)
        
        % Adjust bus index for state variables (voltages)
        lm=nbus+n-ngs(n)-nss(n)-ns;
        
        if kb(n) ~= 1
            delta(n) = delta(n)+DX(nn);  % Update voltage angles
        end
        
        if kb(n) == 0
            Vm(n)=Vm(n)+DX(lm);  % Update voltage magnitudes
        end
    end
    
    maxerror=max(abs(DC));  % Calculate maximum power mismatch
    
end

V = Vm.*cos(delta)+j*Vm.*sin(delta);  % Calculate bus voltages
deltad=180/pi*delta;  % Convert voltage angles to degrees

i=sqrt(-1);
k=0;

for n = 1:nbus
    if kb(n) == 1
        k=k+1;

        % Complex power injection at swing bus
        S(n)= P(n)+j*Q(n);

         % Real power generated at swing bus
        Pg(n) = P(n) + Pd(n);

        % Reactive power generated at swing bus
        Qg(n) = Q(n) + Qd(n) - Qsh(n);
        
        Pgg(k)=Pg(n);% Store real power generated for swing bus

        Qgg(k)=Qg(n);% Store reactive power generated for swing bus

    elseif  kb(n) ==2
        k=k+1;
        S(n)=P(n)+j*Q(n);% Complex power injection at PV bus

        Qg(n) = Q(n) + Qd(n) - Qsh(n);% Reactive power generated at PV bus

        Pgg(k)=Pg(n);% Store real power generated for PV bus

        Qgg(k)=Qg(n);% Store reactive power generated for PV bus
    end
    yload(n) = (Pd(n)- j*Qd(n)+j*Qsh(n))/(Vm(n)^2);  % Calculate the admittance of the load
end

busdata(:,3)=Vm';  % Update voltage magnitudes in the bus data

busdata(:,4)=deltad';  % Update voltage angles in degrees in the bus data

Pgt = sum(Pg);  % Calculate total real power generation
Qgt = sum(Qg);  % Calculate total reactive power generation
Pdt = sum(Pd);  % Calculate total real power demand
Qdt = sum(Qd);  % Calculate total reactive power demand
Qsht = sum(Qsh);  % Calculate total reactive power demand from shunt capacitors and reactors





