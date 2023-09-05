% BARIGYE OSBERT 2021171055 POWER ENGINEERING FINAL PROJEECT
%LINE FLOW ANALYSIS CODE 


SLT = 0; %Total line losses (complex power) in the system initialized to 0.

fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--  Power at bus & line flow    --Line loss--  Transformer\n')
fprintf('     from  to    MW      Mvar     MVA       MW      Mvar      tap\n')

% Loop over each bus in the system
for n = 1:nbus
    busprt = 0;
    for L = 1:nbr
        if busprt == 0
            fprintf('   \n')
            fprintf('%6g', n)  % Bus number
            fprintf('      %9.3f', P(n))  % Real power at the bus
            fprintf('%9.3f', Q(n))  % Reactive power at the bus 
            fprintf('%9.3f\n', abs(S(n)))  % Apparent power at the bus 
            busprt = 1;
        else
            % No action needed if the bus has already been printed
        end
        % Check if the line is connected to the current bus (n)
       if nl(L) == n
            k = nr(L);
            % Current flowing into the line from the bus
            In = (V(n) - a(L)*V(k))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(n);  

            % Current flowing into the line from the other bus
            Ik = (V(k) - V(n)/a(L))*y(L) + Bc(L)*V(k);  

            % Complex power injection at the "from" bus
            Snk = V(n)*conj(In);  

            % Complex power injection at the "to" bus
            Skn = V(k)*conj(Ik);  
            
            % Total complex power loss in the line
            SL  = Snk + Skn;  

            % Accumulate the total line losses
            SLT = SLT + SL;  
        elseif nr(L) == n
            k = nl(L);

            % Current flowing into the line from the bus
            In = (V(n) - V(k)/a(L))*y(L) + Bc(L)*V(n);  

            % Current flowing into the line from the other bus
            Ik = (V(k) - a(L)*V(n))*y(L)/a(L)^2 + Bc(L)/a(L)^2*V(k); 

            % Complex power injection at the "from" bus
            Snk = V(n)*conj(In);  

            % Complex power injection at the "to" bus
            Skn = V(k)*conj(Ik); 

            % Total complex power loss in the line
            SL  = Snk + Skn;  

            % Accumulate the total line losses
            SLT = SLT + SL;  
        else
            % No action needed if the line is not connected to the current bus
       end

       if nl(L) == n || nr(L) == n

            % Bus number connected to the line
            fprintf('%12g', k)  
            
           % Real power injection at the "from" bus in PU
            fprintf('%9.3f', real(Snk))  

            % Reactive power injection at the "from" bus in PU
            fprintf('%9.3f', imag(Snk))  
            
            % Apparent power injection at the "from" bus in PU
            fprintf('%9.3f', abs(Snk))

            % Real power loss in the line in PU
            fprintf('%9.3f', real(SL))  
            if nl(L) == n && a(L) ~= 1

                % Reactive power loss in the line in Mvar (only for "from" bus)
                fprintf('%9.3f', imag(SL))  

                 % Transformer tap ratio
                fprintf('%9.3f\n', a(L))  
            else

                 % Reactive power loss in the line in Mvar (for "to" bus or tap ratio = 1)
                fprintf('%9.3f\n', imag(SL))  
            end
        else
            % No action needed if the line is not connected to the current bus
        end
    end
end

 % Total line losses divided by 2 (since each line contributes to losses twice)
SLT = SLT/2; 
fprintf('   \n')
fprintf('    Total loss                         ')
fprintf('%9.3f', real(SLT))  % Total real power loss in the system in pu
fprintf('%9.3f\n', imag(SLT))  % Total reactive power loss in the system in pu

% Clearing temporary variables
clear Ik In SL SLT Skn Snk