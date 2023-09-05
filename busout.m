
% BARIGYE OSBERT
% PRINTING THE POWER FLOW ANALYSIS RESULTS IN TABLE FORM 
%clc
disp("POWER FLOW ANALYSIS USING NEWTON RAPHSON ")
fprintf('                      Maximum Power Mismatch = %g \n', maxerror)
fprintf('                             No. of Iterations = %g \n\n', iter)
header =['    Bus  Voltage  Angle    ------Load------    ---Generation---   Injected'
       '    No.  Mag.     Degree     PU       PUar       PU       PUar       PUar '
       '                                                                          '];


disp(header)
for n = 1:nbus
    fprintf(' %5g', n), fprintf(' %7.3f', Vm(n)),
    fprintf(' %8.3f', deltad(n)), fprintf(' %9.3f', Pd(n)),
    fprintf(' %9.3f', Qd(n)), fprintf(' %9.3f', Pg(n)),
    fprintf(' %9.3f ', Qg(n)), fprintf(' %8.3f\n', Qsh(n))
end
fprintf('      \n'), fprintf('    Total              ')
fprintf(' %9.3f', Pdt), fprintf(' %9.3f', Qdt),
fprintf(' %9.3f', Pgt), fprintf(' %9.3f', Qgt), fprintf(' %9.3f\n\n', Qsht)

