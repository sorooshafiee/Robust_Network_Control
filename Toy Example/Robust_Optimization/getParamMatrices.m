function [F,G,details] = getParamMatrices(P,xReal)
% Purpose: given an expression P and the parameters x that it depends on
%          return the matrices P = Fx + G (parametric on all the other
%          variables that may be present in your problem)
%          For fast implementation provide depRows (rows that depend on
%          xReal)
% Author:  Darivianakis Georgios
%
% email:   gdarivia@control.ee.ethz.ch

if nargin < 2
    xReal = recover(getvariables(P));
end

realVar = getvariables(xReal);
rVlen = length(realVar);
% Make sure they are sorted
x = recover(realVar);
Qtemp = [];
htemp = [];
for ii = 1:length(P)
    if sum(ismember(depends(P(ii)),realVar)) == 0
        Qtemp = [Qtemp;zeros(1,rVlen)];
        htemp = [htemp;P(ii)];
    else
        
        % much faster code...
        [Q,c,f,allvar,~] = quaddecompGD(P(ii),x);
        if ~isreal(Q) % Numerical noise common on imaginary parts
            Qr = real(Q);
            Qi = imag(Q);
            Qr(abs(Qr)<1e-10) = 0;
            Qi(abs(Qi)<1e-10) = 0;
            cr = real(c);
            ci = imag(c);
            cr(abs(cr)<1e-10) = 0;
            ci(abs(ci)<1e-10) = 0;
            Q = Qr + sqrt(-1)*Qi;
            c = cr + sqrt(-1)*ci;
        end
        
        
        used_variables = getvariables(allvar);
        
        notparameters = find(ismember(used_variables,getvariables(x)));%parameters);
        parameters = setdiff(1:length(used_variables),notparameters);
        if ~isempty(parameters)
            y = recover(used_variables(parameters));
            
            % parametric case (nonlinear Pparam)
            ht = f + c(parameters)'*y + y'*Q(parameters,parameters)*y;
            c(parameters) = [];
            Q2 = Q(notparameters,parameters);
            
            c = c + 2*Q2*y;
        else
            ht = f;
        end
        
        %% create the matrices
        Qtemp = [Qtemp;c'];
        htemp = [htemp;ht];
    end
end

F = Qtemp;
G = htemp;
details = x;

end