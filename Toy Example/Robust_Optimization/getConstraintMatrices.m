

function [Ae,be,A,b] = getConstraintMatrices(constr)
% [Ae,be,A,b] = getConstraintMatrices(constr)
% Returns the problem's matrices
%
% i.e.    Ae x == be
%          A x <= b


% validity check
if isempty(constr)
    constr = [];
end


% problem variables
vars = recover(getvariables(constr));

% get the equality and inequality constaints
F_eq = extractConstraints(constr,'equality');
LFeq = length(F_eq);
if LFeq
    F_in = constr - F_eq;
else
    F_in = constr;
end
LFin = length(F_in);

% You have to split the equalities and inequalies because the sdpvar
% command works for con >= 0 and do not handle the equalities properly
% (please be carefull, the way it is used here is the proper).

if LFeq
    FeqVar = sdpvar(F_eq);
    [Ae,be,~] = getParamMatrices(FeqVar,vars);
else
    Ae = [];
    be = [];
end

if LFin
    FinVar = sdpvar(F_in);
    [A,b,~] = getParamMatrices(FinVar,vars);
else
    A = [];
    b = [];
end



% assign the correct signs, the sdpvar command gives the expression >= 0...
Ae = -Ae;
A = -A;
end







