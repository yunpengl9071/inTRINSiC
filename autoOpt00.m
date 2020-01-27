function [x0,fval,exitflag,output,lambda,grad,hessian] = autoOpt00(x0,t,m,gene,ub,l,l1)
%% This is an auto generated MATLAB file from Optimization Tool.
st = size(t); sm = size(m);
%% Start with the default options
options = optimoptions('fmincon');
%% Modify options setting
options = optimoptions('fmincon','GradObj','off');
options = optimoptions(options,'Display', 'off');
options = optimoptions(options,'Algorithm', 'interior-point');
options = optimoptions(options,'MaxFunEvals', 1000000);
options = optimoptions(options,'MaxIter', 1000);
options = optimoptions(options,'TolFun', 1e-6);
%options = optimoptions(options,'PlotFcns', { @optimplotfval });

exitflag = 100;


while exitflag ~= 1 && exitflag ~= 2
    [x0,fval,exitflag,output,lambda,grad,hessian] = ...
    fmincon(@(x)lsqL1L2(x,t,m,gene,l,l1),x0,[],[],[],[],zeros(st(2)+sm(2)+1,1),ub,[],options);
end
