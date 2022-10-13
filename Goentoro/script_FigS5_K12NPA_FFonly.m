% script_FigS5_K12NPA_FFonly

%
% K12_PA
%
x0 = 1;
x1 = 10;
y0 = x0./(K11 + x0);
y1 = x1./(K11 + x1);
RHS = y1' - y0';
LHS = (1/x0-1/x1) - (1./K22)*(1./(K11' + x1) - 1./(K11' + x0));
K12_PA = repmat(RHS,nK,1)./LHS;

%
% creating K12_NPAplus and minus arrays
%
RHS = K11'.*(y1' - y0'/(1 + vareps));
LHS = 1/(1 + vareps)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAplus = repmat(RHS,nK,1)./LHS;

% K12_NPAminus
vareps1 = -vareps;
RHS = K11'.*(y1' - y0'/(1 + vareps1));
LHS = 1/(1 + vareps1)*(1 + K111/x0 + (1./K22/x0)*(K11'.*y0')) - ...
	(1 + K111/x1 + (1./K22/x1)*(K11'.*y1'));
K12_NPAminus = repmat(RHS,nK,1)./LHS;