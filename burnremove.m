% Burnaby size removal
% Burnaby's technique allows to remove the effect of isometric size from
% multivariate morphometric data by projecting the data points onto a plane 
% orthogonal to a size-related axis (e.g., isometric size vector or PC1).
%
% Syntax: shape=burnremove(data)
%         shape=burnremove(data,vsize)
%
% Input:
%     DATA  is an N-by-M matrix of morphometric data with each row 
%           representing a multivariate observation, and each column
%           representing a variable.
%     VSIZE is a M-component vector representing the size axis.
%
% Output:
%     SHAPE is N-by-M matrix of size-adjusted data.
%
% The formulations are from McCOY et al 2006. Default coefficients for 
% isometric size vector are from JOLICOEUR 1963 and JOLICOEUR 1984.

% Igor Dzeverin 07-May-2010 

function shape=burnremove(data,vsize)

d=size(data);
if nargin<2 vsize=ones(d(2),1)/sqrt(d(2)); end
    
s=size(vsize);
if min(s)>1, error('VSIZE is not a vector.'); end
if s(1)<s(2) vsize=vsize'; end
if s(1)~=d(2), error('The dimensions of DATA and VSIZE must agree.'); end

shape=data*(eye(d(2))-vsize/(vsize'*vsize)*vsize');