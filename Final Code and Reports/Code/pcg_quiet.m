% helper function used to suppress print statements from pcg()
function [x,flag,relres,iter,resvec] = pcg_quiet(varargin)
[x,flag,relres,iter,resvec] = pcg(varargin{:});