function shim = cvx_cbc( shim )

% CVX_SOLVER_SHIM	CBC interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

if ~isempty( shim.solve )
    return
end
if isempty( shim.name )
    fname = 'opti_cbc.m';
    ps = pathsep;
    shim.name = 'CBC';
    shim.dualize = true;
    flen = length(fname);
    fpaths = which( fname, '-all' );
    if ~iscell(fpaths)
      fpaths = { fpaths };
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths)
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) )
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        tshim.location = new_dir;
        if isempty( tshim.error )
            tshim.check = @check;
            tshim.solve = @solve;
            tshim.eargs = {};
            if k ~= 1
                tshim.path = [ new_dir, ps ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim )
        shim = oshim;
        shim.error = 'Could not find a CBC installation.';
    end
else
    shim.check = @check;
    shim.solve = @solve;
end
    
function found_bad = check( nonls ) %#ok
found_bad = false;

function [ x, status, tol, iters, y, z] = solve( At, b, c, nonls, quiet, prec, settings )

n  = length( c );
m  = length( b );
lb = -Inf(n,1);
ub = +Inf(n,1);
vtype = 'C';
vtype = vtype(ones(n,1));
rr = zeros(0,1);
cc = rr; 
vv = rr;
zinv = rr;
is_ip = false;
for k = 1 : length( nonls )
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    tt = nonls( k ).type;
    if strncmp( tt, 'i_', 2 )
      is_ip = true;
      vtype(temp) = 'I';
      if strcmp(tt,'i_binary')
        vtype(temp) = 'B';
      end
    elseif nn == 1 || isequal( tt, 'nonnegative' )
        lb(temp) = 0;
    elseif isequal( tt, 'lorentz' )
        if nn == 2
            rr2  = [ temp ; temp ];
            cc2  = reshape( floor( 1 : 0.5 : 2 * nv + 0.5 ), 4, nv );
            vv2  = [1;1;-1;1]; vv = vv(:,ones(1,nv));
            rr   = [ rr ; rr(:) ];
            cc   = [ cc ; cc(:) ];
            vv   = [ vv ; vv(:) ];
            zinv = [ zinv ; temp(:) ];
        else
            error('CBC does not support nonlinear constraints.' );
        end
    else
      error('CBC does not support nonlinear constraints.' );
    end
end
if ~isempty(rr)
  znorm = [1:n]';
  znorm(zinv) = [];
  rr = [ rr ; znorm ];
  cc = [ cc ; znorm ];
  vv = [ vv ; ones(size(znorm)) ];
  reord = sparse( rr, cc, vv, n, n );
  At = reord' * At;
  c  = reord' * c;
end
H = [];
ru = b;
rl = b;
[xx, fmin, errnum, extra] = cvx_run_solver( @opti_cbc, H, full(c), At', full(rl), full(ru), lb, ub, vtype, 'xx', 'fmin', 'errnum', 'extra', settings, 8);
iters = extra.Nodes;
x = full( xx );
y = zeros(size(At,2), 1);
z = zeros(length(c), 1);
if ~isempty( rr )
  x = reord * x;
end
status = extra.Status;
if (errnum ~= 1)
  tol = Inf;
elseif strncmp(status,'Gap Reached', 3)
  tol = prec(3);
else
  tol = prec(2);
end

% Copyright 2005-2016 CVX Research, Inc.
% See the file LICENSE.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
