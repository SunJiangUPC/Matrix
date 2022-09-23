function [x,flag,relres,iter,resvec] = bicgstabCSR(N,AI,AJ,AA,b,tol,maxit,...
    LI,LJ,LA,UI,UJ,UA,x0,varargin)
%BICGSTAB   BiConjugate Gradients Stabilized Method.

% Check for an acceptable number of input arguments
if nargin < 5
    error(message('MATLAB:bicgstab:NotEnoughInputs'));
end

% Assign default values to unspecified parameters
if nargin < 6 || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:bicgstab:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if nargin < 7 || isempty(maxit)
    maxit = min(N,20);
end
maxit = max(maxit, 0);

% Check for all zero right hand side vector => all zero solution
n2b = norm(b);                     % Norm of rhs vector, b
if (n2b == 0)                      % if    rhs vector is all zeros
    x = zeros(N,1);                % then  solution is all zeros
    flag = 0;                      % a valid solution has been obtained
    relres = 0;                    % the relative residual is actually 0/0
    iter = 0;                      % no iterations need be performed
    resvec = 0;                    % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,NaN);
    end
    return
end

if nargin >= 10
    existL = 1;
    Ltype = 'matrix';
else
    existL = 0;
    Ltype = 'matrix';
end

if nargin >= 13
    existU = 1;
    Utype = 'matrix';
else
    existU = 0;
    Utype = 'matrix';
end

if ((nargin >= 14) && ~isempty(x0))
    if ~isequal(size(x0),[N,1])
        error(message('MATLAB:bicgstab:WrongInitGuessSize', N));
    else
        x = x0;
    end
else
    x = zeros(N,1);
end

if (nargin >= 15)
    error(message('MATLAB:bicgstab:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                          % Iterate which has minimal residual so far
imin = 0;                          % Iteration at which xmin was computed
tolb = tol * n2b;                  % Relative tolerance
% r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
r = b - dgemvCSR(N,AI,AJ,AA,x);
% normr = norm(r);                   % Norm of residual
normr = dnorm(N,r);
normr_act = normr;

if (normr <= tolb)                 % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = 0;
    resvec = normr;
    if (nargout < 2)
        itermsg('bicgstab',tol,maxit,0,flag,iter,relres);
    end
    return
end

rt = r;                            % Shadow residual
resvec = zeros(2*maxit+1,1);       % Preallocate vector for norm of residuals
resvec(1) = normr;                 % resvec(1) = norm(b-A*x0)
normrmin = normr;                  % Norm of residual from xmin
rho = 1;
omega = 1;
stag = 0;                          % stagnation of the method
alpha = [];                        % overshadow any functions named alpha
moresteps = 0;
maxmsteps = min([floor(N/50),10,N-maxit]);
maxstagsteps = 3;
% loop over maxit iterations (unless convergence or failure)

for ii = 1 : maxit
    rho1 = rho;
%     rho = rt' * r;
    rho = ddot(N,rt,r);
    if (rho == 0.0) || isinf(rho)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    if ii == 1
        p = r;
    else
        beta = (rho/rho1)*(alpha/omega);
        if (beta == 0) || ~isfinite(beta)
            flag = 4;
            break
        end
%         p = r + beta * (p - omega * v);
        p = daxpby(N, r, daxpby(N,p,v,beta,-beta*omega),1,1);
        
    end
    if existL
%         ph = iterapp('mldivide',m1fun,Ltype,m1fcnstr,p,varargin{:});
        ph = dsolveLCSR(N,LI,LJ,LA,p);
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    else
        ph = p;
    end
    if existU
%         ph = iterapp('mldivide',m2fun,Utype,m2fcnstr,ph,varargin{:});
        ph = dsolveUCSR(N,UI,UJ,UA,ph);
        if ~all(isfinite(ph))
            flag = 2;
            resvec = resvec(1:2*ii-1);
            break
        end
    end
%     v = iterapp('mtimes',afun,atype,afcnstr,ph,varargin{:});
    v = dgemvCSR(N,AI,AJ,AA,ph);
%     rtv = rt' * v;
    rtv = ddot(N,rt,v);
    if (rtv == 0) || isinf(rtv)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    alpha = rho / rtv;
    if isinf(alpha)
        flag = 4;
        resvec = resvec(1:2*ii-1);
        break
    end
    
    if abs(alpha)*norm(ph) < eps*norm(x)
        stag = stag + 1;
    else
        stag = 0;
    end
    
%     xhalf = x + alpha * ph;        % form the "half" iterate
    xhalf = daxpby(N,x,ph,1,alpha);
%     s = r - alpha * v;             % residual associated with xhalf
    s = daxpby(N,r,v, 1, -alpha);
%     normr = norm(s);
    normr = dnorm(N,s);
    normr_act = normr;
    resvec(2*ii) = normr;
    
    % check for convergence
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
%         s = b - iterapp('mtimes',afun,atype,afcnstr,xhalf,varargin{:});
        s = b - dgemvCSR(N,AI,AJ,AA,xhalf);
%         normr_act = norm(s);
        normr_act = dnorm(N,s);
        resvec(2*ii) = normr_act;
        if normr_act <= tolb
            x = xhalf;
            flag = 0;
            iter = ii - 0.5;
            resvec = resvec(1:2*ii);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                x = xhalf;
                resvec = resvec(1:2*ii);
                break;
            end
        end
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii);
        break
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = xhalf;
        imin = ii - 0.5;
    end
    
    if existL
%         sh = iterapp('mldivide',m1fun,Ltype,m1fcnstr,s,varargin{:});
        sh = dsolveLCSR(N,LI,LJ,LA,s);
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    else
        sh = s;
    end
    if existU
%         sh = iterapp('mldivide',m2fun,Utype,m2fcnstr,sh,varargin{:});
        sh = dsolveUCSR(N,UI,UJ,UA,sh);
        if ~all(isfinite(sh))
            flag = 2;
            resvec = resvec(1:2*ii);
            break
        end
    end
%     t = iterapp('mtimes',afun,atype,afcnstr,sh,varargin{:});
    t = dgemvCSR(N,AI,AJ,AA,sh);
%     tt = t' * t;
    tt = ddot(N,t,t);
    if (tt == 0) || isinf(tt)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
%     omega = (t' * s) / tt;
    omega = ddot(N,t,s) / tt;
    if isinf(omega)
        flag = 4;
        resvec = resvec(1:2*ii);
        break
    end
    
    if abs(omega)*norm(sh) < eps*norm(xhalf)
        stag = stag + 1;
    else
        stag = 0;
    end
    
%     x = xhalf + omega * sh;        % x = (x + alpha * ph) + omega * sh
    x = daxpby(N,xhalf, sh, 1, omega);
%     r = s - omega * t;
    r = daxpby(N,s,t,1,-omega);
%     normr = norm(r);
    normr = dnorm(N,r);
    normr_act = normr;
    resvec(2*ii+1) = normr;
    
    % check for convergence        
    if (normr <= tolb || stag >= maxstagsteps || moresteps)
%         r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        r = b - dgemvCSR(N,AI,AJ,AA,x);
%         normr_act = norm(r);
        normr_act = dnorm(N,r);
        resvec(2*ii+1) = normr_act;
        if normr_act <= tolb
            flag = 0;
            iter = ii;
            resvec = resvec(1:2*ii+1);
            break
        else
            if stag >= maxstagsteps && moresteps == 0
                stag = 0;
            end
            moresteps = moresteps + 1;
            if moresteps >= maxmsteps
                if ~warned
                    warning(message('MATLAB:bicgstab:tooSmallTolerance'));
                end
                flag = 3;
                resvec = resvec(1:2*ii+1);
                break;
            end
        end        
    end
    
    if normr_act < normrmin        % update minimal norm quantities
        normrmin = normr_act;
        xmin = x;
        imin = ii;
    end
    
    if stag >= maxstagsteps
        flag = 3;
        resvec = resvec(1:2*ii+1);
        break
    end
    
end                                % for ii = 1 : maxit

if isempty(ii)
    ii = 0;
end

% returned solution is first with minimal residual
if flag == 0
    relres = normr_act / n2b;
else
%     r = b - iterapp('mtimes',afun,atype,afcnstr,xmin,varargin{:});
    r = b - dgemvCSR(N,AI,AJ,AA,xmin);
    if norm(r) <= normr_act
        x = xmin;
        iter = imin;
%         relres = norm(r) / n2b;
        relres = dnorm(N,r) / n2b;
    else
        iter = ii;
        relres = normr_act / n2b;
    end
end

% only display a message if the output flag is not used
if nargout < 2
    itermsg('bicgstab',tol,maxit,ii,flag,iter,relres);
end
