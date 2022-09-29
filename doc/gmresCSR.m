function [x,flag,relres,iter,resvec] = gmres_CSR(N,AI,AJ,AA,b,restart,tol,maxit,...
    LI,LJ,LA,UI,UJ,UA,x,varargin)
% A,b,restart,tol,maxit,M1,M2,x,varargin

% GMRES   Generalized Minimum Residual Method.
% CSR∏Ò Ω

if (nargin < 5)
    error(message('MATLAB:gmres:NumInputs'));
end

% Check matrix and right hand side vector inputs have appropriate sizes
% [M,N] = size(A);
M = N;
if ~isequal(size(b),[M,1])
    error(message('MATLAB:gmres:VectorSize', M));
end

% Assign default values to unspecified parameters
if (nargin < 6) || isempty(restart) || (restart == N)
    restarted = false;
else
    restarted = true;
    restart = max(restart, 0);
end
if (nargin < 7) || isempty(tol)
    tol = 1e-6;
end
warned = 0;
if tol < eps
    warning(message('MATLAB:gmres:tooSmallTolerance'));
    warned = 1;
    tol = eps;
elseif tol >= 1
    warning(message('MATLAB:gmres:tooBigTolerance'));
    warned = 1;
    tol = 1-eps;
end
if (nargin < 8) || isempty(maxit)
    if restarted
        maxit = min(ceil(N/restart),10);
    else
        maxit = min(N,10);
    end
end
maxit = max(maxit, 0);

if restarted
    outer = maxit;
    if restart > N
        warning(message('MATLAB:gmres:tooManyInnerItsRestart',restart, N));
        restart = N;
    end
    inner = restart;
else
    outer = 1;
    if maxit > N
        warning(message('MATLAB:gmres:tooManyInnerItsMaxit',maxit, N));
        maxit = N;
    end
    inner = maxit;
end

% Check for all zero right hand side vector => all zero solution
% n2b = norm(b);                   % Norm of rhs vector, b
n2b = dnorm(N, b); 
if (n2b == 0)                    % if    rhs vector is all zeros
    x = zeros(N,1);              % then  solution is all zeros
    flag = 0;                    % a valid solution has been obtained
    relres = 0;                  % the relative residual is actually 0/0
    iter = [0 0];                % no iterations need be performed
    resvec = 0;                  % resvec(1) = norm(b-A*x) = norm(0)
    if (nargout < 2)
        itermsg('gmres',tol,maxit,0,flag,iter,NaN);
    end
    return;
end

if ((nargin >= 11) && ~isempty(LA))
    existM1 = 1;
    if ~isequal(size(LA),[M,M])
        error(message('MATLAB:gmres:PreConditioner1Size', M));
    end
else
    existM1 = 0;
    m1type = 'matrix';
end

if ((nargin >= 14) && ~isempty(UA))
    existM2 = 1;
    if ~isequal(size(UA),[M,M])
        error(message('MATLAB:gmres:PreConditioner2Size', M));
    end
else
    existM2 = 0;
    m2type = 'matrix';
end

if ((nargin >= 15) && ~isempty(x))
    if ~isequal(size(x),[N,1])
        error(message('MATLAB:gmres:XoSize', N));
    end
else
    x = zeros(N,1);
end

if ((nargin > 15) && strcmp(atype,'matrix') && ...
        strcmp(m1type,'matrix') && strcmp(m2type,'matrix'))
    error(message('MATLAB:gmres:TooManyInputs'));
end

% Set up for the method
flag = 1;
xmin = x;                        % Iterate which has minimal residual so far
imin = 0;                        % "Outer" iteration at which xmin was computed
jmin = 0;                        % "Inner" iteration at which xmin was computed
tolb = tol * n2b;                % Relative tolerance
evalxm = 0;
stag = 0;
moresteps = 0;
maxmsteps = min([floor(N/50),5,N-maxit]);
maxstagsteps = 3;
minupdated = 0;

% x0iszero = (norm(x) == 0);
x0iszero = (dnorm(N,x) == 0);
% r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
r = b - dgemvCSR(N,AI,AJ,AA,x);
% normr = norm(r);                 % Norm of initial residual
normr = dnorm(N, r);
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2b;
    iter = [0 0];
    resvec = normr;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return
end
minv_b = b;

if existM1
    % r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
    r = dsolveLCSR(N,LI,LJ,LA,r);
    if ~x0iszero
        % minv_b = iterapp('mldivide',m1fun,m1type,m1fcnstr,b,varargin{:});
        minv_b = dsolveLCSR(N,LI,LJ,LA,b);
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

if existM2
    % r = iterapp('mldivide',m2fun,m2type,m2fcnstr,r,varargin{:});
    r = dsolveUCSR(N,UI,UJ,UA,r);
    if ~x0iszero
        % minv_b = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_b,varargin{:});
        minv_b = dsolveUCSR(N,UI,UJ,UA,minv_b);
    else
        minv_b = r;
    end
    if ~all(isfinite(r)) || ~all(isfinite(minv_b))
        flag = 2;
        x = xmin;
        relres = normr / n2b;
        iter = [0 0];
        resvec = normr;
        return
    end
end

% normr = norm(r);                 % norm of the preconditioned residual
normr = dnorm(N, r);
% n2minv_b = norm(minv_b);         % norm of the preconditioned rhs
n2minv_b = dnorm(N, minv_b);
clear minv_b;
tolb = tol * n2minv_b;
if (normr <= tolb)               % Initial guess is a good enough solution
    flag = 0;
    relres = normr / n2minv_b;
    iter = [0 0];
    resvec = n2minv_b;
    if (nargout < 2)
        itermsg('gmres',tol,maxit,[0 0],flag,iter,relres);
    end
    return;
end

resvec = zeros(inner*outer+1,1);  % Preallocate vector for norm of residuals
resvec(1) = normr;                % resvec(1) = norm(b-A*x0)
normrmin = normr;                 % Norm of residual from xmin

%  Preallocate J to hold the Given's rotation constants.
J = zeros(2,inner);

U = zeros(N,inner);
R = zeros(inner,inner);
w = zeros(inner+1,1);

for outiter = 1 : outer
    %  Construct u for Householder reflector.
    %  u = r + sign(r(1))*||r||*e1
    u = r;
    % normr = norm(r);
    normr = dnorm(N,r);
    beta = scalarsign(r(1))*normr;
    u(1) = u(1) + beta;
    % u = u / norm(u);
    u = u / dnorm(N, u);
    
    U(:,1) = u;
    
    %  Apply Householder projection to r.
    %  w = r - 2*u*u'*r;
    w(1) = -beta;
    
    for initer = 1 : inner
        %  Form P1*P2*P3...Pj*ej.
        %  v = Pj*ej = ej - 2*u*u'*ej
        v = -2*(u(initer)')*u;
        v(initer) = v(initer) + 1;
        %  v = P1*P2*...Pjm1*(Pj*ej)
        for k = (initer-1):-1:1
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        %  Explicitly normalize v to reduce the effects of round-off.
        % v = v/norm(v);
        v = v/dnorm(N,v);
        
        %  Apply A to v.
        % v = iterapp('mtimes',afun,atype,afcnstr,v,varargin{:});
        v = dgemvCSR(N,AI,AJ,AA,v);
        %  Apply Preconditioner.
        if existM1
            % v = iterapp('mldivide',m1fun,m1type,m1fcnstr,v,varargin{:});
            v = dsolveLCSR(N,LI,LJ,LA,v);
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        
        if existM2
            % v = iterapp('mldivide',m2fun,m2type,m2fcnstr,v,varargin{:});
            v = dsolveUCSR(N,UI,UJ,UA,v);
            if ~all(isfinite(v))
                flag = 2;
                break
            end
        end
        %  Form Pj*Pj-1*...P1*Av.
        for k = 1:initer
            Utemp = U(:,k);
            v = v - Utemp*(2*(Utemp'*v));
        end
        
        %  Determine Pj+1.
        if (initer ~= length(v))
            %  Construct u for Householder reflector Pj+1.
            u = v;
            u(1:initer) = 0;
            % alpha = norm(u);
            alpha = dnorm(N, u);
            if (alpha ~= 0)
                alpha = scalarsign(v(initer+1))*alpha;
                %  u = v(initer+1:end) +
                %        sign(v(initer+1))*||v(initer+1:end)||*e_{initer+1)
                u(initer+1) = u(initer+1) + alpha;
                % u = u / norm(u);
                u = u / dnorm(N, u);
                U(:,initer+1) = u;
                
                %  Apply Pj+1 to v.
                %  v = v - 2*u*(u'*v);
                v(initer+2:end) = 0;
                v(initer+1) = -alpha;
            end
        end
        
        %  Apply Given's rotations to the newly formed v.
        for colJ = 1:initer-1
            tmpv = v(colJ);
            v(colJ)   = conj(J(1,colJ))*v(colJ) + conj(J(2,colJ))*v(colJ+1);
            v(colJ+1) = -J(2,colJ)*tmpv + J(1,colJ)*v(colJ+1);
        end
        
        %  Compute Given's rotation Jm.
        if ~(initer==length(v))
            % rho = norm(v(initer:initer+1));
            rho = dnorm(2,v(initer:initer+1));
            J(:,initer) = v(initer:initer+1)./rho;
            w(initer+1) = -J(2,initer).*w(initer);
            w(initer) = conj(J(1,initer)).*w(initer);
            %
            conj(J(1,initer))
            %
            v(initer) = rho;
            v(initer+1) = 0;
        end
        
        R(:,initer) = v(1:inner);
        
        normr = abs(w(initer+1));
        resvec((outiter-1)*inner+initer+1) = normr;
        normr_act = normr;
        
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
            if evalxm == 0
                ytmp = R(1:initer,1:initer) \ w(1:initer);
                additive = U(:,initer)*(-2*ytmp(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + ytmp(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + ytmp(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                if norm(additive) < eps*norm(x)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                xm = x + additive;
                evalxm = 1;
            elseif evalxm == 1
                addvc = [-(R(1:initer-1,1:initer-1)\R(1:initer-1,initer))*...
                    (w(initer)/R(initer,initer)); w(initer)/R(initer,initer)];
                if norm(addvc) < eps*norm(xm)
                    stag = stag + 1;
                else
                    stag = 0;
                end
                additive = U(:,initer)*(-2*addvc(initer)*conj(U(initer,initer)));
                additive(initer) = additive(initer) + addvc(initer);
                for k = initer-1 : -1 : 1
                    additive(k) = additive(k) + addvc(k);
                    additive = additive - U(:,k)*(2*(U(:,k)'*additive));
                end
                xm = xm + additive;
            end
            % r = b - iterapp('mtimes',afun,atype,afcnstr,xm,varargin{:});
            r = b - dgemvCSR(N,AI,AJ,AA,xm);
            if norm(r) <= tol*n2b
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            end
            minv_r = r;
            if existM1
                % minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
                minv_r = dsolveLCSR(N,LI,LJ,LA,r);
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            if existM2
                % minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
                minv_r = dsolveUCSR(N,UI,UJ,UA,minv_r);
                if ~all(isfinite(minv_r))
                    flag = 2;
                    break
                end
            end
            
            % normr_act = norm(minv_r);
            normr_act = dnorm(N, minv_r);
            resvec((outiter-1)*inner+initer+1) = normr_act;
            
            if normr_act <= normrmin
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                xmin = xm;
                minupdated = 1;
            end
            
            if normr_act <= tolb
                x = xm;
                flag = 0;
                iter = [outiter, initer];
                break
            else
                if stag >= maxstagsteps && moresteps == 0
                    stag = 0;
                end
                moresteps = moresteps + 1;
                if moresteps >= maxmsteps
                    if ~warned
                        warning(message('MATLAB:gmres:tooSmallTolerance'));
                    end
                    flag = 3;
                    iter = [outiter, initer];
                    break;
                end
            end
        end
        
        if normr_act <= normrmin
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
            minupdated = 1;
        end
        
        if stag >= maxstagsteps
            flag = 3;
            break;
        end
    end         % ends inner loop
    
    if isempty(initer)
        initer = 0;
    end
    
    evalxm = 0;
    
    if flag ~= 0
        if minupdated
            idx = jmin;
        else
            idx = initer;
        end
        if idx > 0 % Allow case inner==0 to flow through
            y = R(1:idx,1:idx) \ w(1:idx);
            additive = U(:,idx)*(-2*y(idx)*conj(U(idx,idx)));
            additive(idx) = additive(idx) + y(idx);
            for k = idx-1 : -1 : 1
                additive(k) = additive(k) + y(k);
                additive = additive - U(:,k)*(2*(U(:,k)'*additive));
            end
            x = x + additive;
        end
        xmin = x;
        % r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
        r = b - dgemvCSR(N,AI,AJ,AA,x);
        minv_r = r;
        if existM1
            % minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
            minv_r = dsolveLCSR(N,LI,LJ,LA,r);
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        if existM2
            % minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
            minv_r = dsolveUCSR(N,UI,UJ,UA,minv_r);
            if ~all(isfinite(minv_r))
                flag = 2;
                break
            end
        end
        % normr_act = norm(minv_r);
        normr_act = dnorm(N, minv_r);
        r = minv_r;
    end
    
    if normr_act <= normrmin
        xmin = x;
        normrmin = normr_act;
        imin = outiter;
        jmin = initer;
    end
    
    if flag == 3
        break;
    end
    if normr_act <= tolb
        flag = 0;
        iter = [outiter, initer];
        break;
    end
    minupdated = 0;
end         % ends outer loop

if isempty(outiter)
    outiter = 0;
    initer = 0;
    normr_act = normrmin;
end

% returned solution is that with minimum residual
if flag == 0
    relres = normr_act / n2minv_b;
else
    x = xmin;
    iter = [imin jmin];
    relres = normr_act / n2minv_b;
end

resvec = resvec(1:max(outiter-1,0)*inner+initer+1);
if flag == 2 && initer ~= 0
    resvec(end) = [];
end

% only display a message if the output flag is not used
if nargout < 2
    if restarted
        itermsg(sprintf('gmres(%d)',restart),tol,maxit,[outiter initer],flag,iter,relres);
    else
        itermsg(sprintf('gmres'),tol,maxit,initer,flag,iter(2),relres);
    end
end

function sgn = scalarsign(d)
sgn = sign(d);
if (sgn == 0)
    sgn = 1;
end
