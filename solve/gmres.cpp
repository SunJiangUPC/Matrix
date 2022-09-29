#include "solve.h"
#include "../blas/blas.h"
#include "../blas/vector.h"
#include "../blas/CSRmatrix.h"
#include "../config/information.h"
#include "../config/config.h"

static int scalarsign(double d)
{
    if (signbit(d))
    {
        return -1;// 为负数
    }
    return 1;// 非负数
}

// gmres algorithm
void gmres(IterInfo* info, double* x, const CSRMatrix* A, const double* b, const double* Tol, const int* MaxIter, 
    const int* Restart, const CSRMatrix* L, const CSRMatrix* U, const double* x0)
{
    // GMRES   Generalized Minimum Residual Method.
    /* 输入的变量一定是已经分配好内存的(info 除外):
     * x: 解, N*1数组, 已经分配好内存
     * A: 系数矩阵, N*N矩阵
     * b: 右边项, N*1数组
     * Restart: 重启动次数
     * Tol: 实数容差
     * MaxIter: 最大迭代次数
     * L: 预处理矩阵, 单位下三角矩阵, N*N矩阵
     * U: 预处理矩阵, 上三角矩阵, N*N矩阵
     * x0: 解初始值
     *
     * flag: 计算结果编号
     *    0 ---- gmres在maxit次外迭代内收敛至所需公差tol。
     *    1 ---- gmres迭代maxit次，但未收敛。
     *    2 ---- 预条件子M是病态的。
     *    3 ---- gmres已停滞。（两次连续迭代相同。）
     */
     // -------- 基本定义 -------- //
    double tol = 0.0;
    int maxit = 0;
    int restart = 0;
    int N = 0;
    double n2b = 0.0;
    int flag = 1;
    double relres = 0;// the relative residual
    int iter[] = { 0,0 };// iterations: outiter, initer
    double* resvec = NULL;// resvec(1) = norm(b - A * x) = norm(0)
    int warned = 0;
    bool restarted = false;
    double* xmin = x;// Iterate which has minimal residual so far
    int imin = 0; // "Outer" iteration at which xmin was computed
    int jmin = 0; // "Inner" iteration at which xmin was computed
    double tolb = tol * n2b;// Relative tolerance
    int evalxm = 0;
    int stag = 0;
    int moresteps = 0;
    int maxmsteps = min(min(N / 50, 5), N - maxit);// min([floor(N/50),5,N-maxit]);
    int maxstagsteps = 3;
    int minupdated = 0;
    bool x0iszero = (dnorm(&N, x) < eps) ? true : false;// x0iszero = (dnorm(N, x) == 0);
    int outer = 0;
    int inner = 0;
    bool existM1 = false;
    bool existM2 = false;
    const double dNegOne = -1.0, dPosOne = 1.0;

    double* r = NULL;
    double normr = 0.0;
    int nSize = 0;
    int nSizeNew = 0;
    
    // Check matrixand right hand side vector inputs have appropriate sizes
    // Assign default values to unspecified parameters
    
    if (Restart == NULL || (*Restart == N)) {
        restarted = false;
    }
    else {
        restarted = true;
        restart = max(*Restart, 0);
    }
    
    if (*Tol == NULL) {
        tol = 1.0e-6;
    }
    else {
        tol = *Tol;
    }

    if (tol < eps) {
        warning("GMRES: Too Small Tolerance.");
        warned = 1;
        tol = eps;
    }
    else if (tol >= 1) {
        warning("GMRES: Too Big Tolerance.");
        warned = 1;
        tol = 1 - eps;
    }
    if (*MaxIter == NULL) {
        if (restarted) {
            // maxit = min(ceil(N / restart), 10);
            if (N % restart == 0) {
                maxit = N / restart;
            }
            else {
                maxit = N / restart + 1;
            }
            maxit = min(maxit, 10);
        }
        else {
            maxit = min(N, 10);
        }
    }
    maxit = max(maxit, 0);

    if (restarted) {
        outer = maxit;
        if (restart > N) {
            warning("GMRES: Too Many Inner Iteration Restart.");
            restart = N;
        }
        inner = restart;
    }
    else {
        outer = 1;
        if (maxit > N) {
            warning("GMRES: Too Many Inner Iteration Maxit.");
            maxit = N;
        }
        inner = maxit;
    }
    // Check for all zero right hand side vector = > all zero solution
    // n2b = norm(b);// Norm of rhs vector, b
    n2b = dnorm(&N, b);
    if (n2b == 0) { // if    rhs vector is all zeros
        //x = zeros(N, 1);// then  solution is all zeros
        dzeros(&N, x);
        flag = 0;// a valid solution has been obtained
        relres = 0;// the relative residual is actually 0 / 0
        //iter = [0 0];// no iterations need be performed
        iter[0] = 0; iter[1] = 0;
        resvec = 0;// resvec(1) = norm(b - A * x) = norm(0)
        return;
    }
    
    if (L !=NULL ) {
        existM1 = true;
        if (L->N != A->N) {
            error("GMRES: PreConditioner L Size is Wrong.");
        }
    }
    else {
        existM1 = false;
    }
    
    if (U != NULL) {
        existM2 = true;
        if (U->N != A->N) {
            error("GMRES: PreConditioner U Size is Wrong.");
        }
    }
    else {
        existM2 = 0;
    }
    
    if (x0 != NULL) {
        dcopy(&N, x0, x);
    }
    else{
        //x = zeros(N, 1);
        dzeros(&N, x);
    }

    
    // Set up for the method
    flag = 1;
    xmin = x;// Iterate which has minimal residual so far
    imin = 0; // "Outer" iteration at which xmin was computed
    jmin = 0; // "Inner" iteration at which xmin was computed
    tolb = tol * n2b;// Relative tolerance
    evalxm = 0;
    stag = 0;
    moresteps = 0;
    //maxmsteps = min([floor(N / 50), 5, N - maxit]);
    maxmsteps = min(min(N / 50, 5), N - maxit);
    maxstagsteps = 3;
    minupdated = 0;

    // x0iszero = (norm(x) == 0);
    if (dnorm(&N, x) < eps)
    {
        x0iszero = true;
    }
    else {
        x0iszero = false;
    }
    // r = b - iterapp('mtimes', afun, atype, afcnstr, x, varargin{ : });
    //r = b - dgemvCSR(N, AI, AJ, AA, x);
    dgemvCSR(A, x, r);
    daxpby(&N, &dPosOne, b, &dNegOne, r);

    // normr = norm(r);% Norm of initial residual
    normr = dnorm(&N, r);
    if (normr <= tolb) {// Initial guess is a good enough solution
        flag = 0;
        relres = normr / n2b;
        //iter = [0 0];
        iter[0] = 0; iter[1] = 0;
        resvec = normr;
        if (nargout < 2) {
            itermsg('gmres', tol, maxit, [0 0], flag, iter, relres);
        }
        return;
    }
    
    minv_b = b;

    if (existM1) {
        // r = iterapp('mldivide', m1fun, m1type, m1fcnstr, r, varargin{ : });
        r = dsolveLCSR(N, LI, LJ, LA, r);
        if (~x0iszero) {
            // minv_b = iterapp('mldivide', m1fun, m1type, m1fcnstr, b, varargin{ : });
            minv_b = dsolveLCSR(N, LI, LJ, LA, b);
        }
        else {
            minv_b = r;
        }
        if (~all(isfinite(r)) || ~all(isfinite(minv_b))) {
            flag = 2;
            x = xmin;
            relres = normr / n2b;
            iter = [0 0];
            resvec = normr;
            return;
        }
    }
    
    if (existM2) {
        // r = iterapp('mldivide', m2fun, m2type, m2fcnstr, r, varargin{ : });
        r = dsolveUCSR(N, UI, UJ, UA, r);
        if (~x0iszero) {
            // minv_b = iterapp('mldivide', m2fun, m2type, m2fcnstr, minv_b, varargin{ : });
            minv_b = dsolveUCSR(N, UI, UJ, UA, minv_b);
        }
        else {
            minv_b = r;
        }
        if (~all(isfinite(r)) || ~all(isfinite(minv_b))) {
            flag = 2;
            x = xmin;
            relres = normr / n2b;
            iter = [0 0];
            resvec = normr;
            return;
        }
    }
    
    // normr = norm(r);% norm of the preconditioned residual
    normr = dnorm(N, r);
    // n2minv_b = norm(minv_b);% norm of the preconditioned rhs
    n2minv_b = dnorm(N, minv_b);
    clear minv_b;
    tolb = tol * n2minv_b;
    if (normr <= tolb) {// Initial guess is a good enough solution
        flag = 0;
        relres = normr / n2minv_b;
        iter = [0 0];
        resvec = n2minv_b;
        if (nargout < 2) {
            itermsg('gmres', tol, maxit, [0 0], flag, iter, relres);
        }
        return;
    }
    
    resvec = zeros(inner * outer + 1, 1);// Preallocate vector for norm of residuals
    resvec(1) = normr;// resvec(1) = norm(b - A * x0)
    normrmin = normr;// Norm of residual from xmin
    // Preallocate J to hold the Given's rotation constants.
    J = zeros(2, inner);
    U = zeros(N, inner);
    R = zeros(inner, inner);
    w = zeros(inner + 1, 1);
    for (outiter = 1 : outer) {
        // Construct u for Householder reflector.
        // u = r + sign(r(1)) * || r || *e1
        u = r;
        // normr = norm(r);
        normr = dnorm(N, r);
        beta = scalarsign(r(1)) * normr;
        u(1) = u(1) + beta;
        // u = u / norm(u);
        u = u / dnorm(N, u);
        U(:, 1) = u;
        // Apply Householder projection to r.
        // w = r - 2 * u * u'*r;
        w(1) = -beta;

        for (initer = 1 : inner) {
            // Form P1 * P2 * P3...Pj * ej.
            // v = Pj * ej = ej - 2 * u * u'*ej
            v = -2 * (u(initer)')*u;
                v(initer) = v(initer) + 1;
            // v = P1 * P2*...Pjm1 * (Pj * ej)
            for (k = (initer - 1) : -1 : 1) {
                Utemp = U(:, k);
                v = v - Utemp * (2 * (Utemp'*v));
            }

            // Explicitly normalize v to reduce the effects of round - off.
            // v = v / norm(v);
            v = v / dnorm(N, v);

            // Apply A to v.
            // v = iterapp('mtimes', afun, atype, afcnstr, v, varargin{ : });
            v = dgemvCSR(N, AI, AJ, AA, v);
            // Apply Preconditioner.
            if (existM1) {
                // v = iterapp('mldivide', m1fun, m1type, m1fcnstr, v, varargin{ : });
                v = dsolveLCSR(N, LI, LJ, LA, v);
                if (~all(isfinite(v))) {
                    flag = 2;
                    break;
                }
            }

            if (existM2) {
                // v = iterapp('mldivide', m2fun, m2type, m2fcnstr, v, varargin{ : });
                v = dsolveUCSR(N, UI, UJ, UA, v);
                if (~all(isfinite(v))) {
                    flag = 2;
                    break;
                }
            }

            // Form Pj* Pj - 1*...P1 * Av.
            for (k = 1 : initer) {
                Utemp = U(:, k);
                v = v - Utemp * (2 * (Utemp'*v));
            }

            // Determine Pj + 1.
            if (initer ~= length(v)) {
                // Construct u for Householder reflector Pj + 1.
                u = v;
                u(1:initer) = 0;
                // alpha = norm(u);
                alpha = dnorm(N, u);
                if (alpha ~= 0) {

                    alpha = scalarsign(v(initer + 1)) * alpha;
                    // u = v(initer + 1:end) +
                        // sign(v(initer + 1)) * || v(initer + 1:end) || *e_{ initer + 1)
                    u(initer + 1) = u(initer + 1) + alpha;
                    // u = u / norm(u);
                    u = u / dnorm(N, u);
                    U(:, initer + 1) = u;

                    // Apply Pj + 1 to v.
                    // v = v - 2 * u * (u'*v);
                    v(initer + 2:end) = 0;
                    v(initer + 1) = -alpha;
                }
            }

            // Apply Given's rotations to the newly formed v.
            for (colJ = 1 : initer - 1) {
                tmpv = v(colJ);
                v(colJ) = conj(J(1, colJ)) * v(colJ) + conj(J(2, colJ)) * v(colJ + 1);
                v(colJ + 1) = -J(2, colJ) * tmpv + J(1, colJ) * v(colJ + 1);
            }

            // Compute Given's rotation Jm.
            if (~(initer == length(v))) {
                // rho = norm(v(initer:initer + 1));
                rho = dnorm(2, v(initer:initer + 1));
                J(:, initer) = v(initer:initer + 1). / rho;
                w(initer + 1) = -J(2, initer).*w(initer);
                w(initer) = conj(J(1, initer)).*w(initer);
                //
                conj(J(1, initer))
                    //
                    v(initer) = rho;
                v(initer + 1) = 0;
            }

            R(:, initer) = v(1:inner);

            normr = abs(w(initer + 1));
            resvec((outiter - 1) * inner + initer + 1) = normr;
            normr_act = normr;

            if (normr <= tolb || stag >= maxstagsteps || moresteps) {
                if (evalxm == 0) {
                    ytmp = R(1:initer, 1 : initer) \ w(1:initer);
                    additive = U(:, initer) * (-2 * ytmp(initer) * conj(U(initer, initer)));
                    additive(initer) = additive(initer) + ytmp(initer);
                    for (k = initer - 1 : -1 : 1) {
                        additive(k) = additive(k) + ytmp(k);
                        additive = additive - U(:, k) * (2 * (U(:, k)'*additive));
                    }
                    if (norm(additive) < eps * norm(x)) {
                        stag = stag + 1;
                    }
                    else {
                        stag = 0;
                    }
                    xm = x + additive;
                    evalxm = 1;
                }
                else if (evalxm == 1) {
                    addvc = [-(R(1:initer - 1, 1 : initer - 1)\R(1:initer - 1, initer)) *
                        (w(initer) / R(initer, initer)); w(initer) / R(initer, initer)];
                    if (norm(addvc) < eps * norm(xm)) {
                        stag = stag + 1;
                    }
                    else {
                        stag = 0;
                    }
                    additive = U(:, initer) * (-2 * addvc(initer) * conj(U(initer, initer)));
                    additive(initer) = additive(initer) + addvc(initer);
                    for (k = initer - 1 : -1 : 1) {
                        additive(k) = additive(k) + addvc(k);
                        additive = additive - U(:, k) * (2 * (U(:, k)'*additive));
                    }
                    xm = xm + additive;
                }

                // r = b - iterapp('mtimes',afun,atype,afcnstr,xm,varargin{:});
                r = b - dgemvCSR(N, AI, AJ, AA, xm);
                if (norm(r) <= tol * n2b) {
                    x = xm;
                    flag = 0;
                    iter = [outiter, initer];
                    break;
                }
                minv_r = r;
                if (existM1) {
                    // minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
                    minv_r = dsolveLCSR(N, LI, LJ, LA, r);
                    if (~all(isfinite(minv_r))) {
                        flag = 2;
                        break;
                    }
                }

                if (existM2) {
                    // minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
                    minv_r = dsolveUCSR(N, UI, UJ, UA, minv_r);
                    if (~all(isfinite(minv_r))) {
                        flag = 2;
                        break;
                    }
                }

                // normr_act = norm(minv_r);
                normr_act = dnorm(N, minv_r);
                resvec((outiter - 1) * inner + initer + 1) = normr_act;

                if (normr_act <= normrmin) {
                    normrmin = normr_act;
                    imin = outiter;
                    jmin = initer;
                    xmin = xm;
                    minupdated = 1;
                }

                if (normr_act <= tolb) {
                    x = xm;
                    flag = 0;
                    iter = [outiter, initer];
                    break;
                }
                else {
                    if (stag >= maxstagsteps && moresteps == 0) {
                        stag = 0;
                    }
                    moresteps = moresteps + 1;
                    if (moresteps >= maxmsteps) {
                        if (~warned) {
                            warning(message('MATLAB:gmres:tooSmallTolerance'));
                        }
                        flag = 3;
                        iter = [outiter, initer];
                        break;
                    }
                }
            }

            if (normr_act <= normrmin) {
                normrmin = normr_act;
                imin = outiter;
                jmin = initer;
                minupdated = 1;
            }

            if (stag >= maxstagsteps) {
                flag = 3;
                break;
            }
        } // ends inner loop

        if (isempty(initer)) {
            initer = 0;
        }

        evalxm = 0;

        if (flag ~= 0) {
            if (minupdated) {
                idx = jmin;
            }
            else {
                idx = initer;
            }
            if (idx > 0) { // Allow case inner == 0 to flow through
                y = R(1:idx, 1 : idx) \ w(1:idx);
                additive = U(:, idx) * (-2 * y(idx) * conj(U(idx, idx)));
                additive(idx) = additive(idx) + y(idx);
                for (k = idx - 1 : -1 : 1) {
                    additive(k) = additive(k) + y(k);
                    additive = additive - U(:, k) * (2 * (U(:, k)'*additive));
                }
                x = x + additive;
            }
            xmin = x;
            // r = b - iterapp('mtimes',afun,atype,afcnstr,x,varargin{:});
            r = b - dgemvCSR(N, AI, AJ, AA, x);
            minv_r = r;
            if (existM1) {
                // minv_r = iterapp('mldivide',m1fun,m1type,m1fcnstr,r,varargin{:});
                minv_r = dsolveLCSR(N, LI, LJ, LA, r);
                if (~all(isfinite(minv_r))) {
                    flag = 2;
                    break;
                }
            }
            if (existM2) {
                // minv_r = iterapp('mldivide',m2fun,m2type,m2fcnstr,minv_r,varargin{:});
                minv_r = dsolveUCSR(N, UI, UJ, UA, minv_r);
                if (~all(isfinite(minv_r))) {
                    flag = 2;
                    break
                }
            }
            // normr_act = norm(minv_r);
            normr_act = dnorm(N, minv_r);
            r = minv_r;
        }

        if (normr_act <= normrmin) {
            xmin = x;
            normrmin = normr_act;
            imin = outiter;
            jmin = initer;
        }

        if (flag == 3) {
            break;
        }
        if (normr_act <= tolb) {
            flag = 0;
            iter = [outiter, initer];
            break;
        }
        minupdated = 0;
    }// ends outer loop

    if (isempty(outiter)) {
        outiter = 0;
        initer = 0;
        normr_act = normrmin;
    }
    
    // returned solution is that with minimum residual
    if (flag == 0) {
        relres = normr_act / n2minv_b;
    }
    else {
        x = xmin;
        iter = [imin jmin];
        relres = normr_act / n2minv_b;
    }
    
    resvec = resvec(1:max(outiter - 1,0) * inner + initer + 1);
    if (flag == 2 && initer ~= 0) {
        resvec(end) = [];
    }
    
    // only display a message if the output flag is not used
    if (nargout < 2) {
        if (restarted) {
            itermsg(sprintf('gmres(%d)', restart), tol, maxit, [outiter initer], flag, iter, relres);
        }
        else {
            itermsg(sprintf('gmres'), tol, maxit, initer, flag, iter(2), relres);
        }
    }
}

