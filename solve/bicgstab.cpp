#include "solve.h"
#include "../blas/blas.h"
#include "../blas/vector.h"
#include "../blas/CSRmatrix.h"
#include "../config/information.h"
#include "../config/config.h"

//#include <mkl.h>

void bicgstab(IterInfo* info, double* x, const CSRMatrix* A, const double* b, const double* Tol, const int* MaxIter,
    const CSRMatrix* L, const CSRMatrix* U, const double* x0)
{
    // BICGSTAB   BiConjugate Gradients Stabilized Method.
    /* 输入的变量一定是已经分配好内存的(info 除外):
    * x: 解, N*1数组, 已经分配好内存
    * A: 系数矩阵, N*N矩阵
    * b: 右边项, N*1数组
    * Tol: 实数容差
    * MaxIter: 最大迭代次数
    * L: 预处理矩阵, 单位下三角矩阵, N*N矩阵
    * U: 预处理矩阵, 上三角矩阵, N*N矩阵
    * x0: 解初始值
    * 
    * flag: 计算结果编号
    *    0 ---- bicgstab在maxit次迭代内收敛至所需容差tol;
    *    1 ---- bicgstab迭代maxit次,但未收敛;
    *    2 ---- 预条件子M是病态的;
    *    3 ---- bicgstab已停滞(两次连续迭代相同);
    *    4 ---- 在执行bicgstab时计算的某个标量太小或太大,以致无法继续计算.
    */
    // -------- 基本定义 -------- //
    int N = 0;
    double tol = 0.0;
    int maxit = 0;
    int flag = 0;// 计算结果编号
    double relres = 0.0;// the relative residual
    float iter = 0.0;// iterations need be performed: 半步
    double* resvec = NULL;// 每个迭代步的残差, resvec = norm(b - A * x)
    bool existL = false, existU = false;
    char Ltype = '0', Utype = '0';
    int warned = 0;// 警告编号
    double n2b = 0.0;// Norm of rhs vector, b

    // -------- 检查输入 -------- //
    if ((A == NULL) && (b == NULL))
    {
        error("BICGSTAB: Not Enough Inputs. Please check coefficient matrix A or right hand side vector b.");
    }
    else
    {
        N = A->N;
    }
    // 输入参数不合适则选用默认值
    if ((Tol == nullptr) || !isfinite(*Tol) || isnan(*Tol))
    {
        
        tol = 1.0e-6;
        warning("BICGSTAB: Tolerance is set to default value.");
    }
    else
    {
        tol = *Tol;
    }
    if (tol < eps)
    {
        warning("BICGSTAB: Too Small Tolerance.");
        warned = 1;
        tol = eps;
    }
    else if (tol > 1.0)
    {
        warning("BICGSTAB: Too Big Tolerance.");
        warned = 1;
        tol = 1 - eps;
    }
    if ((MaxIter == NULL) || (*MaxIter > 1000000))
    {
        maxit = min(N, 20);
    }
    else
    {
        maxit = *MaxIter;
        maxit = max(maxit, 0);
    }
    // 检查是否右边项全为0
    // Check for all zero right hand side vector = > all zero solution
    n2b = dnorm(&N, b);// b的2范数
    if (n2b < eps) // 解全为0
    {
        dzeros(&N, x);// then  solution is all zeros
        flag = 0;// a valid solution has been obtained
        relres = 0;// the relative residual is actually 0 / 0
        iter = 0;// no iterations need be performed
        resvec = (double*)malloc(1 * sizeof(double));// resvec(1) = norm(b - A * x) = norm(0)
        *resvec = 0.0;

        info->initial(1);
        info->flag = flag;
        info->relres = relres;
        info->iter = iter;
        info->resvec[0] = *resvec;
        return;
    }
    // 检查预处理矩阵
    if (L != NULL)
    {
        existL = true;
        Ltype = 'S';
        if (L->N != A->N)
        {
            error("Lower matrix dimension is not equal to coefficient matrix.");
            return;
        }
    }
    else
    {
        existL = false;
        Ltype = '0';
    }
    if (U != NULL)
    {
        existU = true;
        Utype = 'S';
        if (U->N != A->N)
        {
            error("Upper matrix dimension is not equal to coefficient matrix.");
            return;
        }
    }
    else
    {
        existU = false;
        Utype = '0';
    }
    // 检查初始值是否存在
    if (x0 != NULL)
    {
        dcopy(&N, x0, x);
    }
    else
    {
        dzeros(&N, x);
    }

    // ---------------------------------------- //
    // BICGSTAB算法计算求解
    // 变量定义
    const double dPosOne = 1.0, dNegOne = -1.0;
    int Size = 0, NewSize = 0;
    double* xmin = NULL;// Iterate which has minimal residual so far
    float imin = 0;// Iteration at which xmin was computed
    double tolb = 0.0;// Relative tolerance
    double* r = NULL;
    double normr = 0.0;// Norm of residual
    double normr_act = 0.0;

    double* rt = NULL;// Shadow residual
    double normrmin = 0.0;// Norm of residual from xmin
    double rho = 0.0, omega = 0.0;
    int stag = 0;// stagnation of the method
    double alpha = 0.0;// overshadow any functions named alpha
    int moresteps = 0;
    int maxmsteps = 0;
    int maxstagsteps = 0;
    double rho1 = 0.0;
    double beta = 0.0;
    double* p = NULL;
    double* ph = NULL;
    double* v = NULL;
    double* temp = NULL;
    double da = 0.0, db = 0.0;
    double rtv = 0.0;

    double* xhalf = NULL;
    double* s = NULL;
    double* sh = NULL;
    double* t = NULL;
    double tt = 0;

    //
    resvec = (double*)malloc((2 * maxit + 1) * sizeof(double));
    xmin = (double*)malloc(N * sizeof(double));
    r = (double*)malloc(N * sizeof(double));
    if ((resvec == NULL) || (xmin == NULL) || (r == NULL))
    {
        error("内存分配失败!");
        return;
    }

    flag = 1;
    tolb = tol * n2b;
    dcopy(&N, x, xmin);
    // r = b - dgemvCSR(N, AI, AJ, AA, x);
    dgemvCSR(A, x, r);// r = A*x
    daxpby(&N, &dPosOne, b, &dNegOne, r);// r = 1*b - 1 * r
    // normr = norm(r);% Norm of residual
    normr = dnorm(&N, r);
    normr_act = normr;

    if (normr <= tolb) // Initial guess is a good enough solution
    {
        flag = 0;
        relres = normr / n2b;
        iter = 0;
        resvec[0] = normr;

        info->initial(1);
        info->flag = flag;
        info->relres = relres;
        info->iter = iter;
        info->resvec[0] = resvec[0];
        return;
    }

    // ---------------------------------------- //
    rt = (double*)malloc(N * sizeof(double));
    p = (double*)malloc(N * sizeof(double));
    ph = (double*)malloc(N * sizeof(double));
    v = (double*)malloc(N * sizeof(double));
    temp = (double*)malloc(N * sizeof(double));
    if (rt == NULL || p == NULL || ph == NULL || v == NULL || temp == NULL)
    {
        error("内存分配失败!");
        return;
    }
    xhalf = (double*)malloc(N * sizeof(double));
    s = (double*)malloc(N * sizeof(double));
    sh = (double*)malloc(N * sizeof(double));
    t = (double*)malloc(N * sizeof(double));
    if (xhalf == NULL || s == NULL || sh == NULL || t == NULL)
    {
        error("内存分配失败!");
        return;
    }

    dcopy(&N, r, rt);
    Size = 2 * maxit + 1;
    dzeros(&Size, resvec);
    resvec[0] = normr;
    normrmin = normr;
    rho = 1.0;
    omega = 1.0;
    stag = 0;
    alpha = 1.0;
    moresteps = 0;
    maxmsteps = min(min(N / 50, 10), N - maxit);
    maxstagsteps = 3;
    rho1 = 1.0;
    beta = 1.0;
    da = 1.0;
    db = 1.0;

    // loop over maxit iterations(unless convergence or failure)
    int ii = 0;
    for ( ii = 1; ii < maxit; ii++)
    {
        rho1 = rho;
        rho = ddot(&N, rt, r);// rho = rt' * r;
        if ((rho == 0.0) || isinf(rho))
        {
            flag = 4;
            Size = 2 * maxit + 1;
            NewSize = 2 * ii - 1;
            //
            info->initial(NewSize);
            info->flag = flag;
            //info->relres = relres;
            //info->iter = iter;
            dcopy(&NewSize, resvec, info->resvec);//resvec = resvec[1:2 * ii - 1];
            break;
        }

        if (ii == 1)
        {
            dcopy(&N, r, p);
        }
        else
        {
            beta = (rho / rho1) * (alpha / omega);
            if ((abs(beta) < eps) || !isfinite(beta))
            {
                flag = 4;
                break;
            }
            // p = r + beta * (p - omega * v);
            //p = daxpby(N, r, daxpby(N, p, v, beta, -beta * omega), 1, 1);
            da = -beta * omega;
            db = beta;
            daxpby(&N, &da, v, &db, p);
            da = 1;
            db = 1;
            daxpby(&N, &da, r, &db, p);
        }

        if (existL)
        {
            // ph = iterapp('mldivide', m1fun, Ltype, m1fcnstr, p, varargin{ : });
            //ph = solveLCSR(N, LI, LJ, LA, p);
            solveLCSR(L, p, ph);
            if (!isfiniteall(&N, ph))
            {
                flag = 2;
                //resvec = resvec(1:2 * ii - 1);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii - 1;
                info->initial(NewSize);
                info->flag = flag;
                //info->relres = relres;
                //info->iter = iter;
                dcopy(&NewSize, resvec, info->resvec);//resvec = resvec[1:2 * ii - 1];
                break;
            }
        }
        else
        {
            dcopy(&N, p, ph);
        }

        if (existU)
        {
            // ph = iterapp('mldivide', m2fun, Utype, m2fcnstr, ph, varargin{ : });
            //ph = dsolveUCSR(N, UI, UJ, UA, ph);
            solveUCSR(U, ph, temp);
            dcopy(&N, temp, ph);
            if (!isfiniteall(&N, ph))
            {
                flag = 2;
                //resvec = resvec(1:2 * ii - 1);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii - 1;
                info->initial(NewSize);
                info->flag = flag;
                //info->relres = relres;
                //info->iter = iter;
                dcopy(&NewSize, resvec, info->resvec);//resvec = resvec[1:2 * ii - 1];
                break;
            }
        }
        // v = iterapp('mtimes', afun, atype, afcnstr, ph, varargin{ : });

        dgemvCSR(A, ph, v);//v = dgemvCSR(N,AI,AJ,AA,ph);
        // rtv = rt' * v;
        rtv = ddot(&N, rt, v);
        if ((rtv == 0) || isinf(rtv))
        {
            flag = 4;
            //resvec = resvec(1:2 * ii - 1);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii - 1;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }
        alpha = rho / rtv;
        if (isinf(alpha))
        {
            flag = 4;
            //resvec = resvec(1:2 * ii - 1);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii - 1;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }

        if (abs(alpha) * dnorm(&N, ph) < eps * dnorm(&N, x))
        {
            stag += 1;
        }
        else
        {
            stag = 0;
        }

        // xhalf = x + alpha * ph;% form the "half" iterate
        dcopy(&N, ph, xhalf);
        daxpby(&N, &dPosOne, x, &alpha, xhalf);
        // s = r - alpha * v;% residual associated with xhalf
        dscal(&N, &dNegOne, v, s);
        daxpby(&N, &dPosOne, r, &alpha, s);
        // normr = norm(s);
        normr = dnorm(&N, s);
        normr_act = normr;
        resvec[2 * ii - 1] = normr;//resvec(2*ii) = normr;

        // check for convergence
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
        {
            // s = b - iterapp('mtimes', afun, atype, afcnstr, xhalf, varargin{ : });
            dgemvCSR(A, xhalf, s);
            daxpby(&N, &dPosOne, b, &dNegOne, s);
            // normr_act = norm(s);
            normr_act = dnorm(&N, s);
            resvec[2 * ii - 1] = normr_act;
            if (normr_act <= tolb)
            {
                dcopy(&N, xhalf, x);// x = xhalf;
                flag = 0;
                iter = ii - 0.5;
                //resvec = resvec(1:2 * ii);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii;
                info->initial(NewSize);
                info->flag = flag;
                info->iter = iter;
                dcopy(&NewSize, resvec, info->resvec);
                break;
            }
            else
            {
                if ((stag >= maxstagsteps) && (moresteps == 0))
                {
                    stag = 0;
                }
                moresteps += 1;
                if (moresteps >= maxmsteps)
                {
                    if (!warned)
                    {
                        warning("BICGSTAB: Too Small Tolerance.");
                    }
                    flag = 3;
                    dcopy(&N, xhalf, x);// x = xhalf;
                    //resvec = resvec(1:2 * ii);
                    Size = 2 * maxit + 1;
                    NewSize = 2 * ii;
                    info->initial(NewSize);
                    info->flag = flag;
                    dcopy(&NewSize, resvec, info->resvec);
                    break;
                }
            }
        }

        if (stag >= maxstagsteps)
        {
            flag = 3;
            //resvec = resvec(1:2 * ii);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }

        if (normr_act < normrmin) // update minimal norm quantities
        {
            normrmin = normr_act;
            dcopy(&N, xhalf, xmin);// xmin = xhalf;
            imin = ii - 0.5;
        }

        if (existL)
        {
            // sh = iterapp('mldivide', m1fun, Ltype, m1fcnstr, s, varargin{ : });
            solveLCSR(L, s, sh);
            if (!isfiniteall(&N, sh))
            {
                flag = 2;
                //resvec = resvec(1:2 * ii);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii;
                info->initial(NewSize);
                info->flag = flag;
                dcopy(&NewSize, resvec, info->resvec);
                break;
            }
        }
        else
        {
            dcopy(&N, s, sh);// sh = s;
        }

        if (existU)
        {
            // sh = iterapp('mldivide', m2fun, Utype, m2fcnstr, sh, varargin{ : });
            solveUCSR(U, sh, temp);
            dcopy(&N, temp, sh);
            if (!isfiniteall(&N, sh))
            {
                flag = 2;
                //resvec = resvec(1:2 * ii);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii;
                info->initial(NewSize);
                info->flag = flag;
                dcopy(&NewSize, resvec, info->resvec);
                break;
            }
        }
        // t = iterapp('mtimes', afun, atype, afcnstr, sh, varargin{ : });
        dgemvCSR(A, sh, t);
        // tt = t' * t;
        tt = ddot(&N, t, t);
        if ((abs(tt) < eps) || isinf(tt))
        {
            flag = 4;
            //resvec = resvec(1:2 * ii);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }
        // omega = (t' * s) / tt;
        omega = ddot(&N, t, s) / tt;
        if (isinf(omega))
        {
            flag = 4;
            //resvec = resvec(1:2 * ii);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }

        if (abs(omega) * dnorm(&N, sh) < eps * dnorm(&N, xhalf))
        {
            stag += 1;
        }
        else
        {
            stag = 0;
        }

        // x = xhalf + omega * sh;% x = (x + alpha * ph) + omega * sh
        dcopy(&N, sh, x);
        daxpby(&N, &dPosOne, xhalf, &omega, x);
        // r = s - omega * t;
        dscal(&N, &dNegOne, t, r);
        daxpby(&N, &dPosOne, s, &omega, r);
        // normr = norm(r);
        normr = dnorm(&N, r);
        normr_act = normr;
        resvec[2 * ii] = normr;;// resvec(2 * ii + 1) = normr;

        // check for convergence
        if (normr <= tolb || stag >= maxstagsteps || moresteps)
        {
            // r = b - iterapp('mtimes', afun, atype, afcnstr, x, varargin{ : });
            dgemvCSR(A, x, r);
            daxpby(&N, &dPosOne, b, &dNegOne, r);
            // normr_act = norm(r);
            normr_act = dnorm(&N, r);
            resvec[2 * ii + 1] = normr_act; //resvec(2 * ii + 1) = normr_act;
            if (normr_act <= tolb)
            {
                flag = 0;
                iter = ii;
                //resvec = resvec(1:2 * ii + 1);
                Size = 2 * maxit + 1;
                NewSize = 2 * ii + 1;
                info->initial(NewSize);
                info->flag = flag;
                info->iter = iter;
                dcopy(&NewSize, resvec, info->resvec);
                break;
            }
            else
            {
                if (stag >= maxstagsteps && moresteps == 0)
                {
                    stag = 0;
                }
                moresteps += 1;
                if (moresteps >= maxmsteps)
                {
                    if (!warned)
                    {
                        warning("BICGSTAB: Too Small Tolerance.");
                    }
                    flag = 3;
                    //resvec = resvec(1:2 * ii + 1);
                    Size = 2 * maxit + 1;
                    NewSize = 2 * ii + 1;
                    info->initial(NewSize);
                    info->flag = flag;
                    dcopy(&NewSize, resvec, info->resvec);
                    break;
                }
            }
        }

        if (normr_act < normrmin) // update minimal norm quantities
        {
            normrmin = normr_act;
            dcopy(&N, x, xmin);// xmin = x;
            imin = ii;
        }

        if (stag >= maxstagsteps)
        {
            flag = 3;
            //resvec = resvec(1:2 * ii + 1);
            Size = 2 * maxit + 1;
            NewSize = 2 * ii + 1;
            info->initial(NewSize);
            info->flag = flag;
            dcopy(&NewSize, resvec, info->resvec);
            break;
        }

    } // for ii = 1 : maxit

    printf("bicgstab iteration: %d\n", ii);
    
    if (ii < 1)
    {
        ii = 0;
    }
    
    // returned solution is first with minimal residual
    if (flag == 0)
    {
        relres = normr_act / n2b;
        info->flag = flag;
        info->iter = iter;
        info->relres = relres;
        NewSize = 2 * iter + 1;
        info->initial(NewSize);
        dcopy(&NewSize, resvec, info->resvec);
    }
    else
    {
        // r = b - iterapp('mtimes', afun, atype, afcnstr, xmin, varargin{ : });
        dgemvCSR(A, xmin, r);
        daxpby(&N, &dPosOne, b, &dNegOne, r);
        if (dnorm(&N, r) <= normr_act)
        {
            dcopy(&N, xmin, x);// x = xmin;
            iter = imin;
            // relres = norm(r) / n2b;
            relres = dnorm(&N, r) / n2b;
        }
        else
        {
            iter = ii;
            relres = normr_act / n2b;
        }
        info->flag = flag;
        info->iter = iter;
        info->relres = relres;
        NewSize = 2 * ii + 1;
        info->initial(NewSize);
        dcopy(&NewSize, resvec, info->resvec);
    }
    
    // -------- 释放内存 -------- //
    free(resvec);
    free(xmin);
    free(r);
    free(rt);
    free(p);
    free(ph);
    free(v);
    free(temp);
    free(xhalf);
    free(s);
    free(sh);
    free(t);
}


/*
int MKL_fgmres(double relativeErr) 
{
    //
    // Solve sparse equation using GMRES iterative method.
    // Coefficient matrix with CSR format
    //
    // Currently fgmres/ilut supported in MKL only in sequential mode.
    //
    // 2013-May-1
    //
    const int size = 20, N = 40;
    MKL_INT ia[5] = { 1, 4, 7, 10, 13 }, ibilut[5];
    MKL_INT ja[12] = { 1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4 }, jbilut[16];
    double A[12] = { 4., -1., -1., -1., 4., -1., -1., 4., -1., -1., -1., 4. };

    MKL_INT ipar[size];
    double dpar[size], tmp[N * (2 * N + 1) + (N * (N + 9)) / 2 + 1];
    double trvec[N], bilut[12];
    double expected_solution[N] = { 1.0, 1.0, 1.0, 1.0 };
    double rhs[N], b[N];
    double computed_solution[N];
    double residual[N];

    MKL_INT itercount, ierr = 0;
    MKL_INT RCI_request, i, ivar;
    double dvar;
    char cvar, cvar1, cvar2;
    MKL_INT matsize = 10, incx = 1, ref_nit = 4;
    double ref_norm2 = 7.836719E+0, nrm2;
    MKL_INT maxfil;
    double tol;

    printf("--------------------------------------------------------\n");
    printf("The FULLY ADVANCED example RCI FGMRES with ILUT\n");
    printf("preconditioner to solve \n");
    printf("the non-degenerate algebraic system of linear equations\n");
    printf("--------------------------------------------------------\n\n");

    ivar = N;
    cvar = 'N';

    mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, expected_solution, rhs);
    i = 1;
    dcopy(&ivar, rhs, &i, b, &i);
    for (i = 0; i < N; i++)
    {
        computed_solution[i] = 0.0;
    }
    computed_solution[0] = 100.0;

    dfgmres_init(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request != 0)
    {
        goto FAILED;
    }
    ipar[30] = 1;
    dpar[30] = 1.E-5;
    tol = 1.E-6;
    maxfil = 1;

    dcsrilut(&ivar, A, ia, ja, bilut, ibilut, jbilut, &tol, &maxfil, ipar, dpar, &ierr);

    nrm2 = dnrm2(&matsize, bilut, &incx);
    if (ierr != 0)
    {
        printf("Preconditioner dcsrilut has returned the ERROR code % d\n", ierr);
        goto FAILED;
    }
    ipar[14] = 2;
    ipar[7] = 0;
    ipar[10] = 1;
    dpar[0] = 1.0E-3;

    dfgmres_check(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request != 0)
    {
        goto FAILED;
    }

ONE:
    dfgmres(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp);
    if (RCI_request == 0)
    {
        goto COMPLETE;
    }
    if (RCI_request == 1)
    {
        mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
        goto ONE;
    }
    if (RCI_request == 2)
    {
        ipar[12] = 1;
        dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
        mkl_dcsrgemv(&cvar, &ivar, A, ia, ja, b, residual);
        dvar = -1.0E0;
        i = 1;
        daxpy(&ivar, &dvar, rhs, &i, residual, &i);
        dvar = dnrm2(&ivar, residual, &i);
        if (dvar < 1.0E-3)
        {
            goto COMPLETE;
        }
        else
        {
            goto ONE;
        }
    }
    if (RCI_request == 3)
    {
        cvar1 = 'L';
        cvar = 'N';
        cvar2 = 'U';
        mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, bilut, ibilut, jbilut, &tmp[ipar[21] - 1], trvec);
        cvar1 = 'U';
        cvar = 'N';
        cvar2 = 'N';
        mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, bilut, ibilut, jbilut, trvec, &tmp[ipar[22] - 1]);
        goto ONE;
    }
    if (RCI_request == 4)
    {
        if (dpar[6] < 1.0E-12)
        {
            goto COMPLETE;
        }
        else
        {
            goto ONE;
        }
    }
    else
    {
        goto FAILED;
    }

COMPLETE:
    ipar[12] = 0;
    dfgmres_get(&ivar, computed_solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
    printf("The system has been solved \n");
    printf("\nNumber of iterations: %d\n", itercount);
    printf("\n");
    MKL_Free_Buffers();
    if (itercount == ref_nit && fabs(ref_norm2 - nrm2) < 1.e-6)
    {
        printf("--------------------------------------------------------\n");
        printf("C example of FGMRES with ILUT preconditioner \n");
        printf("has successfully PASSED all stages of computations\n");
        printf("--------------------------------------------------------\n");
        return 0;
    }

FAILED:
    printf("Calculate Failed for ILUT+FGMRES.\n");
    MKL_Free_Buffers();
    return 1;
}

*/

