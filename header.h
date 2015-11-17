#ifndef DLIBFUNCS_H_INCLUDED
#define DLIBFUNCS_H_INCLUDED

#include<iostream>
#include <sstream>
#include<fstream>
#include<string>
#include<cmath>
#include<dlib/matrix.h>
//#include<dlib/matrix_la_abstract.h>
#include<dlib/optimization.h>
#include<dlib/rand.h>

#include <scythestat/rng/mersenne.h>
#include <scythestat/distributions.h>
#include <scythestat/ide.h>
#include <scythestat/la.h>
#include <scythestat/matrix.h>
#include <scythestat/rng.h>
#include <scythestat/smath.h>
#include <scythestat/stat.h>
#include <scythestat/optimize.h>
#include <time.h>

using namespace dlib;
using namespace scythe;

const double BIGNUM = std::numeric_limits<double>::infinity();
const double SigLevel = 0.05;
// const double pi=std::asin(1)*2;
// dlib::pi
typedef matrix<double,0,1> dbcolvec;
typedef matrix<double> dbmat;
typedef matrix<int> intmat;


// ############# optimization routines equipped with Convergence check
class bobyqa_implementationConvgFail
    {
        typedef long integer;
        typedef double doublereal;

    public:

        template <
            typename funct,
            typename T,
            typename U
            >
        double find_min (
            const funct& f,
            T& x,
            long npt,
            const U& xl_,
            const U& xu_,
            const double rhobeg,
            const double rhoend,
            const long max_f_evals,  int& ConvgFail
        ) const
        {
            const unsigned long n = x.size();
            const unsigned long w_size = (npt+5)*(npt+n)+3*n*(n+5)/2;
            scoped_ptr<doublereal[]> w(new doublereal[w_size]);

            // make these temporary matrices becuse U might be some
            // kind of matrix_exp that doesn't support taking the address
            // of the first element.
            matrix<double,0,1> xl(xl_);
            matrix<double,0,1> xu(xu_);


            return bobyqa_ (f,
                            x.size(),
                            npt,
                            &x(0),
                            &xl(0),
                            &xu(0),
                            rhobeg,
                            rhoend,
                            max_f_evals,
                            w.get(), ConvgFail );
        }

    private:


        template <typename funct>
        doublereal bobyqa_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const doublereal rhoend,
            const integer maxfun,
            doublereal *w,  int& ConvgFail
        ) const
        {

            /* System generated locals */
            integer i__1;
            doublereal d__1, d__2;

            /* Local variables */
            integer j, id_, np, iw, igo, ihq, ixb, ixa, ifv, isl, jsl, ipq, ivl, ixn, ixo, ixp, isu, jsu, ndim;
            doublereal temp, zero;
            integer ibmat, izmat;


            /*     This subroutine seeks the least value of a function of many variables, */
            /*     by applying a trust region method that forms quadratic models by */
            /*     interpolation. There is usually some freedom in the interpolation */
            /*     conditions, which is taken up by minimizing the Frobenius norm of */
            /*     the change to the second derivative of the model, beginning with the */
            /*     zero matrix. The values of the variables are constrained by upper and */
            /*     lower bounds. The arguments of the subroutine are as follows. */

            /*     N must be set to the number of variables and must be at least two. */
            /*     NPT is the number of interpolation conditions. Its value must be in */
            /*       the interval [N+2,(N+1)(N+2)/2]. Choices that exceed 2*N+1 are not */
            /*       recommended. */
            /*     Initial values of the variables must be set in X(1),X(2),...,X(N). They */
            /*       will be changed to the values that give the least calculated F. */
            /*     For I=1,2,...,N, XL(I) and XU(I) must provide the lower and upper */
            /*       bounds, respectively, on X(I). The construction of quadratic models */
            /*       requires XL(I) to be strictly less than XU(I) for each I. Further, */
            /*       the contribution to a model from changes to the I-th variable is */
            /*       damaged severely by rounding errors if XU(I)-XL(I) is too small. */
            /*     RHOBEG and RHOEND must be set to the initial and final values of a trust */
            /*       region radius, so both must be positive with RHOEND no greater than */
            /*       RHOBEG. Typically, RHOBEG should be about one tenth of the greatest */
            /*       expected change to a variable, while RHOEND should indicate the */
            /*       accuracy that is required in the final values of the variables. An */
            /*       error return occurs if any of the differences XU(I)-XL(I), I=1,...,N, */
            /*       is less than 2*RHOBEG. */
            /*     MAXFUN must be set to an upper bound on the number of calls of CALFUN. */
            /*     The array W will be used for working space. Its length must be at least */
            /*       (NPT+5)*(NPT+N)+3*N*(N+5)/2. */

            /* Parameter adjustments */
            --w;
            --xu;
            --xl;
            --x;

            /* Function Body */
            np = n + 1;

            /*     Return if the value of NPT is unacceptable. */
            if (npt < n + 2 || npt > (n + 2) * np / 2) {
    throw bobyqa_failure("Return from BOBYQA because NPT is not in the required interval");
    // ConvgFail=1;  return 0;
                //goto L40;
            }

            /*     Partition the working space array, so that different parts of it can */
            /*     be treated separately during the calculation of BOBYQB. The partition */
            /*     requires the first (NPT+2)*(NPT+N)+3*N*(N+5)/2 elements of W plus the */
            /*     space that is taken by the last array in the argument list of BOBYQB. */

            ndim = npt + n;
            ixb = 1;
            ixp = ixb + n;
            ifv = ixp + n * npt;
            ixo = ifv + npt;
            igo = ixo + n;
            ihq = igo + n;
            ipq = ihq + n * np / 2;
            ibmat = ipq + npt;
            izmat = ibmat + ndim * n;
            isl = izmat + npt * (npt - np);
            isu = isl + n;
            ixn = isu + n;
            ixa = ixn + n;
            id_ = ixa + n;
            ivl = id_ + n;
            iw = ivl + ndim;

            /*     Return if there is insufficient space between the bounds. Modify the */
            /*     initial X if necessary in order to avoid conflicts between the bounds */
            /*     and the construction of the first quadratic model. The lower and upper */
            /*     bounds on moves from the updated X are set now, in the ISL and ISU */
            /*     partitions of W, in order to provide useful and exact information about */
            /*     components of X that become within distance RHOBEG from their bounds. */

            zero = 0.;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                temp = xu[j] - xl[j];
                if (temp < rhobeg + rhobeg) {
                    throw bobyqa_failure("Return from BOBYQA because one of the differences in x_lower and x_upper is less than 2*rho_begin");
                    //goto L40;
                }
                jsl = isl + j - 1;
                jsu = jsl + n;
                w[jsl] = xl[j] - x[j];
                w[jsu] = xu[j] - x[j];
                if (w[jsl] >= -(rhobeg)) {
                    if (w[jsl] >= zero) {
                        x[j] = xl[j];
                        w[jsl] = zero;
                        w[jsu] = temp;
                    } else {
                        x[j] = xl[j] + rhobeg;
                        w[jsl] = -(rhobeg);
                        /* Computing MAX */
                        d__1 = xu[j] - x[j];
                        w[jsu] = std::max(d__1,rhobeg);
                    }
                } else if (w[jsu] <= rhobeg) {
                    if (w[jsu] <= zero) {
                        x[j] = xu[j];
                        w[jsl] = -temp;
                        w[jsu] = zero;
                    } else {
                        x[j] = xu[j] - rhobeg;
                        /* Computing MIN */
                        d__1 = xl[j] - x[j], d__2 = -(rhobeg);
                        w[jsl] = std::min(d__1,d__2);
                        w[jsu] = rhobeg;
                    }
                }
                /* L30: */
            }

            /*     Make the call of BOBYQB. */

            return bobyqb_(calfun, n, npt, &x[1], &xl[1], &xu[1], rhobeg, rhoend, maxfun, &w[
                    ixb], &w[ixp], &w[ifv], &w[ixo], &w[igo], &w[ihq], &w[ipq], &w[
                    ibmat], &w[izmat], ndim, &w[isl], &w[isu], &w[ixn], &w[ixa], &w[
                    id_], &w[ivl], &w[iw],ConvgFail);
            //L40:
            ;
        } /* bobyqa_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        doublereal bobyqb_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const doublereal rhoend,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *xopt,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *sl,
            doublereal *su,
            doublereal *xnew,
            doublereal *xalt,
            doublereal *d__,
            doublereal *vlag,
            doublereal *w,  int& ConvgFail
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2, i__3;
            doublereal d__1, d__2, d__3, d__4;

            /* Local variables */
            doublereal f = 0;
            integer i__, j, k, ih, nf, jj, nh, ip, jp;
            doublereal dx;
            integer np;
            doublereal den = 0, one = 0, ten = 0, dsq = 0, rho = 0, sum = 0, two = 0, diff = 0, half = 0, beta = 0, gisq = 0;
            integer knew = 0;
            doublereal temp, suma, sumb, bsum, fopt;
            integer kopt = 0, nptm;
            doublereal zero, curv;
            integer ksav;
            doublereal gqsq = 0, dist = 0, sumw = 0, sumz = 0, diffa = 0, diffb = 0, diffc = 0, hdiag = 0;
            integer kbase;
            doublereal alpha = 0, delta = 0, adelt = 0, denom = 0, fsave = 0, bdtol = 0, delsq = 0;
            integer nresc, nfsav;
            doublereal ratio = 0, dnorm = 0, vquad = 0, pqold = 0, tenth = 0;
            integer itest;
            doublereal sumpq, scaden;
            doublereal errbig, cauchy, fracsq, biglsq, densav;
            doublereal bdtest;
            doublereal crvmin, frhosq;
            doublereal distsq;
            integer ntrits;
            doublereal xoptsq;



            /*     The arguments N, NPT, X, XL, XU, RHOBEG, RHOEND, IPRINT and MAXFUN */
            /*       are identical to the corresponding arguments in SUBROUTINE BOBYQA. */
            /*     XBASE holds a shift of origin that should reduce the contributions */
            /*       from rounding errors to values of the model and Lagrange functions. */
            /*     XPT is a two-dimensional array that holds the coordinates of the */
            /*       interpolation points relative to XBASE. */
            /*     FVAL holds the values of F at the interpolation points. */
            /*     XOPT is set to the displacement from XBASE of the trust region centre. */
            /*     GOPT holds the gradient of the quadratic model at XBASE+XOPT. */
            /*     HQ holds the explicit second derivatives of the quadratic model. */
            /*     PQ contains the parameters of the implicit second derivatives of the */
            /*       quadratic model. */
            /*     BMAT holds the last N columns of H. */
            /*     ZMAT holds the factorization of the leading NPT by NPT submatrix of H, */
            /*       this factorization being ZMAT times ZMAT^T, which provides both the */
            /*       correct rank and positive semi-definiteness. */
            /*     NDIM is the first dimension of BMAT and has the value NPT+N. */
            /*     SL and SU hold the differences XL-XBASE and XU-XBASE, respectively. */
            /*       All the components of every XOPT are going to satisfy the bounds */
            /*       SL(I) .LEQ. XOPT(I) .LEQ. SU(I), with appropriate equalities when */
            /*       XOPT is on a constraint boundary. */
            /*     XNEW is chosen by SUBROUTINE TRSBOX or ALTMOV. Usually XBASE+XNEW is the */
            /*       vector of variables for the next call of CALFUN. XNEW also satisfies */
            /*       the SL and SU constraints in the way that has just been mentioned. */
            /*     XALT is an alternative to XNEW, chosen by ALTMOV, that may replace XNEW */
            /*       in order to increase the denominator in the updating of UPDATE. */
            /*     D is reserved for a trial step from XOPT, which is usually XNEW-XOPT. */
            /*     VLAG contains the values of the Lagrange functions at a new point X. */
            /*       They are part of a product that requires VLAG to be of length NDIM. */
            /*     W is a one-dimensional array that is used for working space. Its length */
            /*       must be at least 3*NDIM = 3*(NPT+N). */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --x;
            --xl;
            --xu;
            --xbase;
            --fval;
            --xopt;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --xnew;
            --xalt;
            --d__;
            --vlag;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            ten = 10.;
            tenth = .1;
            two = 2.;
            zero = 0.;
            np = n + 1;
            nptm = npt - np;
            nh = n * np / 2;

            /*     The call of PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
            /*     BMAT and ZMAT for the first iteration, with the corresponding values of */
            /*     of NF and KOPT, which are the number of calls of CALFUN so far and the */
            /*     index of the interpolation point at the trust region centre. Then the */
            /*     initial XOPT is set too. The branch to label 720 occurs if MAXFUN is */
            /*     less than NPT. GOPT will be updated if KOPT is different from KBASE. */

            prelim_(calfun, n, npt, &x[1], &xl[1], &xu[1], rhobeg, maxfun, &xbase[1],
                    &xpt[xpt_offset], &fval[1], &gopt[1], &hq[1], &pq[1], &bmat[bmat_offset],
                    &zmat[zmat_offset], ndim, &sl[1], &su[1], nf, kopt);
            xoptsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xopt[i__] = xpt[kopt + i__ * xpt_dim1];
                /* L10: */
                /* Computing 2nd power */
                d__1 = xopt[i__];
                xoptsq += d__1 * d__1;
            }
            fsave = fval[1];
            if (nf < npt) {
               // throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
ConvgFail=1;  return 0;  //  objective function has been called max_f_evals times.
                //goto L720;
            }
            kbase = 1;

            /*     Complete the settings that are required for the iterative procedure. */

            rho = rhobeg;
            delta = rho;
            nresc = nf;
            ntrits = 0;
            diffa = zero;
            diffb = zero;
            itest = 0;
            nfsav = nf;

            /*     Update GOPT if necessary before the first iteration and after each */
            /*     call of RESCUE that makes a call of CALFUN. */

L20:
            if (kopt != kbase) {
                ih = 0;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    i__2 = j;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        ++ih;
                        if (i__ < j) {
                            gopt[j] += hq[ih] * xopt[i__];
                        }
                        /* L30: */
                        gopt[i__] += hq[ih] * xopt[j];
                    }
                }
                if (nf > npt) {
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        temp = zero;
                        i__1 = n;
                        for (j = 1; j <= i__1; ++j) {
                            /* L40: */
                            temp += xpt[k + j * xpt_dim1] * xopt[j];
                        }
                        temp = pq[k] * temp;
                        i__1 = n;
                        for (i__ = 1; i__ <= i__1; ++i__) {
                            /* L50: */
                            gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
                        }
                    }
                }
            }

            /*     Generate the next point in the trust region that provides a small value */
            /*     of the quadratic model subject to the constraints on the variables. */
            /*     The integer NTRITS is set to the number "trust region" iterations that */
            /*     have occurred since the last "alternative" iteration. If the length */
            /*     of XNEW-XOPT is less than HALF*RHO, however, then there is a branch to */
            /*     label 650 or 680 with NTRITS=-1, instead of calculating F at XNEW. */

L60:
            trsbox_(n, npt, &xpt[xpt_offset], &xopt[1], &gopt[1], &hq[1], &pq[1], &sl[1],
                    &su[1], delta, &xnew[1], &d__[1], &w[1], &w[np], &w[np + n],
                    &w[np + (n << 1)], &w[np + n * 3], &dsq, &crvmin);
            /* Computing MIN */
            d__1 = delta, d__2 = std::sqrt(dsq);
            dnorm = std::min(d__1,d__2);
            if (dnorm < half * rho) {
                ntrits = -1;
                /* Computing 2nd power */
                d__1 = ten * rho;
                distsq = d__1 * d__1;
                if (nf <= nfsav + 2) {
                    goto L650;
                }

                /*     The following choice between labels 650 and 680 depends on whether or */
                /*     not our work with the current RHO seems to be complete. Either RHO is */
                /*     decreased or termination occurs if the errors in the quadratic model at */
                /*     the last three interpolation points compare favourably with predictions */
                /*     of likely improvements to the model within distance HALF*RHO of XOPT. */

                /* Computing MAX */
                d__1 = std::max(diffa,diffb);
                errbig = std::max(d__1,diffc);
                frhosq = rho * .125 * rho;
                if (crvmin > zero && errbig > frhosq * crvmin) {
                    goto L650;
                }
                bdtol = errbig / rho;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    bdtest = bdtol;
                    if (xnew[j] == sl[j]) {
                        bdtest = w[j];
                    }
                    if (xnew[j] == su[j]) {
                        bdtest = -w[j];
                    }
                    if (bdtest < bdtol) {
                        curv = hq[(j + j * j) / 2];
                        i__2 = npt;
                        for (k = 1; k <= i__2; ++k) {
                            /* L70: */
                            /* Computing 2nd power */
                            d__1 = xpt[k + j * xpt_dim1];
                            curv += pq[k] * (d__1 * d__1);
                        }
                        bdtest += half * curv * rho;
                        if (bdtest < bdtol) {
                            goto L650;
                        }
                    }
                    /* L80: */
                }
                goto L680;
            }
            ++ntrits;

            /*     Severe cancellation is likely to occur if XOPT is too far from XBASE. */
            /*     If the following test holds, then XBASE is shifted so that XOPT becomes */
            /*     zero. The appropriate changes are made to BMAT and to the second */
            /*     derivatives of the current model, beginning with the changes to BMAT */
            /*     that do not depend on ZMAT. VLAG is used temporarily for working space. */

L90:
            if (dsq <= xoptsq * .001) {
                fracsq = xoptsq * .25;
                sumpq = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    sumpq += pq[k];
                    sum = -half * xoptsq;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        /* L100: */
                        sum += xpt[k + i__ * xpt_dim1] * xopt[i__];
                    }
                    w[npt + k] = sum;
                    temp = fracsq - half * sum;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        w[i__] = bmat[k + i__ * bmat_dim1];
                        vlag[i__] = sum * xpt[k + i__ * xpt_dim1] + temp * xopt[i__];
                        ip = npt + i__;
                        i__3 = i__;
                        for (j = 1; j <= i__3; ++j) {
                            /* L110: */
                            bmat[ip + j * bmat_dim1] = bmat[ip + j * bmat_dim1] + w[
                                i__] * vlag[j] + vlag[i__] * w[j];
                        }
                    }
                }

                /*     Then the revisions of BMAT that depend on ZMAT are calculated. */

                i__3 = nptm;
                for (jj = 1; jj <= i__3; ++jj) {
                    sumz = zero;
                    sumw = zero;
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        sumz += zmat[k + jj * zmat_dim1];
                        vlag[k] = w[npt + k] * zmat[k + jj * zmat_dim1];
                        /* L120: */
                        sumw += vlag[k];
                    }
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        sum = (fracsq * sumz - half * sumw) * xopt[j];
                        i__1 = npt;
                        for (k = 1; k <= i__1; ++k) {
                            /* L130: */
                            sum += vlag[k] * xpt[k + j * xpt_dim1];
                        }
                        w[j] = sum;
                        i__1 = npt;
                        for (k = 1; k <= i__1; ++k) {
                            /* L140: */
                            bmat[k + j * bmat_dim1] += sum * zmat[k + jj * zmat_dim1];
                        }
                    }
                    i__1 = n;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ip = i__ + npt;
                        temp = w[i__];
                        i__2 = i__;
                        for (j = 1; j <= i__2; ++j) {
                            /* L150: */
                            bmat[ip + j * bmat_dim1] += temp * w[j];
                        }
                    }
                }

                /*     The following instructions complete the shift, including the changes */
                /*     to the second derivative parameters of the quadratic model. */

                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    w[j] = -half * sumpq * xopt[j];
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        w[j] += pq[k] * xpt[k + j * xpt_dim1];
                        /* L160: */
                        xpt[k + j * xpt_dim1] -= xopt[j];
                    }
                    i__1 = j;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ++ih;
                        hq[ih] = hq[ih] + w[i__] * xopt[j] + xopt[i__] * w[j];
                        /* L170: */
                        bmat[npt + i__ + j * bmat_dim1] = bmat[npt + j + i__ *
                            bmat_dim1];
                    }
                }
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    xbase[i__] += xopt[i__];
                    xnew[i__] -= xopt[i__];
                    sl[i__] -= xopt[i__];
                    su[i__] -= xopt[i__];
                    /* L180: */
                    xopt[i__] = zero;
                }
                xoptsq = zero;
            }
            if (ntrits == 0) {
                goto L210;
            }
            goto L230;

            /*     XBASE is also moved to XOPT by a call of RESCUE. This calculation is */
            /*     more expensive than the previous shift, because new matrices BMAT and */
            /*     ZMAT are generated from scratch, which may include the replacement of */
            /*     interpolation points whose positions seem to be causing near linear */
            /*     dependence in the interpolation conditions. Therefore RESCUE is called */
            /*     only if rounding errors have reduced by at least a factor of two the */
            /*     denominator of the formula for updating the H matrix. It provides a */
            /*     useful safeguard, but is not invoked in most applications of BOBYQA. */

L190:
            nfsav = nf;
            kbase = kopt;
            rescue_(calfun, n, npt, &xl[1], &xu[1], maxfun, &xbase[1], &xpt[
                    xpt_offset], &fval[1], &xopt[1], &gopt[1], &hq[1], &pq[1], &bmat[
                    bmat_offset], &zmat[zmat_offset], ndim, &sl[1], &su[1], nf, delta,
                    kopt, &vlag[1], &w[1], &w[n + np], &w[ndim + np]);

            /*     XOPT is updated now in case the branch below to label 720 is taken. */
            /*     Any updating of GOPT occurs after the branch below to label 20, which */
            /*     leads to a trust region iteration as does the branch to label 60. */

            xoptsq = zero;
            if (kopt != kbase) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    xopt[i__] = xpt[kopt + i__ * xpt_dim1];
                    /* L200: */
                    /* Computing 2nd power */
                    d__1 = xopt[i__];
                    xoptsq += d__1 * d__1;
                }
            }
            if (nf < 0) {
                nf = maxfun;
                // throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
     ConvgFail=1;  return 0;
                //goto L720;
            }
            nresc = nf;
            if (nfsav < nf) {
                nfsav = nf;
                goto L20;
            }
            if (ntrits > 0) {
                goto L60;
            }

            /*     Pick two alternative vectors of variables, relative to XBASE, that */
            /*     are suitable as new positions of the KNEW-th interpolation point. */
            /*     Firstly, XNEW is set to the point on a line through XOPT and another */
            /*     interpolation point that minimizes the predicted value of the next */
            /*     denominator, subject to ||XNEW - XOPT|| .LEQ. ADELT and to the SL */
            /*     and SU bounds. Secondly, XALT is set to the best feasible point on */
            /*     a constrained version of the Cauchy step of the KNEW-th Lagrange */
            /*     function, the corresponding value of the square of this function */
            /*     being returned in CAUCHY. The choice between these alternatives is */
            /*     going to be made when the denominator is calculated. */

L210:
            altmov_(n, npt, &xpt[xpt_offset], &xopt[1], &bmat[bmat_offset], &zmat[zmat_offset],
                    ndim, &sl[1], &su[1], kopt, knew, adelt, &xnew[1],
                    &xalt[1], alpha, cauchy, &w[1], &w[np], &w[ndim + 1]);
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                /* L220: */
                d__[i__] = xnew[i__] - xopt[i__];
            }

            /*     Calculate VLAG and BETA for the current choice of D. The scalar */
            /*     product of D with XPT(K,.) is going to be held in W(NPT+K) for */
            /*     use when VQUAD is calculated. */

L230:
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                suma = zero;
                sumb = zero;
                sum = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    suma += xpt[k + j * xpt_dim1] * d__[j];
                    sumb += xpt[k + j * xpt_dim1] * xopt[j];
                    /* L240: */
                    sum += bmat[k + j * bmat_dim1] * d__[j];
                }
                w[k] = suma * (half * suma + sumb);
                vlag[k] = sum;
                /* L250: */
                w[npt + k] = suma;
            }
            beta = zero;
            i__1 = nptm;
            for (jj = 1; jj <= i__1; ++jj) {
                sum = zero;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L260: */
                    sum += zmat[k + jj * zmat_dim1] * w[k];
                }
                beta -= sum * sum;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L270: */
                    vlag[k] += sum * zmat[k + jj * zmat_dim1];
                }
            }
            dsq = zero;
            bsum = zero;
            dx = zero;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                /* Computing 2nd power */
                d__1 = d__[j];
                dsq += d__1 * d__1;
                sum = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L280: */
                    sum += w[k] * bmat[k + j * bmat_dim1];
                }
                bsum += sum * d__[j];
                jp = npt + j;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L290: */
                    sum += bmat[jp + i__ * bmat_dim1] * d__[i__];
                }
                vlag[jp] = sum;
                bsum += sum * d__[j];
                /* L300: */
                dx += d__[j] * xopt[j];
            }
            beta = dx * dx + dsq * (xoptsq + dx + dx + half * dsq) + beta - bsum;
            vlag[kopt] += one;

            /*     If NTRITS is zero, the denominator may be increased by replacing */
            /*     the step D of ALTMOV by a Cauchy step. Then RESCUE may be called if */
            /*     rounding errors have damaged the chosen denominator. */

            if (ntrits == 0) {
                /* Computing 2nd power */
                d__1 = vlag[knew];
                denom = d__1 * d__1 + alpha * beta;
                if (denom < cauchy && cauchy > zero) {
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        xnew[i__] = xalt[i__];
                        /* L310: */
                        d__[i__] = xnew[i__] - xopt[i__];
                    }
                    cauchy = zero;
                    goto L230;
                }
                /* Computing 2nd power */
                d__1 = vlag[knew];
                if (denom <= half * (d__1 * d__1)) {
                    if (nf > nresc) {
                        goto L190;
                    }
          //  throw bobyqa_failure("Return from BOBYQA because of much cancellation in a denominator.");
          ConvgFail=1;  return 0;
                    //goto L720;
                }

                /*     Alternatively, if NTRITS is positive, then set KNEW to the index of */
                /*     the next interpolation point to be deleted to make room for a trust */
                /*     region step. Again RESCUE may be called if rounding errors have damaged */
                /*     the chosen denominator, which is the reason for attempting to select */
                /*     KNEW before calculating the next value of the objective function. */

            } else {
                delsq = delta * delta;
                scaden = zero;
                biglsq = zero;
                knew = 0;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    if (k == kopt) {
                        goto L350;
                    }
                    hdiag = zero;
                    i__1 = nptm;
                    for (jj = 1; jj <= i__1; ++jj) {
                        /* L330: */
                        /* Computing 2nd power */
                        d__1 = zmat[k + jj * zmat_dim1];
                        hdiag += d__1 * d__1;
                    }
                    /* Computing 2nd power */
                    d__1 = vlag[k];
                    den = beta * hdiag + d__1 * d__1;
                    distsq = zero;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L340: */
                        /* Computing 2nd power */
                        d__1 = xpt[k + j * xpt_dim1] - xopt[j];
                        distsq += d__1 * d__1;
                    }
                    /* Computing MAX */
                    /* Computing 2nd power */
                    d__3 = distsq / delsq;
                    d__1 = one, d__2 = d__3 * d__3;
                    temp = std::max(d__1,d__2);
                    if (temp * den > scaden) {
                        scaden = temp * den;
                        knew = k;
                        denom = den;
                    }
                    /* Computing MAX */
                    /* Computing 2nd power */
                    d__3 = vlag[k];
                    d__1 = biglsq, d__2 = temp * (d__3 * d__3);
                    biglsq = std::max(d__1,d__2);
L350:
                    ;
                }
                if (scaden <= half * biglsq) {
                    if (nf > nresc) {
                        goto L190;
                    }
  //  throw bobyqa_failure("Return from BOBYQA because of much cancellation in a denominator.");
                    //goto L720;
   ConvgFail=1;  return 0;
                }
            }

            /*     Put the variables for the next calculation of the objective function */
            /*       in XNEW, with any adjustments for the bounds. */


            /*     Calculate the value of the objective function at XBASE+XNEW, unless */
            /*       the limit on the number of calculations of F has been reached. */

L360:
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* Computing MIN */
                /* Computing MAX */
                d__3 = xl[i__], d__4 = xbase[i__] + xnew[i__];
                d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                x[i__] = std::min(d__1,d__2);
                if (xnew[i__] == sl[i__]) {
                    x[i__] = xl[i__];
                }
                if (xnew[i__] == su[i__]) {
                    x[i__] = xu[i__];
                }
                /* L380: */
            }
            if (nf >= maxfun) {
              //  throw bobyqa_failure("Return from BOBYQA because the objective function has been called max_f_evals times.");
	    ConvgFail=1;  return 0;
                //goto L720;
            }
            ++nf;
            f = calfun(mat(&x[1], n));
            if (ntrits == -1) {
                fsave = f;
                goto L720;
            }

            /*     Use the quadratic model to predict the change in F due to the step D, */
            /*       and set DIFF to the error of this prediction. */

            fopt = fval[kopt];
            vquad = zero;
            ih = 0;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                vquad += d__[j] * gopt[j];
                i__1 = j;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    ++ih;
                    temp = d__[i__] * d__[j];
                    if (i__ == j) {
                        temp = half * temp;
                    }
                    /* L410: */
                    vquad += hq[ih] * temp;
                }
            }
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                /* L420: */
                /* Computing 2nd power */
                d__1 = w[npt + k];
                vquad += half * pq[k] * (d__1 * d__1);
            }
            diff = f - fopt - vquad;
            diffc = diffb;
            diffb = diffa;
            diffa = std::abs(diff);
            if (dnorm > rho) {
                nfsav = nf;
            }

            /*     Pick the next value of DELTA after a trust region step. */

            if (ntrits > 0) {
                if (vquad >= zero) {
  //  throw ("Return from BOBYQA because a trust region step has failed to reduce Q.");
                    //goto L720;
 { return 0; ConvgFail=1; }
                }
                ratio = (f - fopt) / vquad;
                if (ratio <= tenth) {
                    /* Computing MIN */
                    d__1 = half * delta;
                    delta = std::min(d__1,dnorm);
                } else if (ratio <= .7) {
                    /* Computing MAX */
                    d__1 = half * delta;
                    delta = std::max(d__1,dnorm);
                } else {
                    /* Computing MAX */
                    d__1 = half * delta, d__2 = dnorm + dnorm;
                    delta = std::max(d__1,d__2);
                }
                if (delta <= rho * 1.5) {
                    delta = rho;
                }

                /*     Recalculate KNEW and DENOM if the new F is less than FOPT. */

                if (f < fopt) {
                    ksav = knew;
                    densav = denom;
                    delsq = delta * delta;
                    scaden = zero;
                    biglsq = zero;
                    knew = 0;
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        hdiag = zero;
                        i__2 = nptm;
                        for (jj = 1; jj <= i__2; ++jj) {
                            /* L440: */
                            /* Computing 2nd power */
                            d__1 = zmat[k + jj * zmat_dim1];
                            hdiag += d__1 * d__1;
                        }
                        /* Computing 2nd power */
                        d__1 = vlag[k];
                        den = beta * hdiag + d__1 * d__1;
                        distsq = zero;
                        i__2 = n;
                        for (j = 1; j <= i__2; ++j) {
                            /* L450: */
                            /* Computing 2nd power */
                            d__1 = xpt[k + j * xpt_dim1] - xnew[j];
                            distsq += d__1 * d__1;
                        }
                        /* Computing MAX */
                        /* Computing 2nd power */
                        d__3 = distsq / delsq;
                        d__1 = one, d__2 = d__3 * d__3;
                        temp = std::max(d__1,d__2);
                        if (temp * den > scaden) {
                            scaden = temp * den;
                            knew = k;
                            denom = den;
                        }
                        /* L460: */
                        /* Computing MAX */
                        /* Computing 2nd power */
                        d__3 = vlag[k];
                        d__1 = biglsq, d__2 = temp * (d__3 * d__3);
                        biglsq = std::max(d__1,d__2);
                    }
                    if (scaden <= half * biglsq) {
                        knew = ksav;
                        denom = densav;
                    }
                }
            }

            /*     Update BMAT and ZMAT, so that the KNEW-th interpolation point can be */
            /*     moved. Also update the second derivative terms of the model. */

            update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1],
                    beta, denom, knew, &w[1]);
            ih = 0;
            pqold = pq[knew];
            pq[knew] = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                temp = pqold * xpt[knew + i__ * xpt_dim1];
                i__2 = i__;
                for (j = 1; j <= i__2; ++j) {
                    ++ih;
                    /* L470: */
                    hq[ih] += temp * xpt[knew + j * xpt_dim1];
                }
            }
            i__2 = nptm;
            for (jj = 1; jj <= i__2; ++jj) {
                temp = diff * zmat[knew + jj * zmat_dim1];
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L480: */
                    pq[k] += temp * zmat[k + jj * zmat_dim1];
                }
            }

            /*     Include the new interpolation point, and make the changes to GOPT at */
            /*     the old XOPT that are caused by the updating of the quadratic model. */

            fval[knew] = f;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xpt[knew + i__ * xpt_dim1] = xnew[i__];
                /* L490: */
                w[i__] = bmat[knew + i__ * bmat_dim1];
            }
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                suma = zero;
                i__2 = nptm;
                for (jj = 1; jj <= i__2; ++jj) {
                    /* L500: */
                    suma += zmat[knew + jj * zmat_dim1] * zmat[k + jj * zmat_dim1];
                }
                sumb = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L510: */
                    sumb += xpt[k + j * xpt_dim1] * xopt[j];
                }
                temp = suma * sumb;
                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L520: */
                    w[i__] += temp * xpt[k + i__ * xpt_dim1];
                }
            }
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L530: */
                gopt[i__] += diff * w[i__];
            }

            /*     Update XOPT, GOPT and KOPT if the new calculated F is less than FOPT. */

            if (f < fopt) {
                kopt = knew;
                xoptsq = zero;
                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    xopt[j] = xnew[j];
                    /* Computing 2nd power */
                    d__1 = xopt[j];
                    xoptsq += d__1 * d__1;
                    i__1 = j;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        ++ih;
                        if (i__ < j) {
                            gopt[j] += hq[ih] * d__[i__];
                        }
                        /* L540: */
                        gopt[i__] += hq[ih] * d__[j];
                    }
                }
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    temp = zero;
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        /* L550: */
                        temp += xpt[k + j * xpt_dim1] * d__[j];
                    }
                    temp = pq[k] * temp;
                    i__2 = n;
                    for (i__ = 1; i__ <= i__2; ++i__) {
                        /* L560: */
                        gopt[i__] += temp * xpt[k + i__ * xpt_dim1];
                    }
                }
            }

            /*     Calculate the parameters of the least Frobenius norm interpolant to */
            /*     the current data, the gradient of this interpolant at XOPT being put */
            /*     into VLAG(NPT+I), I=1,2,...,N. */

            if (ntrits > 0) {
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    vlag[k] = fval[k] - fval[kopt];
                    /* L570: */
                    w[k] = zero;
                }
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    sum = zero;
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        /* L580: */
                        sum += zmat[k + j * zmat_dim1] * vlag[k];
                    }
                    i__1 = npt;
                    for (k = 1; k <= i__1; ++k) {
                        /* L590: */
                        w[k] += sum * zmat[k + j * zmat_dim1];
                    }
                }
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    sum = zero;
                    i__2 = n;
                    for (j = 1; j <= i__2; ++j) {
                        /* L600: */
                        sum += xpt[k + j * xpt_dim1] * xopt[j];
                    }
                    w[k + npt] = w[k];
                    /* L610: */
                    w[k] = sum * w[k];
                }
                gqsq = zero;
                gisq = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    sum = zero;
                    i__2 = npt;
                    for (k = 1; k <= i__2; ++k) {
                        /* L620: */
                        sum = sum + bmat[k + i__ * bmat_dim1] * vlag[k] + xpt[k + i__
                            * xpt_dim1] * w[k];
                    }
                    if (xopt[i__] == sl[i__]) {
                        /* Computing MIN */
                        d__2 = zero, d__3 = gopt[i__];
                        /* Computing 2nd power */
                        d__1 = std::min(d__2,d__3);
                        gqsq += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = std::min(zero,sum);
                        gisq += d__1 * d__1;
                    } else if (xopt[i__] == su[i__]) {
                        /* Computing MAX */
                        d__2 = zero, d__3 = gopt[i__];
                        /* Computing 2nd power */
                        d__1 = std::max(d__2,d__3);
                        gqsq += d__1 * d__1;
                        /* Computing 2nd power */
                        d__1 = std::max(zero,sum);
                        gisq += d__1 * d__1;
                    } else {
                        /* Computing 2nd power */
                        d__1 = gopt[i__];
                        gqsq += d__1 * d__1;
                        gisq += sum * sum;
                    }
                    /* L630: */
                    vlag[npt + i__] = sum;
                }

                /*     Test whether to replace the new quadratic model by the least Frobenius */
                /*     norm interpolant, making the replacement if the test is satisfied. */

                ++itest;
                if (gqsq < ten * gisq) {
                    itest = 0;
                }
                if (itest >= 3) {
                    i__1 = std::max(npt,nh);
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        if (i__ <= n) {
                            gopt[i__] = vlag[npt + i__];
                        }
                        if (i__ <= npt) {
                            pq[i__] = w[npt + i__];
                        }
                        if (i__ <= nh) {
                            hq[i__] = zero;
                        }
                        itest = 0;
                        /* L640: */
                    }
                }
            }

            /*     If a trust region step has provided a sufficient decrease in F, then */
            /*     branch for another trust region calculation. The case NTRITS=0 occurs */
            /*     when the new interpolation point was reached by an alternative step. */

            if (ntrits == 0) {
                goto L60;
            }
            if (f <= fopt + tenth * vquad) {
                goto L60;
            }

            /*     Alternatively, find out if the interpolation points are close enough */
            /*       to the best point so far. */

            /* Computing MAX */
            /* Computing 2nd power */
            d__3 = two * delta;
            /* Computing 2nd power */
            d__4 = ten * rho;
            d__1 = d__3 * d__3, d__2 = d__4 * d__4;
            distsq = std::max(d__1,d__2);
L650:
            knew = 0;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                sum = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L660: */
                    /* Computing 2nd power */
                    d__1 = xpt[k + j * xpt_dim1] - xopt[j];
                    sum += d__1 * d__1;
                }
                if (sum > distsq) {
                    knew = k;
                    distsq = sum;
                }
                /* L670: */
            }

            /*     If KNEW is positive, then ALTMOV finds alternative new positions for */
            /*     the KNEW-th interpolation point within distance ADELT of XOPT. It is */
            /*     reached via label 90. Otherwise, there is a branch to label 60 for */
            /*     another trust region iteration, unless the calculations with the */
            /*     current RHO are complete. */

            if (knew > 0) {
                dist = std::sqrt(distsq);
                if (ntrits == -1) {
                    /* Computing MIN */
                    d__1 = tenth * delta, d__2 = half * dist;
                    delta = std::min(d__1,d__2);
                    if (delta <= rho * 1.5) {
                        delta = rho;
                    }
                }
                ntrits = 0;
                /* Computing MAX */
                /* Computing MIN */
                d__2 = tenth * dist;
                d__1 = std::min(d__2,delta);
                adelt = std::max(d__1,rho);
                dsq = adelt * adelt;
                goto L90;
            }
            if (ntrits == -1) {
                goto L680;
            }
            if (ratio > zero) {
                goto L60;
            }
            if (std::max(delta,dnorm) > rho) {
                goto L60;
            }

            /*     The calculations with the current value of RHO are complete. Pick the */
            /*       next values of RHO and DELTA. */

L680:
            if (rho > rhoend) {
                delta = half * rho;
                ratio = rho / rhoend;
                if (ratio <= 16.) {
                    rho = rhoend;
                } else if (ratio <= 250.) {
                    rho = std::sqrt(ratio) * rhoend;
                } else {
                    rho = tenth * rho;
                }
                delta = std::max(delta,rho);
                ntrits = 0;
                nfsav = nf;
                goto L60;
            }

            /*     Return from the calculation, after another Newton-Raphson step, if */
            /*       it is too short to have been tried before. */

            if (ntrits == -1) {
                goto L360;
            }
L720:
            if (fval[kopt] <= fsave) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* Computing MIN */
                    /* Computing MAX */
                    d__3 = xl[i__], d__4 = xbase[i__] + xopt[i__];
                    d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                    x[i__] = std::min(d__1,d__2);
                    if (xopt[i__] == sl[i__]) {
                        x[i__] = xl[i__];
                    }
                    if (xopt[i__] == su[i__]) {
                        x[i__] = xu[i__];
                    }
                    /* L730: */
                }
                f = fval[kopt];
            }

            return f;
        } /* bobyqb_ */

    // ----------------------------------------------------------------------------------------

        void altmov_(
            const integer n,
            const integer npt,
            const doublereal *xpt,
            const doublereal *xopt,
            const doublereal *bmat,
            const doublereal *zmat,
            const integer ndim,
            const doublereal *sl,
            const doublereal *su,
            const integer kopt,
            const integer knew,
            const doublereal adelt,
            doublereal *xnew,
            doublereal *xalt,
            doublereal& alpha,
            doublereal& cauchy,
            doublereal *glag,
            doublereal *hcol,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            integer i__, j, k;
            doublereal ha, gw, one, diff, half;
            integer ilbd, isbd;
            doublereal slbd;
            integer iubd;
            doublereal vlag, subd, temp;
            integer ksav = 0;
            doublereal step = 0, zero = 0, curv = 0;
            integer iflag;
            doublereal scale = 0, csave = 0, tempa = 0, tempb = 0, tempd = 0, const__ = 0, sumin = 0,
                       ggfree = 0;
            integer ibdsav = 0;
            doublereal dderiv = 0, bigstp = 0, predsq = 0, presav = 0, distsq = 0, stpsav = 0, wfixsq = 0, wsqsav = 0;


            /*     The arguments N, NPT, XPT, XOPT, BMAT, ZMAT, NDIM, SL and SU all have */
            /*       the same meanings as the corresponding arguments of BOBYQB. */
            /*     KOPT is the index of the optimal interpolation point. */
            /*     KNEW is the index of the interpolation point that is going to be moved. */
            /*     ADELT is the current trust region bound. */
            /*     XNEW will be set to a suitable new position for the interpolation point */
            /*       XPT(KNEW,.). Specifically, it satisfies the SL, SU and trust region */
            /*       bounds and it should provide a large denominator in the next call of */
            /*       UPDATE. The step XNEW-XOPT from XOPT is restricted to moves along the */
            /*       straight lines through XOPT and another interpolation point. */
            /*     XALT also provides a large value of the modulus of the KNEW-th Lagrange */
            /*       function subject to the constraints that have been mentioned, its main */
            /*       difference from XNEW being that XALT-XOPT is a constrained version of */
            /*       the Cauchy step within the trust region. An exception is that XALT is */
            /*       not calculated if all components of GLAG (see below) are zero. */
            /*     ALPHA will be set to the KNEW-th diagonal element of the H matrix. */
            /*     CAUCHY will be set to the square of the KNEW-th Lagrange function at */
            /*       the step XALT-XOPT from XOPT for the vector XALT that is returned, */
            /*       except that CAUCHY is set to zero if XALT is not calculated. */
            /*     GLAG is a working space vector of length N for the gradient of the */
            /*       KNEW-th Lagrange function at XOPT. */
            /*     HCOL is a working space vector of length NPT for the second derivative */
            /*       coefficients of the KNEW-th Lagrange function. */
            /*     W is a working space vector of length 2N that is going to hold the */
            /*       constrained Cauchy step from XOPT of the Lagrange function, followed */
            /*       by the downhill version of XALT when the uphill step is calculated. */

            /*     Set the first NPT components of W to the leading elements of the */
            /*     KNEW-th column of the H matrix. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xopt;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --xnew;
            --xalt;
            --glag;
            --hcol;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            zero = 0.;
            const__ = one + std::sqrt(2.);
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                /* L10: */
                hcol[k] = zero;
            }
            i__1 = npt - n - 1;
            for (j = 1; j <= i__1; ++j) {
                temp = zmat[knew + j * zmat_dim1];
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L20: */
                    hcol[k] += temp * zmat[k + j * zmat_dim1];
                }
            }
            alpha = hcol[knew];
            ha = half * alpha;

            /*     Calculate the gradient of the KNEW-th Lagrange function at XOPT. */

            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L30: */
                glag[i__] = bmat[knew + i__ * bmat_dim1];
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                temp = zero;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    /* L40: */
                    temp += xpt[k + j * xpt_dim1] * xopt[j];
                }
                temp = hcol[k] * temp;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L50: */
                    glag[i__] += temp * xpt[k + i__ * xpt_dim1];
                }
            }

            /*     Search for a large denominator along the straight lines through XOPT */
            /*     and another interpolation point. SLBD and SUBD will be lower and upper */
            /*     bounds on the step along each of these lines in turn. PREDSQ will be */
            /*     set to the square of the predicted denominator for each line. PRESAV */
            /*     will be set to the largest admissible value of PREDSQ that occurs. */

            presav = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                if (k == kopt) {
                    goto L80;
                }
                dderiv = zero;
                distsq = zero;
                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = xpt[k + i__ * xpt_dim1] - xopt[i__];
                    dderiv += glag[i__] * temp;
                    /* L60: */
                    distsq += temp * temp;
                }
                subd = adelt / std::sqrt(distsq);
                slbd = -subd;
                ilbd = 0;
                iubd = 0;
                sumin = std::min(one,subd);

                /*     Revise SLBD and SUBD if necessary because of the bounds in SL and SU. */

                i__2 = n;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    temp = xpt[k + i__ * xpt_dim1] - xopt[i__];
                    if (temp > zero) {
                        if (slbd * temp < sl[i__] - xopt[i__]) {
                            slbd = (sl[i__] - xopt[i__]) / temp;
                            ilbd = -i__;
                        }
                        if (subd * temp > su[i__] - xopt[i__]) {
                            /* Computing MAX */
                            d__1 = sumin, d__2 = (su[i__] - xopt[i__]) / temp;
                            subd = std::max(d__1,d__2);
                            iubd = i__;
                        }
                    } else if (temp < zero) {
                        if (slbd * temp > su[i__] - xopt[i__]) {
                            slbd = (su[i__] - xopt[i__]) / temp;
                            ilbd = i__;
                        }
                        if (subd * temp < sl[i__] - xopt[i__]) {
                            /* Computing MAX */
                            d__1 = sumin, d__2 = (sl[i__] - xopt[i__]) / temp;
                            subd = std::max(d__1,d__2);
                            iubd = -i__;
                        }
                    }
                    /* L70: */
                }

                /*     Seek a large modulus of the KNEW-th Lagrange function when the index */
                /*     of the other interpolation point on the line through XOPT is KNEW. */

                if (k == knew) {
                    diff = dderiv - one;
                    step = slbd;
                    vlag = slbd * (dderiv - slbd * diff);
                    isbd = ilbd;
                    temp = subd * (dderiv - subd * diff);
                    if (std::abs(temp) > std::abs(vlag)) {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    tempd = half * dderiv;
                    tempa = tempd - diff * slbd;
                    tempb = tempd - diff * subd;
                    if (tempa * tempb < zero) {
                        temp = tempd * tempd / diff;
                        if (std::abs(temp) > std::abs(vlag)) {
                            step = tempd / diff;
                            vlag = temp;
                            isbd = 0;
                        }
                    }

                    /*     Search along each of the other lines through XOPT and another point. */

                } else {
                    step = slbd;
                    vlag = slbd * (one - slbd);
                    isbd = ilbd;
                    temp = subd * (one - subd);
                    if (std::abs(temp) > std::abs(vlag)) {
                        step = subd;
                        vlag = temp;
                        isbd = iubd;
                    }
                    if (subd > half) {
                        if (std::abs(vlag) < .25) {
                            step = half;
                            vlag = .25;
                            isbd = 0;
                        }
                    }
                    vlag *= dderiv;
                }

                /*     Calculate PREDSQ for the current line search and maintain PRESAV. */

                temp = step * (one - step) * distsq;
                predsq = vlag * vlag * (vlag * vlag + ha * temp * temp);
                if (predsq > presav) {
                    presav = predsq;
                    ksav = k;
                    stpsav = step;
                    ibdsav = isbd;
                }
L80:
                ;
            }

            /*     Construct XNEW in a way that satisfies the bound constraints exactly. */

            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                temp = xopt[i__] + stpsav * (xpt[ksav + i__ * xpt_dim1] - xopt[i__]);
                /* L90: */
                /* Computing MAX */
                /* Computing MIN */
                d__3 = su[i__];
                d__1 = sl[i__], d__2 = std::min(d__3,temp);
                xnew[i__] = std::max(d__1,d__2);
            }
            if (ibdsav < 0) {
                xnew[-ibdsav] = sl[-ibdsav];
            }
            if (ibdsav > 0) {
                xnew[ibdsav] = su[ibdsav];
            }

            /*     Prepare for the iterative method that assembles the constrained Cauchy */
            /*     step in W. The sum of squares of the fixed components of W is formed in */
            /*     WFIXSQ, and the free components of W are set to BIGSTP. */

            bigstp = adelt + adelt;
            iflag = 0;
L100:
            wfixsq = zero;
            ggfree = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                w[i__] = zero;
                /* Computing MIN */
                d__1 = xopt[i__] - sl[i__], d__2 = glag[i__];
                tempa = std::min(d__1,d__2);
                /* Computing MAX */
                d__1 = xopt[i__] - su[i__], d__2 = glag[i__];
                tempb = std::max(d__1,d__2);
                if (tempa > zero || tempb < zero) {
                    w[i__] = bigstp;
                    /* Computing 2nd power */
                    d__1 = glag[i__];
                    ggfree += d__1 * d__1;
                }
                /* L110: */
            }
            if (ggfree == zero) {
                cauchy = zero;
                goto L200;
            }

            /*     Investigate whether more components of W can be fixed. */

L120:
            temp = adelt * adelt - wfixsq;
            if (temp > zero) {
                wsqsav = wfixsq;
                step = std::sqrt(temp / ggfree);
                ggfree = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    if (w[i__] == bigstp) {
                        temp = xopt[i__] - step * glag[i__];
                        if (temp <= sl[i__]) {
                            w[i__] = sl[i__] - xopt[i__];
                            /* Computing 2nd power */
                            d__1 = w[i__];
                            wfixsq += d__1 * d__1;
                        } else if (temp >= su[i__]) {
                            w[i__] = su[i__] - xopt[i__];
                            /* Computing 2nd power */
                            d__1 = w[i__];
                            wfixsq += d__1 * d__1;
                        } else {
                            /* Computing 2nd power */
                            d__1 = glag[i__];
                            ggfree += d__1 * d__1;
                        }
                    }
                    /* L130: */
                }
                if (wfixsq > wsqsav && ggfree > zero) {
                    goto L120;
                }
            }

            /*     Set the remaining free components of W and all components of XALT, */
            /*     except that W may be scaled later. */

            gw = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (w[i__] == bigstp) {
                    w[i__] = -step * glag[i__];
                    /* Computing MAX */
                    /* Computing MIN */
                    d__3 = su[i__], d__4 = xopt[i__] + w[i__];
                    d__1 = sl[i__], d__2 = std::min(d__3,d__4);
                    xalt[i__] = std::max(d__1,d__2);
                } else if (w[i__] == zero) {
                    xalt[i__] = xopt[i__];
                } else if (glag[i__] > zero) {
                    xalt[i__] = sl[i__];
                } else {
                    xalt[i__] = su[i__];
                }
                /* L140: */
                gw += glag[i__] * w[i__];
            }

            /*     Set CURV to the curvature of the KNEW-th Lagrange function along W. */
            /*     Scale W by a factor less than one if that can reduce the modulus of */
            /*     the Lagrange function at XOPT+W. Set CAUCHY to the final value of */
            /*     the square of this function. */

            curv = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                temp = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    /* L150: */
                    temp += xpt[k + j * xpt_dim1] * w[j];
                }
                /* L160: */
                curv += hcol[k] * temp * temp;
            }
            if (iflag == 1) {
                curv = -curv;
            }
            if (curv > -gw && curv < -const__ * gw) {
                scale = -gw / curv;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    temp = xopt[i__] + scale * w[i__];
                    /* L170: */
                    /* Computing MAX */
                    /* Computing MIN */
                    d__3 = su[i__];
                    d__1 = sl[i__], d__2 = std::min(d__3,temp);
                    xalt[i__] = std::max(d__1,d__2);
                }
                /* Computing 2nd power */
                d__1 = half * gw * scale;
                cauchy = d__1 * d__1;
            } else {
                /* Computing 2nd power */
                d__1 = gw + half * curv;
                cauchy = d__1 * d__1;
            }

            /*     If IFLAG is zero, then XALT is calculated as before after reversing */
            /*     the sign of GLAG. Thus two XALT vectors become available. The one that */
            /*     is chosen is the one that gives the larger value of CAUCHY. */

            if (iflag == 0) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    glag[i__] = -glag[i__];
                    /* L180: */
                    w[n + i__] = xalt[i__];
                }
                csave = cauchy;
                iflag = 1;
                goto L100;
            }
            if (csave > cauchy) {
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    /* L190: */
                    xalt[i__] = w[n + i__];
                }
                cauchy = csave;
            }
L200:
            ;
        } /* altmov_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        void prelim_(
            const funct& calfun,
            const integer n,
            const integer npt,
            doublereal *x,
            const doublereal *xl,
            const doublereal *xu,
            const doublereal rhobeg,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            const doublereal *sl,
            const doublereal *su,
            integer& nf,
            integer& kopt
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            doublereal f;
            integer i__, j, k, ih, np, nfm;
            doublereal one;
            integer nfx = 0, ipt = 0, jpt = 0;
            doublereal two = 0, fbeg = 0, diff = 0, half = 0, temp = 0, zero = 0, recip = 0, stepa = 0, stepb = 0;
            integer itemp;
            doublereal rhosq;



            /*     The arguments N, NPT, X, XL, XU, RHOBEG, IPRINT and MAXFUN are the */
            /*       same as the corresponding arguments in SUBROUTINE BOBYQA. */
            /*     The arguments XBASE, XPT, FVAL, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU */
            /*       are the same as the corresponding arguments in BOBYQB, the elements */
            /*       of SL and SU being set in BOBYQA. */
            /*     GOPT is usually the gradient of the quadratic model at XOPT+XBASE, but */
            /*       it is set by PRELIM to the gradient of the quadratic model at XBASE. */
            /*       If XOPT is nonzero, BOBYQB will change it to its usual value later. */
            /*     NF is maintaned as the number of calls of CALFUN so far. */
            /*     KOPT will be such that the least calculated value of F so far is at */
            /*       the point XPT(KOPT,.)+XBASE in the space of the variables. */

            /*     SUBROUTINE PRELIM sets the elements of XBASE, XPT, FVAL, GOPT, HQ, PQ, */
            /*     BMAT and ZMAT for the first iteration, and it maintains the values of */
            /*     NF and KOPT. The vector X is also changed by PRELIM. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --x;
            --xl;
            --xu;
            --xbase;
            --fval;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;

            /* Function Body */
            half = .5;
            one = 1.;
            two = 2.;
            zero = 0.;
            rhosq = rhobeg * rhobeg;
            recip = one / rhosq;
            np = n + 1;

            /*     Set XBASE to the initial vector of variables, and set the initial */
            /*     elements of XPT, BMAT, HQ, PQ and ZMAT to zero. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                xbase[j] = x[j];
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L10: */
                    xpt[k + j * xpt_dim1] = zero;
                }
                i__2 = ndim;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L20: */
                    bmat[i__ + j * bmat_dim1] = zero;
                }
            }
            i__2 = n * np / 2;
            for (ih = 1; ih <= i__2; ++ih) {
                /* L30: */
                hq[ih] = zero;
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                pq[k] = zero;
                i__1 = npt - np;
                for (j = 1; j <= i__1; ++j) {
                    /* L40: */
                    zmat[k + j * zmat_dim1] = zero;
                }
            }

            /*     Begin the initialization procedure. NF becomes one more than the number */
            /*     of function values so far. The coordinates of the displacement of the */
            /*     next initial interpolation point from XBASE are set in XPT(NF+1,.). */

            nf = 0;
L50:
            nfm = nf;
            nfx = nf - n;
            ++(nf);
            if (nfm <= n << 1) {
                if (nfm >= 1 && nfm <= n) {
                    stepa = rhobeg;
                    if (su[nfm] == zero) {
                        stepa = -stepa;
                    }
                    xpt[nf + nfm * xpt_dim1] = stepa;
                } else if (nfm > n) {
                    stepa = xpt[nf - n + nfx * xpt_dim1];
                    stepb = -(rhobeg);
                    if (sl[nfx] == zero) {
                        /* Computing MIN */
                        d__1 = two * rhobeg, d__2 = su[nfx];
                        stepb = std::min(d__1,d__2);
                    }
                    if (su[nfx] == zero) {
                        /* Computing MAX */
                        d__1 = -two * rhobeg, d__2 = sl[nfx];
                        stepb = std::max(d__1,d__2);
                    }
                    xpt[nf + nfx * xpt_dim1] = stepb;
                }
            } else {
                itemp = (nfm - np) / n;
                jpt = nfm - itemp * n - n;
                ipt = jpt + itemp;
                if (ipt > n) {
                    itemp = jpt;
                    jpt = ipt - n;
                    ipt = itemp;
                }
                xpt[nf + ipt * xpt_dim1] = xpt[ipt + 1 + ipt * xpt_dim1];
                xpt[nf + jpt * xpt_dim1] = xpt[jpt + 1 + jpt * xpt_dim1];
            }

            /*     Calculate the next value of F. The least function value so far and */
            /*     its index are required. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                /* Computing MIN */
                /* Computing MAX */
                d__3 = xl[j], d__4 = xbase[j] + xpt[nf + j * xpt_dim1];
                d__1 = std::max(d__3,d__4), d__2 = xu[j];
                x[j] = std::min(d__1,d__2);
                if (xpt[nf + j * xpt_dim1] == sl[j]) {
                    x[j] = xl[j];
                }
                if (xpt[nf + j * xpt_dim1] == su[j]) {
                    x[j] = xu[j];
                }
                /* L60: */
            }
            f = calfun(mat(&x[1],n));
            fval[nf] = f;
            if (nf == 1) {
                fbeg = f;
                kopt = 1;
            } else if (f < fval[kopt]) {
                kopt = nf;
            }

            /*     Set the nonzero initial elements of BMAT and the quadratic model in the */
            /*     cases when NF is at most 2*N+1. If NF exceeds N+1, then the positions */
            /*     of the NF-th and (NF-N)-th interpolation points may be switched, in */
            /*     order that the function value at the first of them contributes to the */
            /*     off-diagonal second derivative terms of the initial quadratic model. */

            if (nf <= (n << 1) + 1) {
                if (nf >= 2 && nf <= n + 1) {
                    gopt[nfm] = (f - fbeg) / stepa;
                    if (npt < nf + n) {
                        bmat[nfm * bmat_dim1 + 1] = -one / stepa;
                        bmat[nf + nfm * bmat_dim1] = one / stepa;
                        bmat[npt + nfm + nfm * bmat_dim1] = -half * rhosq;
                    }
                } else if (nf >= n + 2) {
                    ih = nfx * (nfx + 1) / 2;
                    temp = (f - fbeg) / stepb;
                    diff = stepb - stepa;
                    hq[ih] = two * (temp - gopt[nfx]) / diff;
                    gopt[nfx] = (gopt[nfx] * stepb - temp * stepa) / diff;
                    if (stepa * stepb < zero) {
                        if (f < fval[nf - n]) {
                            fval[nf] = fval[nf - n];
                            fval[nf - n] = f;
                            if (kopt == nf) {
                                kopt = nf - n;
                            }
                            xpt[nf - n + nfx * xpt_dim1] = stepb;
                            xpt[nf + nfx * xpt_dim1] = stepa;
                        }
                    }
                    bmat[nfx * bmat_dim1 + 1] = -(stepa + stepb) / (stepa * stepb);
                    bmat[nf + nfx * bmat_dim1] = -half / xpt[nf - n + nfx *
                        xpt_dim1];
                    bmat[nf - n + nfx * bmat_dim1] = -bmat[nfx * bmat_dim1 + 1] -
                        bmat[nf + nfx * bmat_dim1];
                    zmat[nfx * zmat_dim1 + 1] = std::sqrt(two) / (stepa * stepb);
                    zmat[nf + nfx * zmat_dim1] = std::sqrt(half) / rhosq;
                    zmat[nf - n + nfx * zmat_dim1] = -zmat[nfx * zmat_dim1 + 1] -
                        zmat[nf + nfx * zmat_dim1];
                }

                /*     Set the off-diagonal second derivatives of the Lagrange functions and */
                /*     the initial quadratic model. */

            } else {
                ih = ipt * (ipt - 1) / 2 + jpt;
                zmat[nfx * zmat_dim1 + 1] = recip;
                zmat[nf + nfx * zmat_dim1] = recip;
                zmat[ipt + 1 + nfx * zmat_dim1] = -recip;
                zmat[jpt + 1 + nfx * zmat_dim1] = -recip;
                temp = xpt[nf + ipt * xpt_dim1] * xpt[nf + jpt * xpt_dim1];
                hq[ih] = (fbeg - fval[ipt + 1] - fval[jpt + 1] + f) / temp;
            }
            if (nf < npt && nf < maxfun) {
                goto L50;
            }

        } /* prelim_ */

    // ----------------------------------------------------------------------------------------

        template <typename funct>
        void rescue_ (
            const funct& calfun,
            const integer n,
            const integer npt,
            const doublereal *xl,
            const doublereal *xu,
            const integer maxfun,
            doublereal *xbase,
            doublereal *xpt,
            doublereal *fval,
            doublereal *xopt,
            doublereal *gopt,
            doublereal *hq,
            doublereal *pq,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *sl,
            doublereal *su,
            integer& nf,
            const doublereal delta,
            integer& kopt,
            doublereal *vlag,
            doublereal * ptsaux,
            doublereal *ptsid,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, bmat_dim1, bmat_offset, zmat_dim1,
            zmat_offset, i__1, i__2, i__3;
            doublereal d__1, d__2, d__3, d__4;


            /* Local variables */
            doublereal f;
            integer i__, j, k, ih, jp, ip, iq, np, iw;
            doublereal xp = 0, xq = 0, den = 0;
            integer ihp = 0;
            doublereal one;
            integer ihq, jpn, kpt;
            doublereal sum = 0, diff = 0, half = 0, beta = 0;
            integer kold;
            doublereal winc;
            integer nrem, knew;
            doublereal temp, bsum;
            integer nptm;
            doublereal zero = 0, hdiag = 0, fbase = 0, sfrac = 0, denom = 0, vquad = 0, sumpq = 0;
            doublereal dsqmin, distsq, vlmxsq;



            /*     The arguments N, NPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL, XOPT, */
            /*       GOPT, HQ, PQ, BMAT, ZMAT, NDIM, SL and SU have the same meanings as */
            /*       the corresponding arguments of BOBYQB on the entry to RESCUE. */
            /*     NF is maintained as the number of calls of CALFUN so far, except that */
            /*       NF is set to -1 if the value of MAXFUN prevents further progress. */
            /*     KOPT is maintained so that FVAL(KOPT) is the least calculated function */
            /*       value. Its correct value must be given on entry. It is updated if a */
            /*       new least function value is found, but the corresponding changes to */
            /*       XOPT and GOPT have to be made later by the calling program. */
            /*     DELTA is the current trust region radius. */
            /*     VLAG is a working space vector that will be used for the values of the */
            /*       provisional Lagrange functions at each of the interpolation points. */
            /*       They are part of a product that requires VLAG to be of length NDIM. */
            /*     PTSAUX is also a working space array. For J=1,2,...,N, PTSAUX(1,J) and */
            /*       PTSAUX(2,J) specify the two positions of provisional interpolation */
            /*       points when a nonzero step is taken along e_J (the J-th coordinate */
            /*       direction) through XBASE+XOPT, as specified below. Usually these */
            /*       steps have length DELTA, but other lengths are chosen if necessary */
            /*       in order to satisfy the given bounds on the variables. */
            /*     PTSID is also a working space array. It has NPT components that denote */
            /*       provisional new positions of the original interpolation points, in */
            /*       case changes are needed to restore the linear independence of the */
            /*       interpolation conditions. The K-th point is a candidate for change */
            /*       if and only if PTSID(K) is nonzero. In this case let p and q be the */
            /*       integer parts of PTSID(K) and (PTSID(K)-p) multiplied by N+1. If p */
            /*       and q are both positive, the step from XBASE+XOPT to the new K-th */
            /*       interpolation point is PTSAUX(1,p)*e_p + PTSAUX(1,q)*e_q. Otherwise */
            /*       the step is PTSAUX(1,p)*e_p or PTSAUX(2,q)*e_q in the cases q=0 or */
            /*       p=0, respectively. */
            /*     The first NDIM+NPT elements of the array W are used for working space. */
            /*     The final elements of BMAT and ZMAT are set in a well-conditioned way */
            /*       to the values that are appropriate for the new interpolation points. */
            /*     The elements of GOPT, HQ and PQ are also revised to the values that are */
            /*       appropriate to the final quadratic model. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xl;
            --xu;
            --xbase;
            --fval;
            --xopt;
            --gopt;
            --hq;
            --pq;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --sl;
            --su;
            --vlag;
            ptsaux -= 3;
            --ptsid;
            --w;

            /* Function Body */
            half = .5;
            one = 1.;
            zero = 0.;
            np = n + 1;
            sfrac = half / (doublereal) np;
            nptm = npt - np;

            /*     Shift the interpolation points so that XOPT becomes the origin, and set */
            /*     the elements of ZMAT to zero. The value of SUMPQ is required in the */
            /*     updating of HQ below. The squares of the distances from XOPT to the */
            /*     other interpolation points are set at the end of W. Increments of WINC */
            /*     may be added later to these squares to balance the consideration of */
            /*     the choice of point that is going to become current. */

            sumpq = zero;
            winc = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                distsq = zero;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    xpt[k + j * xpt_dim1] -= xopt[j];
                    /* L10: */
                    /* Computing 2nd power */
                    d__1 = xpt[k + j * xpt_dim1];
                    distsq += d__1 * d__1;
                }
                sumpq += pq[k];
                w[ndim + k] = distsq;
                winc = std::max(winc,distsq);
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    /* L20: */
                    zmat[k + j * zmat_dim1] = zero;
                }
            }

            /*     Update HQ so that HQ and PQ define the second derivatives of the model */
            /*     after XBASE has been shifted to the trust region centre. */

            ih = 0;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                w[j] = half * sumpq * xopt[j];
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L30: */
                    w[j] += pq[k] * xpt[k + j * xpt_dim1];
                }
                i__1 = j;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    ++ih;
                    /* L40: */
                    hq[ih] = hq[ih] + w[i__] * xopt[j] + w[j] * xopt[i__];
                }
            }

            /*     Shift XBASE, SL, SU and XOPT. Set the elements of BMAT to zero, and */
            /*     also set the elements of PTSAUX. */

            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                xbase[j] += xopt[j];
                sl[j] -= xopt[j];
                su[j] -= xopt[j];
                xopt[j] = zero;
                /* Computing MIN */
                d__1 = delta, d__2 = su[j];
                ptsaux[(j << 1) + 1] = std::min(d__1,d__2);
                /* Computing MAX */
                d__1 = -(delta), d__2 = sl[j];
                ptsaux[(j << 1) + 2] = std::max(d__1,d__2);
                if (ptsaux[(j << 1) + 1] + ptsaux[(j << 1) + 2] < zero) {
                    temp = ptsaux[(j << 1) + 1];
                    ptsaux[(j << 1) + 1] = ptsaux[(j << 1) + 2];
                    ptsaux[(j << 1) + 2] = temp;
                }
                if ((d__2 = ptsaux[(j << 1) + 2], std::abs(d__2)) < half * (d__1 = ptsaux[(
                            j << 1) + 1], std::abs(d__1))) {
                    ptsaux[(j << 1) + 2] = half * ptsaux[(j << 1) + 1];
                }
                i__2 = ndim;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    /* L50: */
                    bmat[i__ + j * bmat_dim1] = zero;
                }
            }
            fbase = fval[kopt];

            /*     Set the identifiers of the artificial interpolation points that are */
            /*     along a coordinate direction from XOPT, and set the corresponding */
            /*     nonzero elements of BMAT and ZMAT. */

            ptsid[1] = sfrac;
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                jp = j + 1;
                jpn = jp + n;
                ptsid[jp] = (doublereal) j + sfrac;
                if (jpn <= npt) {
                    ptsid[jpn] = (doublereal) j / (doublereal) np + sfrac;
                    temp = one / (ptsaux[(j << 1) + 1] - ptsaux[(j << 1) + 2]);
                    bmat[jp + j * bmat_dim1] = -temp + one / ptsaux[(j << 1) + 1];
                    bmat[jpn + j * bmat_dim1] = temp + one / ptsaux[(j << 1) + 2];
                    bmat[j * bmat_dim1 + 1] = -bmat[jp + j * bmat_dim1] - bmat[jpn +
                        j * bmat_dim1];
                    zmat[j * zmat_dim1 + 1] = std::sqrt(2.) / (d__1 = ptsaux[(j << 1) + 1]
                                                          * ptsaux[(j << 1) + 2], std::abs(d__1));
                    zmat[jp + j * zmat_dim1] = zmat[j * zmat_dim1 + 1] * ptsaux[(j <<
                                                                                 1) + 2] * temp;
                    zmat[jpn + j * zmat_dim1] = -zmat[j * zmat_dim1 + 1] * ptsaux[(j
                                                                                   << 1) + 1] * temp;
                } else {
                    bmat[j * bmat_dim1 + 1] = -one / ptsaux[(j << 1) + 1];
                    bmat[jp + j * bmat_dim1] = one / ptsaux[(j << 1) + 1];
                    /* Computing 2nd power */
                    d__1 = ptsaux[(j << 1) + 1];
                    bmat[j + npt + j * bmat_dim1] = -half * (d__1 * d__1);
                }
                /* L60: */
            }

            /*     Set any remaining identifiers with their nonzero elements of ZMAT. */

            if (npt >= n + np) {
                i__2 = npt;
                for (k = np << 1; k <= i__2; ++k) {
                    iw = (integer) (((doublereal) (k - np) - half) / (doublereal) (n)
                    );
                    ip = k - np - iw * n;
                    iq = ip + iw;
                    if (iq > n) {
                        iq -= n;
                    }
                    ptsid[k] = (doublereal) ip + (doublereal) iq / (doublereal) np +
                        sfrac;
                    temp = one / (ptsaux[(ip << 1) + 1] * ptsaux[(iq << 1) + 1]);
                    zmat[(k - np) * zmat_dim1 + 1] = temp;
                    zmat[ip + 1 + (k - np) * zmat_dim1] = -temp;
                    zmat[iq + 1 + (k - np) * zmat_dim1] = -temp;
                    /* L70: */
                    zmat[k + (k - np) * zmat_dim1] = temp;
                }
            }
            nrem = npt;
            kold = 1;
            knew = kopt;

            /*     Reorder the provisional points in the way that exchanges PTSID(KOLD) */
            /*     with PTSID(KNEW). */

L80:
            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                temp = bmat[kold + j * bmat_dim1];
                bmat[kold + j * bmat_dim1] = bmat[knew + j * bmat_dim1];
                /* L90: */
                bmat[knew + j * bmat_dim1] = temp;
            }
            i__2 = nptm;
            for (j = 1; j <= i__2; ++j) {
                temp = zmat[kold + j * zmat_dim1];
                zmat[kold + j * zmat_dim1] = zmat[knew + j * zmat_dim1];
                /* L100: */
                zmat[knew + j * zmat_dim1] = temp;
            }
            ptsid[kold] = ptsid[knew];
            ptsid[knew] = zero;
            w[ndim + knew] = zero;
            --nrem;
            if (knew != kopt) {
                temp = vlag[kold];
                vlag[kold] = vlag[knew];
                vlag[knew] = temp;

                /*     Update the BMAT and ZMAT matrices so that the status of the KNEW-th */
                /*     interpolation point can be changed from provisional to original. The */
                /*     branch to label 350 occurs if all the original points are reinstated. */
                /*     The nonnegative values of W(NDIM+K) are required in the search below. */

                update_(n, npt, &bmat[bmat_offset], &zmat[zmat_offset], ndim, &vlag[1],
                        beta, denom, knew, &w[1]);
                if (nrem == 0) {
                    goto L350;
                }
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L110: */
                    w[ndim + k] = (d__1 = w[ndim + k], std::abs(d__1));
                }
            }

            /*     Pick the index KNEW of an original interpolation point that has not */
            /*     yet replaced one of the provisional interpolation points, giving */
            /*     attention to the closeness to XOPT and to previous tries with KNEW. */

L120:
            dsqmin = zero;
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                if (w[ndim + k] > zero) {
                    if (dsqmin == zero || w[ndim + k] < dsqmin) {
                        knew = k;
                        dsqmin = w[ndim + k];
                    }
                }
                /* L130: */
            }
            if (dsqmin == zero) {
                goto L260;
            }

            /*     Form the W-vector of the chosen original interpolation point. */

            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                /* L140: */
                w[npt + j] = xpt[knew + j * xpt_dim1];
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                sum = zero;
                if (k == kopt) {
                } else if (ptsid[k] == zero) {
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L150: */
                        sum += w[npt + j] * xpt[k + j * xpt_dim1];
                    }
                } else {
                    ip = (integer) ptsid[k];
                    if (ip > 0) {
                        sum = w[npt + ip] * ptsaux[(ip << 1) + 1];
                    }
                    iq = (integer) ((doublereal) np * ptsid[k] - (doublereal) (ip *
                                                                               np));
                    if (iq > 0) {
                        iw = 1;
                        if (ip == 0) {
                            iw = 2;
                        }
                        sum += w[npt + iq] * ptsaux[iw + (iq << 1)];
                    }
                }
                /* L160: */
                w[k] = half * sum * sum;
            }

            /*     Calculate VLAG and BETA for the required updating of the H matrix if */
            /*     XPT(KNEW,.) is reinstated in the set of interpolation points. */

            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                sum = zero;
                i__1 = n;
                for (j = 1; j <= i__1; ++j) {
                    /* L170: */
                    sum += bmat[k + j * bmat_dim1] * w[npt + j];
                }
                /* L180: */
                vlag[k] = sum;
            }
            beta = zero;
            i__2 = nptm;
            for (j = 1; j <= i__2; ++j) {
                sum = zero;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L190: */
                    sum += zmat[k + j * zmat_dim1] * w[k];
                }
                beta -= sum * sum;
                i__1 = npt;
                for (k = 1; k <= i__1; ++k) {
                    /* L200: */
                    vlag[k] += sum * zmat[k + j * zmat_dim1];
                }
            }
            bsum = zero;
            distsq = zero;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                sum = zero;
                i__2 = npt;
                for (k = 1; k <= i__2; ++k) {
                    /* L210: */
                    sum += bmat[k + j * bmat_dim1] * w[k];
                }
                jp = j + npt;
                bsum += sum * w[jp];
                i__2 = ndim;
                for (ip = npt + 1; ip <= i__2; ++ip) {
                    /* L220: */
                    sum += bmat[ip + j * bmat_dim1] * w[ip];
                }
                bsum += sum * w[jp];
                vlag[jp] = sum;
                /* L230: */
                /* Computing 2nd power */
                d__1 = xpt[knew + j * xpt_dim1];
                distsq += d__1 * d__1;
            }
            beta = half * distsq * distsq + beta - bsum;
            vlag[kopt] += one;

            /*     KOLD is set to the index of the provisional interpolation point that is */
            /*     going to be deleted to make way for the KNEW-th original interpolation */
            /*     point. The choice of KOLD is governed by the avoidance of a small value */
            /*     of the denominator in the updating calculation of UPDATE. */

            denom = zero;
            vlmxsq = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                if (ptsid[k] != zero) {
                    hdiag = zero;
                    i__2 = nptm;
                    for (j = 1; j <= i__2; ++j) {
                        /* L240: */
                        /* Computing 2nd power */
                        d__1 = zmat[k + j * zmat_dim1];
                        hdiag += d__1 * d__1;
                    }
                    /* Computing 2nd power */
                    d__1 = vlag[k];
                    den = beta * hdiag + d__1 * d__1;
                    if (den > denom) {
                        kold = k;
                        denom = den;
                    }
                }
                /* L250: */
                /* Computing MAX */
                /* Computing 2nd power */
                d__3 = vlag[k];
                d__1 = vlmxsq, d__2 = d__3 * d__3;
                vlmxsq = std::max(d__1,d__2);
            }
            if (denom <= vlmxsq * .01) {
                w[ndim + knew] = -w[ndim + knew] - winc;
                goto L120;
            }
            goto L80;

            /*     When label 260 is reached, all the final positions of the interpolation */
            /*     points have been chosen although any changes have not been included yet */
            /*     in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart */
            /*     from the shift of XBASE, the updating of the quadratic model remains to */
            /*     be done. The following cycle through the new interpolation points begins */
            /*     by putting the new point in XPT(KPT,.) and by setting PQ(KPT) to zero, */
            /*     except that a RETURN occurs if MAXFUN prohibits another value of F. */

L260:
            i__1 = npt;
            for (kpt = 1; kpt <= i__1; ++kpt) {
                if (ptsid[kpt] == zero) {
                    goto L340;
                }
                if (nf >= maxfun) {
                    nf = -1;
                    goto L350;
                }
                ih = 0;
                i__2 = n;
                for (j = 1; j <= i__2; ++j) {
                    w[j] = xpt[kpt + j * xpt_dim1];
                    xpt[kpt + j * xpt_dim1] = zero;
                    temp = pq[kpt] * w[j];
                    i__3 = j;
                    for (i__ = 1; i__ <= i__3; ++i__) {
                        ++ih;
                        /* L270: */
                        hq[ih] += temp * w[i__];
                    }
                }
                pq[kpt] = zero;
                ip = (integer) ptsid[kpt];
                iq = (integer) ((doublereal) np * ptsid[kpt] - (doublereal) (ip * np))
                    ;
                if (ip > 0) {
                    xp = ptsaux[(ip << 1) + 1];
                    xpt[kpt + ip * xpt_dim1] = xp;
                }
                if (iq > 0) {
                    xq = ptsaux[(iq << 1) + 1];
                    if (ip == 0) {
                        xq = ptsaux[(iq << 1) + 2];
                    }
                    xpt[kpt + iq * xpt_dim1] = xq;
                }

                /*     Set VQUAD to the value of the current model at the new point. */

                vquad = fbase;
                if (ip > 0) {
                    ihp = (ip + ip * ip) / 2;
                    vquad += xp * (gopt[ip] + half * xp * hq[ihp]);
                }
                if (iq > 0) {
                    ihq = (iq + iq * iq) / 2;
                    vquad += xq * (gopt[iq] + half * xq * hq[ihq]);
                    if (ip > 0) {
                        iw = std::max(ihp,ihq) - (i__3 = ip - iq, std::abs(i__3));
                        vquad += xp * xq * hq[iw];
                    }
                }
                i__3 = npt;
                for (k = 1; k <= i__3; ++k) {
                    temp = zero;
                    if (ip > 0) {
                        temp += xp * xpt[k + ip * xpt_dim1];
                    }
                    if (iq > 0) {
                        temp += xq * xpt[k + iq * xpt_dim1];
                    }
                    /* L280: */
                    vquad += half * pq[k] * temp * temp;
                }

                /*     Calculate F at the new interpolation point, and set DIFF to the factor */
                /*     that is going to multiply the KPT-th Lagrange function when the model */
                /*     is updated to provide interpolation to the new function value. */

                i__3 = n;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    /* Computing MIN */
                    /* Computing MAX */
                    d__3 = xl[i__], d__4 = xbase[i__] + xpt[kpt + i__ * xpt_dim1];
                    d__1 = std::max(d__3,d__4), d__2 = xu[i__];
                    w[i__] = std::min(d__1,d__2);
                    if (xpt[kpt + i__ * xpt_dim1] == sl[i__]) {
                        w[i__] = xl[i__];
                    }
                    if (xpt[kpt + i__ * xpt_dim1] == su[i__]) {
                        w[i__] = xu[i__];
                    }
                    /* L290: */
                }
                ++(nf);
                f = calfun(mat(&w[1],n));
                fval[kpt] = f;
                if (f < fval[kopt]) {
                    kopt = kpt;
                }
                diff = f - vquad;

                /*     Update the quadratic model. The RETURN from the subroutine occurs when */
                /*     all the new interpolation points are included in the model. */

                i__3 = n;
                for (i__ = 1; i__ <= i__3; ++i__) {
                    /* L310: */
                    gopt[i__] += diff * bmat[kpt + i__ * bmat_dim1];
                }
                i__3 = npt;
                for (k = 1; k <= i__3; ++k) {
                    sum = zero;
                    i__2 = nptm;
                    for (j = 1; j <= i__2; ++j) {
                        /* L320: */
                        sum += zmat[k + j * zmat_dim1] * zmat[kpt + j * zmat_dim1];
                    }
                    temp = diff * sum;
                    if (ptsid[k] == zero) {
                        pq[k] += temp;
                    } else {
                        ip = (integer) ptsid[k];
                        iq = (integer) ((doublereal) np * ptsid[k] - (doublereal) (ip
                                                                                   * np));
                        ihq = (iq * iq + iq) / 2;
                        if (ip == 0) {
                            /* Computing 2nd power */
                            d__1 = ptsaux[(iq << 1) + 2];
                            hq[ihq] += temp * (d__1 * d__1);
                        } else {
                            ihp = (ip * ip + ip) / 2;
                            /* Computing 2nd power */
                            d__1 = ptsaux[(ip << 1) + 1];
                            hq[ihp] += temp * (d__1 * d__1);
                            if (iq > 0) {
                                /* Computing 2nd power */
                                d__1 = ptsaux[(iq << 1) + 1];
                                hq[ihq] += temp * (d__1 * d__1);
                                iw = std::max(ihp,ihq) - (i__2 = iq - ip, std::abs(i__2));
                                hq[iw] += temp * ptsaux[(ip << 1) + 1] * ptsaux[(iq <<
                                                                                 1) + 1];
                            }
                        }
                    }
                    /* L330: */
                }
                ptsid[kpt] = zero;
L340:
                ;
            }
L350:
            ;
        } /* rescue_ */

    // ----------------------------------------------------------------------------------------

        void trsbox_(
            const integer n,
            const integer npt,
            const doublereal *xpt,
            const doublereal *xopt,
            const doublereal *gopt,
            const doublereal *hq,
            const doublereal *pq,
            const doublereal *sl,
            const doublereal *su,
            const doublereal delta,
            doublereal *xnew,
            doublereal *d__,
            doublereal *gnew,
            doublereal *xbdi,
            doublereal *s,
            doublereal *hs,
            doublereal *hred,
            doublereal *dsq,
            doublereal *crvmin
        ) const
        {
            /* System generated locals */
            integer xpt_dim1, xpt_offset, i__1, i__2;
            doublereal d__1, d__2, d__3, d__4;

            /* Local variables */
            integer i__, j, k, ih;
            doublereal ds;
            integer iu;
            doublereal dhd, dhs, cth, one, shs, sth, ssq, half, beta, sdec, blen;
            integer iact = 0, nact = 0;
            doublereal angt, qred;
            integer isav;
            doublereal temp = 0, zero = 0, xsav = 0, xsum = 0, angbd = 0, dredg = 0, sredg = 0;
            integer iterc;
            doublereal resid = 0, delsq = 0, ggsav = 0, tempa = 0, tempb = 0,
                       redmax = 0, dredsq = 0, redsav = 0, onemin = 0, gredsq = 0, rednew = 0;
            integer itcsav = 0;
            doublereal rdprev = 0, rdnext = 0, stplen = 0, stepsq = 0;
            integer itermax = 0;


            /*     The arguments N, NPT, XPT, XOPT, GOPT, HQ, PQ, SL and SU have the same */
            /*       meanings as the corresponding arguments of BOBYQB. */
            /*     DELTA is the trust region radius for the present calculation, which */
            /*       seeks a small value of the quadratic model within distance DELTA of */
            /*       XOPT subject to the bounds on the variables. */
            /*     XNEW will be set to a new vector of variables that is approximately */
            /*       the one that minimizes the quadratic model within the trust region */
            /*       subject to the SL and SU constraints on the variables. It satisfies */
            /*       as equations the bounds that become active during the calculation. */
            /*     D is the calculated trial step from XOPT, generated iteratively from an */
            /*       initial value of zero. Thus XNEW is XOPT+D after the final iteration. */
            /*     GNEW holds the gradient of the quadratic model at XOPT+D. It is updated */
            /*       when D is updated. */
            /*     XBDI is a working space vector. For I=1,2,...,N, the element XBDI(I) is */
            /*       set to -1.0, 0.0, or 1.0, the value being nonzero if and only if the */
            /*       I-th variable has become fixed at a bound, the bound being SL(I) or */
            /*       SU(I) in the case XBDI(I)=-1.0 or XBDI(I)=1.0, respectively. This */
            /*       information is accumulated during the construction of XNEW. */
            /*     The arrays S, HS and HRED are also used for working space. They hold the */
            /*       current search direction, and the changes in the gradient of Q along S */
            /*       and the reduced D, respectively, where the reduced D is the same as D, */
            /*       except that the components of the fixed variables are zero. */
            /*     DSQ will be set to the square of the length of XNEW-XOPT. */
            /*     CRVMIN is set to zero if D reaches the trust region boundary. Otherwise */
            /*       it is set to the least curvature of H that occurs in the conjugate */
            /*       gradient searches that are not restricted by any constraints. The */
            /*       value CRVMIN=-1.0D0 is set, however, if all of these searches are */
            /*       constrained. */

            /*     A version of the truncated conjugate gradient is applied. If a line */
            /*     search is restricted by a constraint, then the procedure is restarted, */
            /*     the values of the variables that are at their bounds being fixed. If */
            /*     the trust region boundary is reached, then further changes may be made */
            /*     to D, each one being in the two dimensional space that is spanned */
            /*     by the current D and the gradient of Q at XOPT+D, staying on the trust */
            /*     region boundary. Termination occurs when the reduction in Q seems to */
            /*     be close to the greatest reduction that can be achieved. */

            /*     Set some constants. */

            /* Parameter adjustments */
            xpt_dim1 = npt;
            xpt_offset = 1 + xpt_dim1;
            xpt -= xpt_offset;
            --xopt;
            --gopt;
            --hq;
            --pq;
            --sl;
            --su;
            --xnew;
            --d__;
            --gnew;
            --xbdi;
            --s;
            --hs;
            --hred;

            /* Function Body */
            half = .5;
            one = 1.;
            onemin = -1.;
            zero = 0.;

            /*     The sign of GOPT(I) gives the sign of the change to the I-th variable */
            /*     that will reduce Q from its value at XOPT. Thus XBDI(I) shows whether */
            /*     or not to fix the I-th variable at one of its bounds initially, with */
            /*     NACT being set to the number of fixed variables. D and GNEW are also */
            /*     set for the first iteration. DELSQ is the upper bound on the sum of */
            /*     squares of the free variables. QRED is the reduction in Q so far. */

            iterc = 0;
            nact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                xbdi[i__] = zero;
                if (xopt[i__] <= sl[i__]) {
                    if (gopt[i__] >= zero) {
                        xbdi[i__] = onemin;
                    }
                } else if (xopt[i__] >= su[i__]) {
                    if (gopt[i__] <= zero) {
                        xbdi[i__] = one;
                    }
                }
                if (xbdi[i__] != zero) {
                    ++nact;
                }
                d__[i__] = zero;
                /* L10: */
                gnew[i__] = gopt[i__];
            }
            delsq = delta * delta;
            qred = zero;
            *crvmin = onemin;

            /*     Set the next search direction of the conjugate gradient method. It is */
            /*     the steepest descent direction initially and when the iterations are */
            /*     restarted because a variable has just been fixed by a bound, and of */
            /*     course the components of the fixed variables are zero. ITERMAX is an */
            /*     upper bound on the indices of the conjugate gradient iterations. */

L20:
            beta = zero;
L30:
            stepsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] != zero) {
                    s[i__] = zero;
                } else if (beta == zero) {
                    s[i__] = -gnew[i__];
                } else {
                    s[i__] = beta * s[i__] - gnew[i__];
                }
                /* L40: */
                /* Computing 2nd power */
                d__1 = s[i__];
                stepsq += d__1 * d__1;
            }
            if (stepsq == zero) {
                goto L190;
            }
            if (beta == zero) {
                gredsq = stepsq;
                itermax = iterc + n - nact;
            }
            if (gredsq * delsq <= qred * 1e-4 * qred) {
                goto L190;
            }

            /*     Multiply the search direction by the second derivative matrix of Q and */
            /*     calculate some scalars for the choice of steplength. Then set BLEN to */
            /*     the length of the the step to the trust region boundary and STPLEN to */
            /*     the steplength, ignoring the simple bounds. */

            goto L210;
L50:
            resid = delsq;
            ds = zero;
            shs = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    resid -= d__1 * d__1;
                    ds += s[i__] * d__[i__];
                    shs += s[i__] * hs[i__];
                }
                /* L60: */
            }
            if (resid <= zero) {
                goto L90;
            }
            temp = std::sqrt(stepsq * resid + ds * ds);
            if (ds < zero) {
                blen = (temp - ds) / stepsq;
            } else {
                blen = resid / (temp + ds);
            }
            stplen = blen;
            if (shs > zero) {
                /* Computing MIN */
                d__1 = blen, d__2 = gredsq / shs;
                stplen = std::min(d__1,d__2);
            }

            /*     Reduce STPLEN if necessary in order to preserve the simple bounds, */
            /*     letting IACT be the index of the new constrained variable. */

            iact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (s[i__] != zero) {
                    xsum = xopt[i__] + d__[i__];
                    if (s[i__] > zero) {
                        temp = (su[i__] - xsum) / s[i__];
                    } else {
                        temp = (sl[i__] - xsum) / s[i__];
                    }
                    if (temp < stplen) {
                        stplen = temp;
                        iact = i__;
                    }
                }
                /* L70: */
            }

            /*     Update CRVMIN, GNEW and D. Set SDEC to the decrease that occurs in Q. */

            sdec = zero;
            if (stplen > zero) {
                ++iterc;
                temp = shs / stepsq;
                if (iact == 0 && temp > zero) {
                    *crvmin = std::min(*crvmin,temp);
                    if (*crvmin == onemin) {
                        *crvmin = temp;
                    }
                }
                ggsav = gredsq;
                gredsq = zero;
                i__1 = n;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    gnew[i__] += stplen * hs[i__];
                    if (xbdi[i__] == zero) {
                        /* Computing 2nd power */
                        d__1 = gnew[i__];
                        gredsq += d__1 * d__1;
                    }
                    /* L80: */
                    d__[i__] += stplen * s[i__];
                }
                /* Computing MAX */
                d__1 = stplen * (ggsav - half * stplen * shs);
                sdec = std::max(d__1,zero);
                qred += sdec;
            }

            /*     Restart the conjugate gradient method if it has hit a new bound. */

            if (iact > 0) {
                ++nact;
                xbdi[iact] = one;
                if (s[iact] < zero) {
                    xbdi[iact] = onemin;
                }
                /* Computing 2nd power */
                d__1 = d__[iact];
                delsq -= d__1 * d__1;
                if (delsq <= zero) {
                    goto L90;
                }
                goto L20;
            }

            /*     If STPLEN is less than BLEN, then either apply another conjugate */
            /*     gradient iteration or RETURN. */

            if (stplen < blen) {
                if (iterc == itermax) {
                    goto L190;
                }
                if (sdec <= qred * .01) {
                    goto L190;
                }
                beta = gredsq / ggsav;
                goto L30;
            }
L90:
            *crvmin = zero;

            /*     Prepare for the alternative iteration by calculating some scalars */
            /*     and by multiplying the reduced D by the second derivative matrix of */
            /*     Q, where S holds the reduced D in the call of GGMULT. */

L100:
            if (nact >= n - 1) {
                goto L190;
            }
            dredsq = zero;
            dredg = zero;
            gredsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    dredsq += d__1 * d__1;
                    dredg += d__[i__] * gnew[i__];
                    /* Computing 2nd power */
                    d__1 = gnew[i__];
                    gredsq += d__1 * d__1;
                    s[i__] = d__[i__];
                } else {
                    s[i__] = zero;
                }
                /* L110: */
            }
            itcsav = iterc;
            goto L210;

            /*     Let the search direction S be a linear combination of the reduced D */
            /*     and the reduced G that is orthogonal to the reduced D. */

L120:
            ++iterc;
            temp = gredsq * dredsq - dredg * dredg;
            if (temp <= qred * 1e-4 * qred) {
                goto L190;
            }
            temp = std::sqrt(temp);
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    s[i__] = (dredg * d__[i__] - dredsq * gnew[i__]) / temp;
                } else {
                    s[i__] = zero;
                }
                /* L130: */
            }
            sredg = -temp;

            /*     By considering the simple bounds on the variables, calculate an upper */
            /*     bound on the tangent of half the angle of the alternative iteration, */
            /*     namely ANGBD, except that, if already a free variable has reached a */
            /*     bound, there is a branch back to label 100 after fixing that variable. */

            angbd = one;
            iact = 0;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    tempa = xopt[i__] + d__[i__] - sl[i__];
                    tempb = su[i__] - xopt[i__] - d__[i__];
                    if (tempa <= zero) {
                        ++nact;
                        xbdi[i__] = onemin;
                        goto L100;
                    } else if (tempb <= zero) {
                        ++nact;
                        xbdi[i__] = one;
                        goto L100;
                    }
                    /* Computing 2nd power */
                    d__1 = d__[i__];
                    /* Computing 2nd power */
                    d__2 = s[i__];
                    ssq = d__1 * d__1 + d__2 * d__2;
                    /* Computing 2nd power */
                    d__1 = xopt[i__] - sl[i__];
                    temp = ssq - d__1 * d__1;
                    if (temp > zero) {
                        temp = std::sqrt(temp) - s[i__];
                        if (angbd * temp > tempa) {
                            angbd = tempa / temp;
                            iact = i__;
                            xsav = onemin;
                        }
                    }
                    /* Computing 2nd power */
                    d__1 = su[i__] - xopt[i__];
                    temp = ssq - d__1 * d__1;
                    if (temp > zero) {
                        temp = std::sqrt(temp) + s[i__];
                        if (angbd * temp > tempb) {
                            angbd = tempb / temp;
                            iact = i__;
                            xsav = one;
                        }
                    }
                }
                /* L140: */
            }

            /*     Calculate HHD and some curvatures for the alternative iteration. */

            goto L210;
L150:
            shs = zero;
            dhs = zero;
            dhd = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                if (xbdi[i__] == zero) {
                    shs += s[i__] * hs[i__];
                    dhs += d__[i__] * hs[i__];
                    dhd += d__[i__] * hred[i__];
                }
                /* L160: */
            }

            /*     Seek the greatest reduction in Q for a range of equally spaced values */
            /*     of ANGT in [0,ANGBD], where ANGT is the tangent of half the angle of */
            /*     the alternative iteration. */

            redmax = zero;
            isav = 0;
            redsav = zero;
            iu = (integer) (angbd * 17. + 3.1);
            i__1 = iu;
            for (i__ = 1; i__ <= i__1; ++i__) {
                angt = angbd * (doublereal) i__ / (doublereal) iu;
                sth = (angt + angt) / (one + angt * angt);
                temp = shs + angt * (angt * dhd - dhs - dhs);
                rednew = sth * (angt * dredg - sredg - half * sth * temp);
                if (rednew > redmax) {
                    redmax = rednew;
                    isav = i__;
                    rdprev = redsav;
                } else if (i__ == isav + 1) {
                    rdnext = rednew;
                }
                /* L170: */
                redsav = rednew;
            }

            /*     Return if the reduction is zero. Otherwise, set the sine and cosine */
            /*     of the angle of the alternative iteration, and calculate SDEC. */

            if (isav == 0) {
                goto L190;
            }
            if (isav < iu) {
                temp = (rdnext - rdprev) / (redmax + redmax - rdprev - rdnext);
                angt = angbd * ((doublereal) isav + half * temp) / (doublereal) iu;
            }
            cth = (one - angt * angt) / (one + angt * angt);
            sth = (angt + angt) / (one + angt * angt);
            temp = shs + angt * (angt * dhd - dhs - dhs);
            sdec = sth * (angt * dredg - sredg - half * sth * temp);
            if (sdec <= zero) {
                goto L190;
            }

            /*     Update GNEW, D and HRED. If the angle of the alternative iteration */
            /*     is restricted by a bound on a free variable, that variable is fixed */
            /*     at the bound. */

            dredg = zero;
            gredsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                gnew[i__] = gnew[i__] + (cth - one) * hred[i__] + sth * hs[i__];
                if (xbdi[i__] == zero) {
                    d__[i__] = cth * d__[i__] + sth * s[i__];
                    dredg += d__[i__] * gnew[i__];
                    /* Computing 2nd power */
                    d__1 = gnew[i__];
                    gredsq += d__1 * d__1;
                }
                /* L180: */
                hred[i__] = cth * hred[i__] + sth * hs[i__];
            }
            qred += sdec;
            if (iact > 0 && isav == iu) {
                ++nact;
                xbdi[iact] = xsav;
                goto L100;
            }

            /*     If SDEC is sufficiently small, then RETURN after setting XNEW to */
            /*     XOPT+D, giving careful attention to the bounds. */

            if (sdec > qred * .01) {
                goto L120;
            }
L190:
            *dsq = zero;
            i__1 = n;
            for (i__ = 1; i__ <= i__1; ++i__) {
                /* Computing MAX */
                /* Computing MIN */
                d__3 = xopt[i__] + d__[i__], d__4 = su[i__];
                d__1 = std::min(d__3,d__4), d__2 = sl[i__];
                xnew[i__] = std::max(d__1,d__2);
                if (xbdi[i__] == onemin) {
                    xnew[i__] = sl[i__];
                }
                if (xbdi[i__] == one) {
                    xnew[i__] = su[i__];
                }
                d__[i__] = xnew[i__] - xopt[i__];
                /* L200: */
                /* Computing 2nd power */
                d__1 = d__[i__];
                *dsq += d__1 * d__1;
            }
            return;
            /*     The following instructions multiply the current S-vector by the second */
            /*     derivative matrix of the quadratic model, putting the product in HS. */
            /*     They are reached from three different parts of the software above and */
            /*     they can be regarded as an external subroutine. */

L210:
            ih = 0;
            i__1 = n;
            for (j = 1; j <= i__1; ++j) {
                hs[j] = zero;
                i__2 = j;
                for (i__ = 1; i__ <= i__2; ++i__) {
                    ++ih;
                    if (i__ < j) {
                        hs[j] += hq[ih] * s[i__];
                    }
                    /* L220: */
                    hs[i__] += hq[ih] * s[j];
                }
            }
            i__2 = npt;
            for (k = 1; k <= i__2; ++k) {
                if (pq[k] != zero) {
                    temp = zero;
                    i__1 = n;
                    for (j = 1; j <= i__1; ++j) {
                        /* L230: */
                        temp += xpt[k + j * xpt_dim1] * s[j];
                    }
                    temp *= pq[k];
                    i__1 = n;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        /* L240: */
                        hs[i__] += temp * xpt[k + i__ * xpt_dim1];
                    }
                }
                /* L250: */
            }
            if (*crvmin != zero) {
                goto L50;
            }
            if (iterc > itcsav) {
                goto L150;
            }
            i__2 = n;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L260: */
                hred[i__] = hs[i__];
            }
            goto L120;
        } /* trsbox_ */

    // ----------------------------------------------------------------------------------------

        void update_(
            const integer n,
            const integer npt,
            doublereal *bmat,
            doublereal *zmat,
            const integer ndim,
            doublereal *vlag,
            const doublereal beta,
            const doublereal denom,
            const integer knew,
            doublereal *w
        ) const
        {
            /* System generated locals */
            integer bmat_dim1, bmat_offset, zmat_dim1, zmat_offset, i__1, i__2;
            doublereal d__1, d__2, d__3;

            /* Local variables */
            integer i__, j, k, jp;
            doublereal one, tau, temp;
            integer nptm;
            doublereal zero, alpha, tempa, tempb, ztest;


            /*     The arrays BMAT and ZMAT are updated, as required by the new position */
            /*     of the interpolation point that has the index KNEW. The vector VLAG has */
            /*     N+NPT components, set on entry to the first NPT and last N components */
            /*     of the product Hw in equation (4.11) of the Powell (2006) paper on */
            /*     NEWUOA. Further, BETA is set on entry to the value of the parameter */
            /*     with that name, and DENOM is set to the denominator of the updating */
            /*     formula. Elements of ZMAT may be treated as zero if their moduli are */
            /*     at most ZTEST. The first NDIM elements of W are used for working space. */

            /*     Set some constants. */

            /* Parameter adjustments */
            zmat_dim1 = npt;
            zmat_offset = 1 + zmat_dim1;
            zmat -= zmat_offset;
            bmat_dim1 = ndim;
            bmat_offset = 1 + bmat_dim1;
            bmat -= bmat_offset;
            --vlag;
            --w;

            /* Function Body */
            one = 1.;
            zero = 0.;
            nptm = npt - n - 1;
            ztest = zero;
            i__1 = npt;
            for (k = 1; k <= i__1; ++k) {
                i__2 = nptm;
                for (j = 1; j <= i__2; ++j) {
                    /* L10: */
                    /* Computing MAX */
                    d__2 = ztest, d__3 = (d__1 = zmat[k + j * zmat_dim1], std::abs(d__1));
                    ztest = std::max(d__2,d__3);
                }
            }
            ztest *= 1e-20;

            /*     Apply the rotations that put zeros in the KNEW-th row of ZMAT. */

            i__2 = nptm;
            for (j = 2; j <= i__2; ++j) {
                if ((d__1 = zmat[knew + j * zmat_dim1], std::abs(d__1)) > ztest) {
                    /* Computing 2nd power */
                    d__1 = zmat[knew + zmat_dim1];
                    /* Computing 2nd power */
                    d__2 = zmat[knew + j * zmat_dim1];
                    temp = std::sqrt(d__1 * d__1 + d__2 * d__2);
                    tempa = zmat[knew + zmat_dim1] / temp;
                    tempb = zmat[knew + j * zmat_dim1] / temp;
                    i__1 = npt;
                    for (i__ = 1; i__ <= i__1; ++i__) {
                        temp = tempa * zmat[i__ + zmat_dim1] + tempb * zmat[i__ + j *
                            zmat_dim1];
                        zmat[i__ + j * zmat_dim1] = tempa * zmat[i__ + j * zmat_dim1]
                            - tempb * zmat[i__ + zmat_dim1];
                        /* L20: */
                        zmat[i__ + zmat_dim1] = temp;
                    }
                }
                zmat[knew + j * zmat_dim1] = zero;
                /* L30: */
            }

            /*     Put the first NPT components of the KNEW-th column of HLAG into W, */
            /*     and calculate the parameters of the updating formula. */

            i__2 = npt;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w[i__] = zmat[knew + zmat_dim1] * zmat[i__ + zmat_dim1];
                /* L40: */
            }
            alpha = w[knew];
            tau = vlag[knew];
            vlag[knew] -= one;

            /*     Complete the updating of ZMAT. */

            temp = std::sqrt(denom);
            tempb = zmat[knew + zmat_dim1] / temp;
            tempa = tau / temp;
            i__2 = npt;
            for (i__ = 1; i__ <= i__2; ++i__) {
                /* L50: */
                zmat[i__ + zmat_dim1] = tempa * zmat[i__ + zmat_dim1] - tempb * vlag[
                    i__];
            }

            /*     Finally, update the matrix BMAT. */

            i__2 = n;
            for (j = 1; j <= i__2; ++j) {
                jp = npt + j;
                w[jp] = bmat[knew + j * bmat_dim1];
                tempa = (alpha * vlag[jp] - tau * w[jp]) / denom;
                tempb = (-(beta) * w[jp] - tau * vlag[jp]) / denom;
                i__1 = jp;
                for (i__ = 1; i__ <= i__1; ++i__) {
                    bmat[i__ + j * bmat_dim1] = bmat[i__ + j * bmat_dim1] + tempa *
                        vlag[i__] + tempb * w[i__];
                    if (i__ > npt) {
                        bmat[jp + (i__ - npt) * bmat_dim1] = bmat[i__ + j *
                            bmat_dim1];
                    }
                    /* L60: */
                }
            }
        } /* update_ */
    };

// ----------------------------------------------------------------------------------------

    template <
        typename funct,
        typename T,
        typename U
        >
    double find_min_bobyqaConvgFail (
        const funct& f,
        T& x,
        long npt,
        const U& x_lower,
        const U& x_upper,
        const double rho_begin,
        const double rho_end,
        const long max_f_evals,  int& ConvgFail
    )
    {
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        // check the requirements.  Also split the assert up so that the error message isn't huge.
        DLIB_CASSERT(is_col_vector(x) && is_col_vector(x_lower) && is_col_vector(x_upper) &&
                    x.size() == x_lower.size() && x_lower.size() == x_upper.size() &&
                    x.size() > 1 && max_f_evals > 1,
            "\tvoid find_min_bobyqa()"
            << "\n\t Invalid arguments have been given to this function"
            << "\n\t is_col_vector(x):       " << is_col_vector(x)
            << "\n\t is_col_vector(x_lower): " << is_col_vector(x_lower)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t x.size():               " << x.size()
            << "\n\t x_lower.size():         " << x_lower.size()
            << "\n\t x_upper.size():         " << x_upper.size()
            << "\n\t max_f_evals:            " << max_f_evals
        );


        DLIB_CASSERT(x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2 &&
                    0 < rho_end && rho_end < rho_begin &&
                    min(x_upper - x_lower) > 2*rho_begin &&
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0,
            "\tvoid find_min_bobyqa()"
            << "\n\t Invalid arguments have been given to this function"
            << "\n\t ntp in valid range: " << (x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2)
            << "\n\t npt:                " << npt
            << "\n\t rho_begin:          " << rho_begin
            << "\n\t rho_end:            " << rho_end
            << "\n\t min(x_upper - x_lower) > 2*rho_begin:           " << (min(x_upper - x_lower) > 2*rho_begin)
            << "\n\t min(x - x_lower) >= 0 && min(x_upper - x) >= 0: " << (min(x - x_lower) >= 0 && min(x_upper - x) >= 0)
        );
       /*

        if(x.size() + 2 <= npt && npt <= (x.size()+1)*(x.size()+2)/2 &&
                    0 < rho_end && rho_end < rho_begin &&
                    min(x_upper - x_lower) > 2*rho_begin &&
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0)
        {
          std::cout<<"Invalid arguments given for find_min_bobyqaConvgFail!!!\n";
          std::cout<<"0 < rho_end && rho_end < rho_begin &&\
                    min(x_upper - x_lower) > 2*rho_begin &&\
                    min(x - x_lower) >= 0 && min(x_upper - x) >= 0 fails\n";
        ConvgFail = 1; return 0;
        }
        */


        bobyqa_implementationConvgFail impl;
        return impl.find_min(f, x, npt, x_lower, x_upper, rho_begin, rho_end, max_f_evals,ConvgFail);
    }



  template <
        typename search_strategy_type,
        typename stop_strategy_type,
        typename funct,
        typename funct_der,
        typename T,
        typename EXP1,
        typename EXP2
        >
    double find_min_box_constrainedConvFail (
        search_strategy_type search_strategy,
        stop_strategy_type stop_strategy,
        const funct& f,
        const funct_der& der,
        T& x,
        const matrix_exp<EXP1>& x_lower,
        const matrix_exp<EXP2>& x_upper, int& ConvgFail, const int& MaxIteration
    )
    {
        /*
            The implementation of this function is more or less based on the discussion in
            the paper Projected Newton-type Methods in Machine Learning by Mark Schmidt, et al.
        */

        // make sure the requires clause is not violated
        COMPILE_TIME_ASSERT(is_matrix<T>::value);
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        DLIB_CASSERT (
            is_col_vector(x) && is_col_vector(x_lower) && is_col_vector(x_upper) &&
            x.size() == x_lower.size() && x.size() == x_upper.size(),
            "\tdouble find_min_box_constrained()"
            << "\n\t The inputs to this function must be equal length column vectors."
            << "\n\t is_col_vector(x):       " << is_col_vector(x)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t is_col_vector(x_upper): " << is_col_vector(x_upper)
            << "\n\t x.size():               " << x.size()
            << "\n\t x_lower.size():         " << x_lower.size()
            << "\n\t x_upper.size():         " << x_upper.size()
        );
        DLIB_ASSERT (
            min(x_upper-x_lower) > 0,
            "\tdouble find_min_box_constrained()"
            << "\n\t You have to supply proper box constraints to this function."
            << "\n\r min(x_upper-x_lower): " << min(x_upper-x_lower)
        );


        T g, s;
        double f_value = f(x);
        g = der(x);

        if (!is_finite(f_value))
          // throw error("The objective function generated non-finite outputs");
		{  ConvgFail=1;   return 0; }
        if (!is_finite(g))
         //   throw error("The objective function generated non-finite outputs");
      {  ConvgFail=1;   return 0; }

        // gap_eps determines how close we have to get to a bound constraint before we
        // start basically dropping it from the optimization and consider it to be an
        // active constraint.
        const double gap_eps = 1e-8;

        double last_alpha = 1;  int IterDone=0;
        while(stop_strategy.should_continue_search(x, f_value, g))
        {
if(IterDone>MaxIteration)  {  ConvgFail=1;   return 0; }

            s = search_strategy.get_next_direction(x, f_value, zero_bounded_variables(gap_eps, g, x, g, x_lower, x_upper));
            s = gap_step_assign_bounded_variables(gap_eps, s, x, g, x_lower, x_upper);

            double alpha = backtracking_line_search(
                        make_line_search_function(clamp_function(f,x_lower,x_upper), x, s, f_value),
                        f_value,
                        dot(g,s), // compute gradient for the line search
                        last_alpha,
                        search_strategy.get_wolfe_rho(),
                        search_strategy.get_max_line_search_iterations());

            // Do a trust region style thing for alpha.  The idea is that if we take a
            // small step then we are likely to take another small step.  So we reuse the
            // alpha from the last iteration unless the line search didn't shrink alpha at
            // all, in that case, we start with a bigger alpha next time.
            if (alpha == last_alpha)
                last_alpha = std::min(last_alpha*10,1.0);
            else
                last_alpha = alpha;

            // Take the search step indicated by the above line search
            x = clamp(x + alpha*s, x_lower, x_upper);
            g = der(x);

  if (!is_finite(f_value))
            // throw error("The objective function generated non-finite outputs");
			 { ConvgFail=1; return 0; }
 if (!is_finite(g))
             // throw error("The objective function generated non-finite outputs");
  { ConvgFail=1; return 0; }

	      IterDone=IterDone+1;
        }

        return f_value;
    }


    template <
        typename search_strategy_type,
        typename stop_strategy_type,
        typename funct,
        typename funct_der,
        typename T
        >
    double find_minConvFail(
        search_strategy_type search_strategy,
        stop_strategy_type stop_strategy,
        const funct& f,
        const funct_der& der,
        T& x,
        double min_f, int& ConvFail, const int& MaxIter
    )
    {
        COMPILE_TIME_ASSERT(is_matrix<T>::value);
        // The starting point (i.e. x) must be a column vector.
        COMPILE_TIME_ASSERT(T::NC <= 1);

        DLIB_CASSERT (
            is_col_vector(x),
            "\tdouble find_min()"
            << "\n\tYou have to supply column vectors to this function"
            << "\n\tx.nc():    " << x.nc()
        );

        T g, s;

        double f_value = f(x);
        g = der(x);

        if (!is_finite(f_value))
        // throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
        if (!is_finite(g))
         //  throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }

      int IterCounts=0;  int MaxIterPlus1=MaxIter+1;
        while(stop_strategy.should_continue_search(x, f_value, g) && f_value > min_f )
        {
            IterCounts = IterCounts+1;
       if( IterCounts > MaxIterPlus1)  { ConvFail=1; return 0;  }
            s = search_strategy.get_next_direction(x, f_value, g);

            double alpha = line_search(
                        make_line_search_function(f,x,s, f_value),
                        f_value,
                        make_line_search_function(der,x,s, g),
                        dot(g,s), // compute initial gradient for the line search
                        search_strategy.get_wolfe_rho(), search_strategy.get_wolfe_sigma(), min_f,
                        search_strategy.get_max_line_search_iterations());

            // Take the search step indicated by the above line search
            x += alpha*s;

            if (!is_finite(f_value))
             // throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
            if (!is_finite(g))
            //  throw error("The objective function generated non-finite outputs");
{  ConvFail=1; return 0;  }
        }

        return f_value;
    }




   template <
        typename stop_strategy_type,
        typename funct_model
        >
    double find_min_trust_regionConvFail (
        stop_strategy_type stop_strategy,
        const funct_model& model,
        typename funct_model::column_vector& x,\
         int& ConvFail, const int& MaxIteration,\
        double radius = 1
    )
    {
        /*
            This is an implementation of algorithm 4.1(Trust Region)
            from the book Numerical Optimization by Nocedal and Wright.
        */

        // make sure requires clause is not broken
        DLIB_ASSERT(is_col_vector(x) && radius > 0,
            "\t double find_min_trust_region()"
            << "\n\t invalid arguments were given to this function"
            << "\n\t is_col_vector(x): " << is_col_vector(x)
            << "\n\t radius:           " << radius
            );

        const double initial_radius = radius;

        typedef typename funct_model::column_vector T;
        typedef typename T::type type;

        typename funct_model::general_matrix h;
        typename funct_model::column_vector g, p, d;
        type f_value = model(x);

        model.get_derivative_and_hessian(x,g,h);

 // DLIB_ASSERT(is_finite(x), "The objective function generated non-finite outputs");

 if( !is_finite(x) )  {ConvFail=1; return 0; }

 // DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");

 if( !is_finite(g) )  {ConvFail=1; return 0; }

//  DLIB_ASSERT(is_finite(h), "The objective function generated non-finite outputs");

 if( !is_finite(h) )  {ConvFail=1; return 0; }

        // Sometimes the loop below won't modify x because the trust region step failed.
        // This bool tells us when we are in that case.
        bool stale_x = false;
   int IterDone=0;
        while(stale_x || stop_strategy.should_continue_search(x, f_value, g))
        {

     if(IterDone>MaxIteration) {ConvFail=1; return 0; }
            const unsigned long iter = solve_trust_region_subproblem(h,
                                                                     g,
                                                                     radius,
                                                                     p,
                                                                     0.1,
                                                                     20);


            const type new_f_value = model(x+p);
            const type predicted_improvement = -0.5*trans(p)*h*p - trans(g)*p;
            const type measured_improvement = (f_value - new_f_value);

            // If the sub-problem can't find a way to improve then stop.  This only happens when p is essentially 0.
            if (std::abs(predicted_improvement) <= std::abs(measured_improvement)*std::numeric_limits<type>::epsilon())
                break;

            // predicted_improvement shouldn't be negative but it might be if something went
            // wrong in the trust region solver.  So put abs() here to guard against that.  This
            // way the sign of rho is determined only by the sign of measured_improvement.
            const type rho = measured_improvement/std::abs(predicted_improvement);



            if (rho < 0.25)
            {
                radius *= 0.25;

                // something has gone horribly wrong if the radius has shrunk to zero.  So just
                // give up if that happens.
                if (radius <= initial_radius*std::numeric_limits<double>::epsilon())
                    break;
            }
            else
            {
                // if rho > 0.75 and we are being checked by the radius
                if (rho > 0.75 && iter > 1)
                {
                    radius = std::min<type>(1000,  2*radius);
                }
            }

            if (rho > 0)
            {
                x = x + p;
                f_value = new_f_value;
                model.get_derivative_and_hessian(x,g,h);
                stale_x = false;
            }
            else
            {
                stale_x = true;
            }


//  DLIB_ASSERT(is_finite(x), "The objective function generated non-finite outputs");

 if( !is_finite(x) )  {ConvFail=1; return 0; }


//  DLIB_ASSERT(is_finite(g), "The objective function generated non-finite outputs");

if( !is_finite(g) )  {ConvFail=1; return 0; }


// DLIB_ASSERT(is_finite(h), "The objective function generated non-finite outputs");

if( !is_finite(h) )  {ConvFail=1; return 0; }

    IterDone=IterDone+1;

        }

        return f_value;
    }



/*  Compute the machine epsilon.
   The epsilon function returns the machine epsilon: the smallest number that, when summed with 1, produces a value greater than one.
  */
  // alternatively use std::numeric_limits::epsilon
template <typename T>
T epsilonGet()
{
    T eps, del, neweps;
    del    = (T) 0.5;
    eps    = (T) 0.0;
    neweps = (T) 1.0;

    while ( del > 0 )
    {
        if ( 1 + neweps > 1 )    /* Then the value might be too large */
        {
            eps = neweps;    /* ...save the current value... */
            neweps -= del;    /* ...and decrement a bit */
        }
        else          /* Then the value is too small */
        {
            neweps += del;    /* ...so increment it */
        }
        del *= 0.5;      /* Reduce the adjustment by half */
    }
    return eps;
}

const double h2 = std::sqrt(epsilonGet<double>());
const double h = std::sqrt(h2);

// numerical gradient
dbcolvec gradcdif(double (*FUNCTOR) (const dbcolvec& ), const dbcolvec& th)
{

    int dim=th.nr();
    dbcolvec grad(dim);
    dbcolvec ei(dim);
    dbcolvec eiF(dim);
    dbcolvec eiB(dim);
    for(int i=0; i<dim; ++i)
    {
        ei=0;
        ei(i)=h;
        eiF=th+ei;
        eiB=th-ei;
        grad(i)= (FUNCTOR(eiF)-FUNCTOR(eiB) ) / (2*h) ;
    }


    return grad;
}

// numerical hessian
dbmat HESScdif(double (*fun) (const dbcolvec& ), const dbcolvec& th)
{
    int dim=th.nr();
    dbmat hess(dim,dim);
    dbcolvec ei(dim);
    dbcolvec ej(dim);
    double fval=fun(th);
    for ( int i = 0; i < dim; ++i)
    {
        ei=0;
        ei(i)=h;
        for ( int j = 0; j < dim; ++j)
        {
            ej = 0;
            ej(j) = h;

            if (i == j)
            {
                hess(i,i) = ( -fun(th + 2.0 * ei) + 16.0 *
                              fun(th + ei) - 30.0 * fval + 16.0 *
                              fun(th - ei) -
                              fun(th - 2.0 * ei)) / (12.0 * h2);
            }
            else
            {
                hess(i,j) = ( fun(th + ei + ej) - fun(th + ei - ej)
                              - fun(th - ei + ej) + fun(th - ei - ej))
                            / (4.0 * h2);
            }
        }
    }


    return hess;
}

// sample from multivariate normal with mean mu and covarance matrix Sig
// output: Out
void mvrnorm(dlib::rand& rnd, const int& sampleSize,const int& dim, const dbcolvec& mu, const dbmat& Sig, dbmat& Out)
{
    dbmat ZM(dim,sampleSize);
    dbmat muM(dim,sampleSize);
    for(int j=0; j < sampleSize; ++j )
    {
        set_subm( muM,range(0,dim-1), range(j,j) ) = mu;
        for(int i=0; i<dim; ++i)
        ZM(i,j)=rnd.get_random_gaussian();
    }

    dbmat L(dim,dim);
    L=chol(Sig);
    Out=trans(muM+L*ZM);
}

// perform complete case action action;
// DatOut initially has row DatinSz
// on output: the first NumComCases rows of DatOut contains cases without NA in any Var
void NAomit(const dbmat& DatIn, const int& DatinSz, const int& dim,\
    const dbmat& MisInd, dbmat& DatOut, int& NumComCases)
{
    NumComCases=0;
    for(int i=0; i<DatinSz; ++i)
    {
        if( sum(subm(MisInd,range(i,i),range(0,dim-1) ) ) <0.5) //
        {
            set_subm( DatOut,range(NumComCases,NumComCases),range(0,dim-1) ) \
                =subm(DatIn,range(i,i),range(0,dim-1) );
            NumComCases=NumComCases+1;
        }
    }
}


// perform complete case action action; also compute missing rates
void NAomit_MisRates(const dbmat& Dat, const dbmat& MisInd, dbmat& DatOut, int& NumComCases, dbcolvec& MissingRates)
{
    NumComCases = 0;  int Sz = Dat.nr(), dim = Dat.nc(),j; MissingRates = 0;
    for(int i=0; i<Sz; ++i)
    {
        if( sum(subm(MisInd,range(i,i),range(0,dim-1) ) ) <0.5) // no NA
        {
        set_subm( DatOut,range(NumComCases,NumComCases),range(0,dim-1) )  \
    =subm(Dat,range(i,i),range(0,dim-1) );
            NumComCases = NumComCases+1;
        }

      for( j=0; j<dim; ++j )
       if( MisInd(i,j) > 0.2 ) MissingRates(j) = MissingRates(j) +1;
    }

      MissingRates = MissingRates/( double(Sz) );
}





double expit(const double& in)
{
    return std::exp(in)/ (  1+std::exp(in) );
}

double logit(const double& in)
{
    return std::log(  in/(1-in) );
}

// generate missing data indicator:
// MCAR:  Pr(rij=1)=Pr(xij is missing)=expit(btj); btj=logit(Misratej)
void GenMCARMisInd(dlib::rand& rnd, const int&Sz, \
    const int& dim, const dbcolvec& misrates, dbmat& MisInd)
{
    double temp;
    double misratej;
    for(int j=0; j<dim; ++j)
    {
        misratej=misrates(j);

        for(int i=0; i<Sz; ++i)
        {
            temp=rnd.get_random_double(); // generate uniform(0,1)
            if(   temp<misratej  )  MisInd(i,j)=1;
            // xij is missing with probability misratej
            else MisInd(i,j)=0;
        }
    }
}

  //  xi1,....,  xi(p-1) are all MCAR
  //  MAR:  Pr(rip=1)=Pr(xip is missing)=expit(bta0+bta1*xi1*(1-ri1))
  //   xi1*(1-ri1) is always observable
  // bta1 is fixed at 0.1; chose bta0 such that E ( Pr(rip=1) )= misrate_p
void GenMARMisInd(dlib::rand& rnd, const int&Sz, \
    const int& dim, const dbcolvec& misrates,const double& beta0,\
   const double& beta1, const dbmat& Dat, const double& mu1, dbmat& MisInd)
{
    double temp,lincom,prob;
    double misratej;

    for(int i=0; i<Sz; ++i)
    {
        for(int j=0; j<(dim-1); ++j)
        {
         misratej=misrates(j);
            temp=rnd.get_random_double(); // generate uniform(0,1)
            if(   temp<misratej  )  MisInd(i,j)=1;
            // xij is missing with probability misratej
            else MisInd(i,j)=0;
        }
        lincom=beta0+beta1*Dat(i,0)*( 1-MisInd(i,0) );
        prob=expit(lincom);
        temp=rnd.get_random_double();
        if(temp<prob) MisInd(i,(dim-1) )=1;
         else  MisInd(i,(dim-1) )=0;

    }
}


  //  xi1,  xi(p-1) are all MCAR
  //  MAR:  Pr(rip=1)=Pr(xip is missing)=expit(bt0+bt1*xi1*(1-ri1))
  //   xi1*(1-ri1) is always observable
// for a given missing rate,  decide beta so that
// E ( expit(bt0+bt1*xi1*(1-ri1)) )= misrates_p
//  simulate N1=8000 iid xi1 from N(mu1, 1 );
//  simulate  N1=8000  ri1 from Pr(ri1=1)=misrates_1
//  average over N1  expit(bt0+bt1*xi1*(1-ri1))
// // loop  for each choice of bt0 from -6 to 0 by=0.001 and b1=0.1
// In R: check expit=function(x) exp(x)/( 1+exp(x) ); curve(expit,from=-6,to=0)
// ### Call this before simulation loop
void MARRatesToBeta(dlib::rand& rnd,const double& mu1,\
  const dbcolvec& misrates,const double& b0cadStepSz, \
   const double& b0halfrange, const int& SimuNumxi1, \
    const double& beta1,double& beta0, double& ActualRate)
{
  int Numb0Cad; Numb0Cad=int(2*b0halfrange/b0cadStepSz)+1;
  // std::cout<<"Numb0Cad: "<<Numb0Cad<<"\n";
  dbcolvec xi1Sam(SimuNumxi1),ri1Sam(SimuNumxi1),linTerm(SimuNumxi1);
  int k,i;
  double temp, b0Cad;
  for(k=0;k<SimuNumxi1;++k)
  {
    xi1Sam(k) = mu1+rnd.get_random_gaussian();
    temp = rnd.get_random_double();
    if( temp < misrates(0) ) ri1Sam(k)=1;
  }

  double   probmisi; ActualRate=1.0;
  // int whichmin=0;
  for(i=0;i<Numb0Cad;++i)
  {
    b0Cad=-b0halfrange+i*b0cadStepSz;
    probmisi=0;
    for(k=0;k<SimuNumxi1;++k)
      { probmisi=probmisi+expit( b0Cad+beta1*xi1Sam(k)*( 1-ri1Sam(k) ) );
      // linTerm(k)=beta1*xi1Sam(k)*( 1-ri1Sam(k) );
      }
      probmisi=probmisi/SimuNumxi1;
    if( std::abs( probmisi-misrates(1) ) < std::abs( ActualRate-misrates(1) )  )
      {  ActualRate=probmisi;  beta0=b0Cad;
      // whichmin=i;
      }

  }
  // std::cout<<"i at: "<<whichmin<<"\n";
  //std::cout<<"max linterm "<<max(linTerm)<<"\n"<<"min linterm "<<min(linTerm)<<"\n";

}


// xi1 rates r1;
// xi2 rates r2,  r2= E( logit{beta02+beta12*xi2}  ) w.r.t xi2~N(mu2,sig2^2) sig2=1;
void MNARRatesToBeta0vd2(dlib::rand& rnd,const double& mu2, const dbcolvec& misrates,\
const double& beta1, double& b0, double& b02,  double& ActualRate)
{
 int SimuNumxi1=8000;  double b0cadStepSz=0.005; double b0halfrange=6;
  b0=std::log( misrates(0)/( 1-misrates(0) ) ) ;  //logit{r1}=beta0;
 int Numb0Cad; Numb0Cad=int(2*b0halfrange/b0cadStepSz)+1;
  double xijSam;
  int k,i;
  double linTerm,b0Cad,ProbMissb0;
  ActualRate=1;

    for(i=0;i<Numb0Cad;++i)
    {
      b0Cad=-b0halfrange+i*b0cadStepSz;
      ProbMissb0=0;
      for(k=0;k<SimuNumxi1;++k)
      {
        xijSam=mu2+rnd.get_random_gaussian();
        linTerm=b0Cad+beta1*xijSam;
        ProbMissb0=ProbMissb0+1/( 1+ std::exp( -1.0*linTerm) );
      }
      ProbMissb0=ProbMissb0/SimuNumxi1;

      if( std::abs( ProbMissb0-misrates(1) ) < \
        std::abs(ActualRate-misrates(1) )  )
      {  ActualRate=ProbMissb0;  b02=b0Cad;   }
    }

}


// xi1 is MCAR:  logit{ Pr(xi1 is missing) }=beta0;
// xi2 is MNAR:  logit{ Pr(xi2 is missing) }=beta02+beta1*xi2;
void GenMNARMisIndd2(dlib::rand& rnd, const int&Sz, const double& beta1,\
 const double& b0, const double& b02, const dbmat& Dat, dbmat& MisInd)
{
    double temp,lincom,prob;
    for(int i=0; i<Sz; ++i)
    {
        for(int j=0;j<2;++j)
        {
        if(j==0) lincom=b0;
        else  lincom=b02 + beta1*Dat(i,1);
        prob=1/( 1+ std::exp( -lincom ) );
        temp=rnd.get_random_double();
        if(temp<prob) MisInd(i,j )=1;
        else  MisInd(i,j )=0;
        }
    }
}


// input: Dat; row i for obs i (1*dim);
// output: sample mean vector mu; sample Covariance matrix Sig
void SamMeanCov(const dbmat& Dat, const int& SamSz, const int& dim,\
    dbcolvec& mu, dbmat& Sig)
{
    for(int i=0; i<dim; ++i) mu(i)=mean( subm(Dat,range(0,SamSz-1),range(i,i)) );
    dbmat All1(SamSz,1);
    All1=1;
    dbmat SamMeanExt(dim,SamSz);
    SamMeanExt=All1*trans(mu);
    Sig=trans(Dat-SamMeanExt)*(Dat-SamMeanExt) /SamSz;
}


dbmat* pMQP;
// functions required for ChibarTtsch2
double QPObj ( const dbcolvec& arg)
{
    return trans(arg)*(*pMQP)*arg;
}

dbcolvec QPObjGrad ( const dbcolvec& arg)
{
    return 2*(*pMQP)*arg;
}

int NegaComp(dbcolvec& arg, const int& dim)
{
    int  Cts=0;
    for(int k=0; k<dim; ++k) if( arg(k)<0 )  Cts=Cts+1;
    return Cts;
}  // test: const int dim=3; dbcolvec z1(dim); z1=-1.1,-0.2,-1e-20; cout<<NegaComp(z1,dim)<<endl;

//############ Important notes:
// Wts are NOT AFFECTED if V is multiplied by a constant; In substitution method, V should be R_ThetaZeroInv_Rt_Est = Sz*CovMuHat; We can use CovMuHat
//############
// input
// output: Wts (of size dim+1) entry i is  wi_pV_minusO (chapter 2)
/*
compute wi_pV_minusO= Pr{ Proj(z,V,O^{-p}) has exactly i negative components,  z~N(0,V) }; where Proj(z,V,O^{-p})=argmin (z-th)^T V^{-1} (z-th) over th \in O^{-p}
 simulation-based calculation:
# (a) generate z~N(0,V) (V=L*trans(L); z1~N(0,I_p), z=L*z1)
# (b) find Proj(z,V,O^{-r})
# (c) repeat (a)-(b) N times, wi_pV_minusO is approximated by ni/N, where ni is the #times  Proj(z,V,O^{-r})  in (b)  has exactly i negative components.

In (b)  set z-th=b;   min trans(b)*invV*b s.t.  b>=z;  solution is bhat; then  Proj(z,V,O^{-r})=z-bhat
*/
// to speed up, use zmat as a paramter (no memeory-allocation)
void ChiBarWtsDimGt2ch2_BOBYQA(const dbmat& V,  dbmat& InvV, const int& dim,\
  const int& SimSz, dlib::rand& rndCBwts, dbcolvec& Wts,int& WtsFail,\
  const int&MaxFunEvalLong)
{

    dbcolvec mu(dim);
    mu=0;

   // mvrnorm(rndCBwts,SimSz,dim,mu,V,ChSam);
    //using namespace std; ofstream out;
   // out.open("mvsample.txt");  out<<ChSam<<endl;  out.close();

    pMQP=&InvV;  long int npt = 2*dim + 1;  int k,j;

    dbcolvec z(dim), z1(dim), proj(dim), ub(dim);
    ub=BIGNUM;
    dbcolvec b(dim);  dbmat L(dim,dim); L=chol(V);

    // matrix<int,0,1> projNegComp(Sz);  projNegComp=0;
    matrix<int,0,1> FreqNegComp(dim+1);
    FreqNegComp=0;
    // ofstream outNegComp; outNegComp.open("NegComp.txt");
    for(k=0; k<SimSz; ++k)
    {

        for(j = 0; j<dim; ++j) { z1(j) = rndCBwts.get_random_gaussian(); } //z1~N(0,I_dim)
          z = L*z1;    // z~N(0,V)
        b=z; // initial value of th is 0

     //find_min_box_constrained(lbfgs_search_strategy(10),
              // objective_delta_stop_strategy(1e-9),
                    //  QPObj, QPObjGrad, b, z,ub );
        find_min_bobyqaConvgFail(QPObj,b,npt,z,ub,10,1e-8,  // stopping trust region radius
            MaxFunEvalLong, // max number of objective function evaluations
            WtsFail);
        proj=z-b;
        // projNegComp(k)=NegaComp(proj,dim);
        FreqNegComp( NegaComp(proj,dim) ) ++;
    }
    // outNegComp<<projNegComp<<endl; outNegComp.close();

    for(k=0; k<(dim+1); ++k)
    {
        Wts(k)=double(FreqNegComp(k))/SimSz;
    }
    pMQP=0;  // RESET global data pointer

}



void ChiBarWtsDimGt2ch2_BFGS(const dbmat& V,  dbmat& InvV, const int& dim,\
  const int& SimSz, dlib::rand& rndCBwts, dbcolvec& Wts,int&WtsFail,int& MaxItera)
{

    dbcolvec mu(dim);   mu=0;

    // mvrnorm(rndCBwts,SimSz,dim,mu,V,ChSam);
   // using namespace std; ofstream out;
  //  out.open("mvsample.txt");

    pMQP=&InvV;

    dbcolvec z(dim), z1(dim), proj(dim), ub(dim);
    ub=BIGNUM;  dbmat L(dim,dim); L=chol(V);
    dbcolvec b(dim);

    // matrix<int,0,1> projNegComp(Sz);  projNegComp=0;
    matrix<int,0,1> FreqNegComp(dim+1);
    FreqNegComp=0;  int k,j;
    // ofstream outNegComp; outNegComp.open("NegComp.txt");
    for(k=0; k<SimSz; ++k)
    {
        for(j = 0; j<dim; ++j) { z1(j) = rndCBwts.get_random_gaussian(); } //z1~N(0,I_dim)
          z = L*z1;    // z~N(0,V)
        b=z; // initial value of th is 0
 //out<<trans(z);
        find_min_box_constrainedConvFail(bfgs_search_strategy(),
                                 objective_delta_stop_strategy(1e-5),
                                 QPObj, QPObjGrad, b, z,ub,WtsFail,MaxItera );
             if( WtsFail == 1)  { return; }
        proj=z-b;
        // projNegComp(k)=NegaComp(proj,dim);
        FreqNegComp( NegaComp(proj,dim) ) ++;
    }
    // outNegComp<<projNegComp<<endl; outNegComp.close();

    for(k=0; k<(dim+1); ++k)
    {
        Wts(k)=double(FreqNegComp(k))/SimSz;
    }
    pMQP=0;  // RESET global data pointer

}



/* when dim=2
w0_pV_minusO=  acos( V[1,2]/sqrt(V[1,1]*V[2,2]) ) /(2*pi);
w1_pV_minusO=0.5; w2_pV_minusO=0.5-w0_pV_minusO
*/
//############ Important notes:
// Wts are NOT AFFECTED if V is multiplied by a constant; In substitution method, V should be R_ThetaZeroInv_Rt_Est = Sz*CovMuHat; We can use CovMuHat
//############
void ChiBarWtsDim2ch2(const dbmat& V,dbcolvec& Wts)
{
    Wts(1)=0.5;
    Wts(0)=std::acos( V(0,1)/std::sqrt( V(0,0)*V(1,1) ) ) /(2*dlib::pi);
    Wts(2)=0.5-Wts(0);
}

// return Pr(ChBar>=x) given x and dim/ weights determing ChBar
// Wts is of length (dim+1)
// Chi-Bar-Upper Tail Prob (pval): cdf= \sum_i^{i=0}{p} w(p-i)_pV_minusO * Pr( \chi_i<=c )= \sum_i^{i=1}{p} w(p-i)_pV_minusO * Pr( \chi_i<=c ) +  wp_pV_minusO; since  Pr( \chi_{0}<=c )=1 for any c>=0
// so cdf=sum_{k=0}^{dim-1} Wts(k)*pchisq(x, double(dim-k) ) + Wts(dim)
double ChiBarPvalue(const int&dim, const dbcolvec& Wts, const double& x)
{

    double cbcdfatx;
    cbcdfatx=Wts(dim);
    for(int k=1; k<=dim; ++k)
        cbcdfatx = cbcdfatx+ Wts(dim-k)*scythe::pchisq(x, double(k));
    return 1-cbcdfatx;
}


/*
 Eta is dim*(dim+1)/2 by 1; L is dim*dim
 correspondence between Eta and L (L is such that Sig=L*t(L) ):
  Eta(a_s,0) corresponds to  L(s,s)   then a_s - a_{s-1}= p+1-s; a0=0;
 a_s=dim*s+s-s(s+1)/2;  for example, dim=3; a1=3
 //input: Eta is dim*(dim+1)/2 by 1;
//output: L is dim*dim;
 */
void EtaToL( const dbcolvec& Eta, const int& dim, dbmat& L)
{
    L=0;
    for(int s=0; s<dim; ++s)
    {
        int a_s=dim*s+s-s*(s+1)/2;
        L(s,s)=exp(Eta(a_s));
        if(s<dim-1)
            set_subm(L,range(s+1,dim-1), range(s,s)) = subm(Eta,range(a_s+1,a_s+dim-s-1), range(0,0));
        // L(s+1,s,dim-1,s)=Eta(a_s+1,0,a_s+dim-s-1,0);
        // note a_{s+1} - a_s= dim-s
    }
}

//L to Eta
void LToEta(const dbmat& L, const int& dim, dbcolvec& Eta )
{
    for(int s=0; s<dim; ++s)
    {
        int a_s=dim*s+s-s*(s+1)/2;
        Eta(a_s)=log( L(s,s));
        if(s<dim-1)
            set_subm(Eta,range(a_s+1,a_s+dim-s-1), range(0,0))=subm(L,range(s+1,dim-1), range(s,s));
    }
}

// compute multivariate normal pdf given det(Sig), inv(Sig), mu and a single observation x
// mu and x are dim*1 column vectors
double logDmvnormD(const double& detSig, const dbcolvec& mu,\
                   const dbmat& InvSig, const int& dim, const dbcolvec& x)
{
    double pdf;
    if(dim==1)
    {
        pdf=-0.5*log(detSig)-0.5* ( InvSig(0,0) ) * (x(0)- mu(0) )* (x(0)-mu(0) );
    }
    else
    {
        pdf=-0.5*log(detSig)-0.5*( trans(x-mu)*InvSig*(x-mu) )(0,0);
    }
    return pdf;
}  // test:  int dim=3; dbmat Sig(dim,dim); Sig=0.2; for(int i=0;i<dim;++i) Sig(i,i)=1;
//  dbmat InvSig(dim,dim); InvSig=inv(Sig); dbcolvec x(dim), mu(dim); mu=1.1,2.2,3.3;  x=1,1,1; double detSig=det(Sig); cout<<logDmvnorm(detSig,mu,InvSig,dim,x);


int MLkOSamComCtr=0;
dbmat* pDat;
// minus log-likelihood of one-sample complete data using dlib LA utilities (less stable than scythe)
// ASSUME pDat already pointed to the dataset
double MLkOSamCom(const dbcolvec& beta)
{
    MLkOSamComCtr=MLkOSamComCtr+1;
    int Sz = (*pDat).nr();
    int dim=(*pDat).nc();
    int SzEta = dim*(dim+1)/2;

    dbcolvec mu(dim);
    mu=subm(beta,range(0,dim-1),range(0,0));

    dbcolvec Eta(SzEta);
    Eta=subm(beta,range(dim,dim+SzEta-1),range(0,0));

    dbmat L(dim,dim);
    EtaToL(Eta,dim,L);

    dbmat InvSig(dim,dim);
    dbmat InvL=inv(L);
    InvSig=trans(InvL)*InvL;
    double detL=det(L);
    double detSig=detL*detL;

    // the loglikelihood
    double loglike = 0.0;
    dbcolvec x(dim);
    x=0;

    for (int i=0; i<Sz; ++i)
    {
        x=trans( subm( (*pDat),range(i,i),range(0,dim-1) ) );
        //std::cout<<" row "<<i<<": "<<trans(x)<<std::endl;
        // std::cout<<"loglike from row: "<<i<<" ="<<logDmvnorm(detSig,mu,InvSig,dim, x)<<std::endl;
        loglike += logDmvnormD(detSig,mu,InvSig,dim, x);
        //  std::cout<<"loglike=: "<<loglike<<std::endl;
    }

    return -1.0 * loglike;
}

dbcolvec MLkOSamComGrad(const dbcolvec& beta)
{
    return gradcdif(MLkOSamCom,beta);
}

dbmat MLkOSamComHessian(const dbcolvec& beta)
{
    dbmat res=HESScdif(MLkOSamCom, beta);
    return res;
}


// GetRi:  let  x(i,_) be the i-th row in the dataset of size 1*dim
// Ri is a  NumObs*dim matrix such that Ri*t( x(i,_) ) is x_{obs,i}
void GetRi(dbmat& Ri, const int&dim, const int&NumObs, const dbcolvec& MisInd)
{

    Ri=0;
    if(NumObs==dim)
    {
        for(int j=0; j<dim; ++j) Ri(j,j)=1;
    }
    else
    {
        int ObsCount=0;
        for(int j=0; j<dim; ++j)
        {
            if( abs( MisInd(j) )<0.01 )
            {
                Ri(ObsCount,j)=1;
                ObsCount=ObsCount+1;
            }
        }
    }

}


 dbmat* PMisInd;  int MLkOSamMARCtr;
  // typedef  std::vector<dbmat> dbmatVec;
  // typedef  std::vector<dbcolvec> dbcolvecVec;

 double MLkOSamMAR(const dbcolvec& theta)
 // ASSUMED that pDat pointed to the dataset and
 // pMisInd pointed to the mis-indicator matrix

 // first loop i=1,Sz obtain how max(NumObs row i) over i;
 //for each possible NumObs, pre-allocate space Sigi, mui,Ri,Li  (using vector of matrices)
 // second loop i=1,Sz  adds log-like from obs_i

  {
    MLkOSamMARCtr=MLkOSamMARCtr+1;
    int Sz=(*pDat).nr();
    int dim=(*pDat).nc();
    int dimEta=dim*(dim+1)/2;
    dbcolvec x(dim),MisIndxi(dim);
    dbcolvec muFull(dim);
    dbcolvec EtaFull(dimEta);
    dbmat LFull(dim,dim),SigFull(dim,dim);

    muFull=subm(theta,range(0,dim-1),range(0,0) );
    EtaFull=subm(theta,range(dim,dim+dimEta-1),range(0,0) );
    EtaToL(EtaFull,dim,LFull);  SigFull=LFull*trans(LFull);

 // if(std::abs( max(SigFull) )>1e6 || std::abs( min(SigFull) )<1e-6 )
// std::cout<<"Sig Reconstructed:\n"<<SigFull<<"\n";

    int i;

    double detLi,detSigi,pdf;
    double loglike = 0.0;
    for ( i=0; i<Sz; ++i)
    {
        MisIndxi=trans( subm( (*PMisInd) ,range(i,i),range(0,dim-1)) );

        int NumObs=dim-int( sum(MisIndxi) ) ;

        if(NumObs==0) continue; // no observed data from row i

        dbmat Ri(NumObs,dim);
        GetRi(Ri,dim,NumObs,MisIndxi);
        dbcolvec mui(NumObs);   mui=Ri*muFull;
        dbmat Sigi(NumObs,NumObs);
        Sigi=Ri*SigFull*trans(Ri);

        dbcolvec xi(NumObs);
        xi=Ri*trans( ( subm((*pDat),range(i,i),range(0,dim-1) )  ) );
        //std::cout<<"observed data in row "<<i<<": "<<trans(xi)<<"\n";
        //std::cout<<"indicator in row "<<i<<": "<<trans(MisIndxi)<<"\n";

        if(NumObs==1)  // if NumObs==1, mui, Sigi, xi are scalars
        {
        pdf=-0.5*log( Sigi(0,0) )-0.5 * (xi(0)- mui(0) )* (xi(0)-mui(0) ) / ( Sigi(0,0) );
        }
        else  // if NumObs>1, mui, xi are colvec and Sigi is matrix
        {
        dbmat Li(NumObs,NumObs);  dbmat InvLi(NumObs,NumObs);
        dbmat InvSigi(NumObs,NumObs);
        Li=chol(Sigi);
        InvLi=inv(Li);
        InvSigi=trans(InvLi)*InvLi;

            detLi=det(Li);
        detSigi=detLi*detLi;

        pdf= -0.5*log(detSigi)-0.5*( trans(xi-mui)*InvSigi*(xi-mui) )(0,0);
        }
       // std::cout<<"pdf from row "<<i<<" "<<pdf<<"\n\n";
        loglike =loglike+pdf;
    }
    return -1.0 * loglike;


  }


dbcolvec MLkOSamMARGrad(const dbcolvec& beta)
 {
 return gradcdif(MLkOSamMAR,beta);
 }

dbmat MLkOSamMARHessian (const dbcolvec& beta)
{
  dbmat res=HESScdif(MLkOSamMAR, beta);
  return res;
}


class D2MAR_model
{
public:
    typedef ::dbcolvec column_vector;
    typedef matrix<double> general_matrix;
    double operator() (
        const column_vector& x
    ) const { return MLkOSamMAR(x); }

    void get_derivative_and_hessian (
        const column_vector& x,
        column_vector& der,
        general_matrix& hess
    ) const
    {
        der = gradcdif(MLkOSamMAR,x);
        hess = HESScdif(MLkOSamMAR,x);
    }
};




//============ for MNAR dim=2
void GetQHPts(std::ifstream& QHfile, int& m, dbmat& QHpts)
  {
    int i,j;
    for(i=0;i<m;++i)
    {
     for(j=0;j<2;++j)
    QHfile>>QHpts(i,j);
    }
  }


dbmat* pGHqpts;
//########## Approximate the loglikehood via GH

// evaluate logged \int ( 1/( sqrt(2pi)*sig) )*(  1/ (1+exp(-beta0-beta1*x) ) exp(-(x-mu)^2/(2*sig*sig) dx   using GH
// theta: mu,sig1,beta0,beta1,  mu depends on i; sig=exp(sig1)
// BEFORE CALLING, ensure pGHqpts is ready
void LogExpectedMissProbUnivGH(const double& mu, const double& sig,\
    const double& beta0, const double& beta1, double& out)
    {
    int mQHpts=(*pGHqpts).nr(); out=0;
    for(int s=0; s<mQHpts; ++s)
    {
    out=out + ( (*pGHqpts)(s,1) ) / ( 1+std::exp( -beta0 \
    -beta1* ( std::sqrt(2) *sig* (*pGHqpts)(s,0) +mu  )  ) );
    }
    out=std::log( out)- 0.5*std::log( dlib::pi );
    }

// compute \int f(ri|xi) f(xi) dxi, where xi~N(mu,Sig) and f(ri2|xi2)= f(ri2|xi2)
//  f(ri1|xi1)=1/( 1+exp(-beta0) )
// using  GQ
// theta= mu, Eta, betav
// BEFORE CALLING, ensure pGHqpts is ready
void LogExpectedMissProbBivd2GH(const dbcolvec& mu, const dbmat& L, \
 const double& beta0, const double& beta02, const double& beta1,  double&out)
 {
   int mQHpts=(*pGHqpts).nr();
// std::cout<<"mQHpts="<<mQHpts<<"\n";
   dbcolvec zj(2),qj(2); double wj1,wj2;
   out=0;
   for(int j1=0;j1<mQHpts;++j1)
    {
        zj(0)=std::sqrt(2)*(*pGHqpts)(j1,0);
        wj1=(*pGHqpts)(j1,1);

        for(int j2=0;j2<mQHpts;++j2)
        {
        zj(1)=std::sqrt(2)*(*pGHqpts)(j2,0);
        wj2=(*pGHqpts)(j2,1);
        qj=mu+L*zj;
        out=out + (wj1*wj2)/ (  ( 1+std::exp( -beta0 )  ) \
        *( 1+std::exp( -beta02-beta1*qj(1) )  )   );
        }
    }
    out=std::log(out)-std::log(dlib::pi);
}



// ensure pDat   pMisInd  pGHqpts are all ready
// xi1 is MCAR, xi2 is MNAR with logit{ Pr(xi2 is missing) }= beta02+beta1*xi2
// theta: mu, Eta, beta0,beta02,beta1
int D2MLkOSamMNARCtrd2GH=0;
double D2MLkOSamMNARd2mGH(const dbcolvec& theta)
{
D2MLkOSamMNARCtrd2GH=D2MLkOSamMNARCtrd2GH+1;
int Sz=(*pDat).nr();   int dim=(*pDat).nc(),dimEta=dim*(dim+1)/2;
    dbmat L(dim,dim),Sig(dim,dim),Linv(dim,dim); dbcolvec Eta(dimEta);
    Eta=subm(theta,range(dim,dim+dimEta-1),range(0,0));
    EtaToL(Eta,dim,L);  Sig=L*trans(L);  Linv=inv(L);  double detL=det(L);
    dbcolvec mu(dim); mu=subm(theta,range(0,dim-1),range(0,0) );

//if( std::abs( max(Sig))>1e8 || std::abs( min(Sig))<1e-8  )
//std::cout<<"reconstructed Sig: \n"<<Sig<<"\n";

    double sig1=std::sqrt( Sig(0,0) ),sig2=std::sqrt( Sig(1,1) ),\
    rou=Sig(0,1)/std::sqrt( Sig(0,0)*Sig(1,1)  );

    double beta0,beta02,beta1 ;
    beta0=theta(5); beta02=theta(6); beta1=theta(7);

    int i;  double logExpectMissProbBiv,logExpectMissProbUni;

    int IfExpectedMissProbUnivGHRequired=0;
    //set IfExpectedMissProbUnivGHRequired=1
    // if there is an i such that xi1,xi2 are both missing
// int printedBivExp=0;

    double loglike = 0.0,loglikei; dbcolvec xi(2);
    double meanCond,SigCond;
    for ( i=0; i<Sz; ++i)
    {
        if( std::abs( (*PMisInd)(i,0) )<0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1,xi2 observed;  bivariate normal pdf use choleskey factor for Sig
        { xi(0)=(*pDat)(i,0);  xi(1)=(*pDat)(i,1);
        loglikei=-std::log(detL)-0.5*( trans(xi-mu)*trans(Linv)*Linv*(xi-mu) )(0,0)\
        -std::log( 1+std::exp(beta0 ) )\
        -std::log( 1+std::exp(beta02+beta1*xi(1) ) );
        }

        else if( std::abs( (*PMisInd)(i,0) )>0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1 missing, xi2 observed
        { loglikei=-std::log( 1+std::exp(beta02+beta1*(*pDat)(i,1) ) ) \
        -std::log( 1+std::exp(-beta0) ) \
        -std::log(sig2)-0.5*( xi(1)-mu(1) )*( xi(1)-mu(1) )/(sig2*sig2);
        }

        else if( std::abs( (*PMisInd)(i,0) )<0.2\
            && std::abs( (*PMisInd)(i,1) )>0.2 )
        //xi1 observed, xi2 missing
        {   // xi2 conditional on xi1
        meanCond=mu(1)+rou*(sig2/sig1)*( (*pDat)(i,0)-mu(0) );
        SigCond= sig2*sig2*(1-rou*rou) ;
        LogExpectedMissProbUnivGH(meanCond,SigCond,beta02,beta1,logExpectMissProbUni);

//std::cout<<"ExpectMissProbUni due to mis xi2 at row "<<i<<" :"<<std::exp(logExpectMissProbUni)<<"\n";
// std::cout<<trans(Condi_musigbetav)<<"\n";
        loglikei=-std::log(sig1)-0.5* ( (*pDat)(i,0)-mu(0) )\
         * ( (*pDat)(i,0)-mu(0) )/(sig1*sig1) \
        +logExpectMissProbUni  -std::log( 1+std::exp(beta0 ) );
        }

        else
        //xi1 missing, xi2 missing
        {
            if(IfExpectedMissProbUnivGHRequired==0)
            {
            LogExpectedMissProbBivd2GH(mu,L,beta0,beta02,beta1,logExpectMissProbBiv);
            IfExpectedMissProbUnivGHRequired=1;
//if(printedBivExp==0) {std::cout<<" ExpectMissProbBiv used:\n"; printedBivExp=1;}
            }

        loglikei=logExpectMissProbBiv;
        }
      // std::cout<<"loglikei: "<<loglikei<<"\n";

        loglike =loglike+loglikei;
    }

  return (-1.0)*loglike;
}


dbcolvec D2MLkOSamMNARd2mGHGrad(const dbcolvec& theta)
{
  return  gradcdif(D2MLkOSamMNARd2mGH,theta);
  	}

dbmat D2MLkOSamMNARd2mGHHess(const dbcolvec& theta)
{
  return  HESScdif(D2MLkOSamMNARd2mGH,theta);
  	}


class D2MNARind2_model
{
public:
    typedef ::dbcolvec column_vector;
    typedef matrix<double> general_matrix;
    double operator() (
        const column_vector& x
    ) const { return D2MLkOSamMNARd2mGH(x); }

    void get_derivative_and_hessian (
        const column_vector& x,
        column_vector& der,
        general_matrix& hess
    ) const
    {
        der = gradcdif(D2MLkOSamMNARd2mGH,x);
        hess = HESScdif(D2MLkOSamMNARd2mGH,x);
    }
};






void BootstrapCovMuHatMNARd2(dbmat& CovMuHat,int & NumSkippedBootSample, \
const int& BootRepNum, const dbcolvec& thetaEst,\
 int& SamSz,dlib::rand& rnd, dbmat& AQHPts,const dbcolvec& thetaTrue,\
 const int& MaxFunEval,const double& StopCriteria)
// draw BootRepNum bootstrap samples from the model indexed by theta to estimate CovMuHat, where MuHat is the MLE of mu under H0 (constrained MLE)
// NumSkippedBootSample returns the number of bootstrap samples dropped due to numerical non-convergence
//
// thetaTrue are true parameter values used as starting values
{
    int dim=2; int dimEta=dim*(dim+1)/2;int dimTheta=dim+dimEta+dim+1;
    long int npt=2*dimTheta+1;

    double b0Hat, b02Hat, beta1Hat;
    dbcolvec MuHat(dim); dbmat SigHat(dim,dim);
    dbcolvec EtaHat(dimEta); dbmat LHat(dim,dim);

    MuHat=subm(thetaEst,range(0,dim-1),range(0,0));
    EtaHat=subm(thetaEst,range(dim,dim+dimEta-1),range(0,0));
    EtaToL(EtaHat,dim,LHat); SigHat=LHat*trans(LHat);
    b0Hat=thetaEst(5);  b02Hat=thetaEst(6); beta1Hat=thetaEst(7);

    dbmat SumMuhatRep(dim,dim); SumMuhatRep=0;
    dbmat SumMuhatMuhatTransRep(dim,dim); SumMuhatMuhatTransRep=0;

    dbcolvec muhatBoot(dim),thetaBoot(dimTheta);
    dbcolvec thetaBootUb(dimTheta),thetaBootLb(dimTheta);

    dbmat DatBoot(SamSz,dim); dbmat MisIndBoot(SamSz,dim);

  int bootDone=0;
    while(bootDone<BootRepNum)
     {
     // generate dat
      mvrnorm(rnd,SamSz,dim,MuHat,SigHat,DatBoot);

     //generate pMisInd
     GenMNARMisIndd2(rnd,SamSz,beta1Hat,b0Hat,b02Hat, DatBoot, MisIndBoot);

     // compute restricted MLE
     pDat = &DatBoot; PMisInd = &MisIndBoot;   pGHqpts = &AQHPts;
     int ConvgFail=0;  thetaBootUb=BIGNUM; thetaBootLb=-BIGNUM;
     thetaBootUb(0)=0; thetaBootUb(1)=0;  thetaBoot=thetaTrue;
     double minObjec= find_min_bobyqaConvgFail(D2MLkOSamMNARd2mGH,
            thetaBoot, npt,    // number of interpolation points
            thetaBootLb,  // lower bound constraint
            thetaBootUb,   // upper bound constraint
            5,    // initial trust region radius
            StopCriteria,  // stopping trust region radius
            MaxFunEval, // max number of objective function evaluations
            ConvgFail
            );

    if(ConvgFail==1 || !is_finite(minObjec) )  { ++NumSkippedBootSample; continue;  }

    muhatBoot=subm(thetaBoot,range(0,1),range(0,0));

    SumMuhatRep=SumMuhatRep + muhatBoot;
    SumMuhatMuhatTransRep=SumMuhatMuhatTransRep + muhatBoot*trans(muhatBoot);
    ++bootDone;
     }

    CovMuHat=SumMuhatMuhatTransRep/BootRepNum-SumMuhatRep*trans(SumMuhatRep)/(BootRepNum*BootRepNum);

}


bool IsWeird_d2MNAR(const dbcolvec& theta, const double& MinNegLogLike,\
 const double&BigPos, const double&BigNeg )
{
    bool outcome;
    if( !is_finite(MinNegLogLike) ) { outcome=true; return outcome;  }

    dbcolvec Eta(3),betaVec(3); dbmat L(2,2);
    int dim=2; Eta=subm(theta,range(2,4),range(0,0) );
    EtaToL(Eta,dim,L); L=L*trans(L);

    betaVec=subm(theta,range(5,7),range(0,0) );

    if( max(L)> 1e2 || min(betaVec)<BigNeg || max(betaVec)>BigPos )
         { outcome=true; return outcome;  }
    else { outcome=false;  return outcome;  }

}


// calculate Tw based on the full dataset, using true Cov(MuHat) in the calculation of Tw and p-value;
// idea: the instances with numerical non-convergence in the calculation of constrained or unconstrained are NOT systematically different, that is, the proportion of RejH0 is also roughly 5% among thoseinstances
void D2_FullDatTwResult(dbcolvec& SamMeanFull, const dbmat&Shat, dbmat& ShatInv,\
 double&Tw, double& StopCriteria, bool& RejH0)
{

     int dim=2;
    dbcolvec ub(2),lb(2);   lb=SamMeanFull; ub=BIGNUM;
 pMQP=&ShatInv;

    Tw=find_min_box_constrained(lbfgs_search_strategy(10),\
    objective_delta_stop_strategy(StopCriteria), \
    QPObj, QPObjGrad, SamMeanFull, lb, ub);

    dbcolvec CbWts(dim+1);  CbWts=0;
    ChiBarWtsDim2ch2(Shat,CbWts); double pval=ChiBarPvalue(dim,CbWts,Tw);
    if( pval<SigLevel ) RejH0=true;
    else RejH0=false;
}



// convert thetahat (mu, eta) to mu sig
void MuEta_to_mu_sig(const int&dim, dbcolvec& theta, dbcolvec& mu, dbmat& Sig)
{
int dimEta = dim*(dim+1)/2;  dbcolvec Eta(dimEta);
mu = subm(theta,range(0,dim-1),range(0,0));
Eta = subm(theta,range(dim,dim+dimEta-1),range(0,0));
dbmat L(dim,dim);   EtaToL(Eta,dim,L);  Sig=L*trans(L);
}

// compute Tw
void get_Tw(double&Tw, const int&dim, dbcolvec& muHat, \
 dbmat& InvCovMuHat,const double&StopCriteria)
{
// Tw= min trans(muhat-proj)*inv(CovMuHat)*(muhat-proj)
// proj =argmin (trans(muhat-x)*inv(CovMuHat)*(muhat-x)) over x\in O_minus
// set  th=muhat-x; proj =argmin trans(th)*inv(CovMuHat)*th  over th >=muhat
 pMQP=&InvCovMuHat;  dbcolvec TwLb(dim),TwUb(dim);

 TwUb=BIGNUM;  TwLb=muHat;
Tw=find_min_box_constrained(lbfgs_search_strategy(10),\
        objective_delta_stop_strategy(StopCriteria), \
        QPObj, QPObjGrad, muHat, TwLb,TwUb);
}


//
dbmat *pQHPts; dbcolvec* pXdat; dbcolvec* pRdat;
  double UniNormMNARMLogLikGH(const dbcolvec& theta)
   {
    int Sz=(*pXdat).nr(),i,s;

   // std::cout<<"Sz: "<<Sz<<"\n";
    double mu,sig,beta0,beta1;
    mu=theta(0); sig=std::exp( theta(1) ); beta0=theta(2); beta1=theta(3);
    double logLike=0; double temp;
    int mQHpts=(*pQHPts).nr();

    double  LogLikeMis=0;

    for(s=0; s<mQHpts; ++s)
    {
    LogLikeMis=LogLikeMis + ( (*pQHPts)(s,1) ) / ( 1+std::exp( -beta0 \
    -beta1* ( std::sqrt(2) *sig* (*pQHPts)(s,0) +mu  )  ) );
    }

    LogLikeMis=std::log( LogLikeMis)- 0.5*std::log( dlib::pi );

    for( i=0; i<Sz; ++i)
    {
      if( std::abs( (*pRdat)(i)-1 ) <0.05 )  // xi is missing
      {
      temp=LogLikeMis;
      // std::cout<<"missing indicator of obs "<<i<<": "<<(*pRdat)(i)<<"\n";
      // std::cout<<"log-like obs "<<i<<": "<<temp<<"\n";
      }
      else  // xi is observed
      {
        // std::cout<<"Pri at row "<<i<<": "<<Pri<<"\n";
      temp=-theta(1) -0.5* ( 1/(sig*sig) )* ( (*pXdat)(i)-mu )*( (*pXdat)(i)-mu )  \
      - std::log( 1+std::exp(beta0+beta1*(*pXdat)(i) ) );
      //std::cout<<"missing indicator of obs "<<i<<": "<<(*pRdat)(i)<<"\n";
      //std::cout<<"log-like obs "<<i<<": "<<temp<<"\n";
      }
      logLike = logLike + temp;
    }

     return -1.0*logLike;
    }


 dbcolvec UniNormMNARMLogLikGradGH(const dbcolvec& beta)
 {
 return gradcdif(UniNormMNARMLogLikGH,beta);
 }

dbmat UniNormMNARMLogLikHessianGH(const dbcolvec& beta)
{
  return HESScdif(UniNormMNARMLogLikGH, beta);
}


//####################
// model 1:
// logit{ xi1 }= alpha01 + alpha11*xi1; logit{ xi2 }= alpha02 + alpha12*xi2;
//####################

// compute log{ \int f(rij|xij)f(xij) dxij }
//  \int f(ri1=1|xi1) f(xi1|xi2)  dxi1  xi1  missing, xi2 observed
//  \int f(ri2=1|xi2) f(xi2|xi1)  dxi2  xi2  missing, xi1 observed
//==== Ensure pDat, PMisInd,  pQHPts are ready
void LogExpectedProbMisXij_mod1(const double& mu,const double& sig,const double& alpha0,\
const double& alpha1,double& LogExpectedProbMisXij,const double&meanXij,const double& sdXij)
{
    int mQHpts=(*pQHPts).nr();   LogExpectedProbMisXij=0;

    for(int s=0; s<mQHpts; ++s)
    {
    LogExpectedProbMisXij=LogExpectedProbMisXij + ( (*pQHPts)(s,1) ) / ( 1+std::exp( -alpha0 \
    - alpha1* ( std::sqrt(2) *sig* (*pQHPts)(s,0) +mu - meanXij  )/sdXij  ) );
    }
    LogExpectedProbMisXij=std::log( LogExpectedProbMisXij );
}

// logit{ xi1 }= alpha01 + alpha11*xi1; logit{ xi2 }= alpha02 + alpha12*xi2;
// logExpectMissProbBiv= log{ \int f(ri1,ri2|xi1,xi2)f(xi1,xi2) dxi1 dxi2}
//==== Ensure pDat, PMisInd,  pQHPts are ready
void LogExpectedProbMisXi1ProbMisXi2_mod1(const dbcolvec& mu, const dbmat& L, \
 const double& alpha01, const double& alpha02, const double& alpha11, const double& alpha12,\
double&out,\
double& meanXi1, double&sdXi1, double& meanXi2, double&sdXi2)
{
   int mQHpts = (*pQHPts).nr();
   // std::cout<<"mQHpts="<<mQHpts<<"\n";
   dbcolvec zj(2),qj(2); double wj1,wj2;
   out = 0;
   for(int j1=0;j1<mQHpts;++j1)
    {
        zj(0) = std::sqrt(2)*(*pQHPts)(j1,0);
        wj1 = (*pQHPts)(j1,1);

        for(int j2=0;j2<mQHpts;++j2)
        {
        zj(1) = std::sqrt(2)*(*pQHPts)(j2,0);
        wj2 = (*pQHPts)(j2,1);
        qj = mu + L*zj;
        out = out + (wj1*wj2)/ (  ( 1+std::exp( -alpha01-alpha11*(qj(0)-meanXi1)/sdXi1 )  ) \
        *( 1+std::exp( -alpha02-alpha12*(qj(1) - meanXi1)/sdXi2 )  )   );
        }
    }
    out=std::log(out);  // -std::log(dlib::pi);
}


// log-likelihood for missing model 1:
// logit{ xi1 }= alpha01 + alpha11*xi1; logit{ xi2 }= alpha02 + alpha12*xi2;
//     alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
//==== Ensure pDat, PMisInd,  pQHPts are ready
double mod_NoRi1_MNARd2_MLogLikGH(const dbcolvec& theta)
{
    int Sz = (*pDat).nr();   int dim = 2,dimEta=dim*(dim+1)/2;
    dbmat L(dim,dim),Sig(dim,dim),Linv(dim,dim); dbcolvec Eta(dimEta);
    Eta = subm(theta,range(dim,dim+dimEta-1),range(0,0));
    EtaToL(Eta,dim,L);  Sig = L*trans(L);  Linv = inv(L);  double detL = det(L);
    dbcolvec mu(dim); mu = subm(theta,range(0,dim-1),range(0,0) );

//if( std::abs( max(Sig))>1e8 || std::abs( min(Sig))<1e-8  )
//std::cout<<"reconstructed Sig: \n"<<Sig<<"\n";

    double var1 =   Sig(0,0), var2 =  Sig(1,1), cov12 = Sig(0,1);   // rou = Sig(0,1)/std::sqrt( Sig(0,0)*Sig(1,1)  );
    double sd1 = std::sqrt(var1), sd2 = std::sqrt(var2);
    double mean1 = mu(0), mean2 = mu(1);

    double alpha01,alpha11,alpha02,alpha12 ;
    alpha01 = theta(5); alpha11 = theta(6); alpha02 = theta(7); alpha12 = theta(8);

    int i;  double logExpectMissProbBiv,logExpectMissProbUni;

    int IfExpectedMissProbUnivGHRequired = 0;
    //set IfExpectedMissProbUnivGHRequired=1
    // if there is an i such that xi1,xi2 are both missing
// int printedBivExp=0;

    double loglike = 0.0,loglikei; dbcolvec xi(2);
    double meanCond,SigCond;
    for ( i=0; i<Sz; ++i)
    {
        if( std::abs( (*PMisInd)(i,0) )<0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1,xi2 observed;  bivariate normal pdf use choleskey factor for Sig
        { xi(0) = (*pDat)(i,0);  xi(1) = (*pDat)(i,1);
        loglikei = -std::log(detL)-0.5*( trans(xi-mu)*trans(Linv)*Linv*(xi-mu) )(0,0)\
        -std::log( 1+std::exp( alpha01 + alpha11*(xi(0)-mean1)/sd1 ) )\
        -std::log( 1+std::exp( alpha02 + alpha12*(xi(1)-mean2)/sd2 ) );
        }

        else if( std::abs( (*PMisInd)(i,0) )>0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1 missing, xi2 observed
        {
            // xi1 conditional on xi2
        meanCond = mu(0)+ cov12* ( (*pDat)(i,1)-mu(1) )/var2;
        SigCond = std::sqrt(var1 - cov12*cov12/var2);
        LogExpectedProbMisXij_mod1(meanCond,SigCond,alpha01, alpha11,logExpectMissProbUni,mean1,sd1);

        loglikei = -0.5*std::log(var2)-0.5* ( (*pDat)(i,1)-mu(1) )\
         * ( (*pDat)(i,1)-mu(1) )/var2 \
        +logExpectMissProbUni  -std::log( 1+std::exp( alpha02 + alpha12*(xi(1)-mean2)/sd2 ) );
        }

        else if( std::abs( (*PMisInd)(i,0) )<0.2\
            && std::abs( (*PMisInd)(i,1) )>0.2 )
        //xi1 observed, xi2 missing
        {   // xi2 conditional on xi1
        meanCond = mu(1)+ cov12* ( (*pDat)(i,0)-mu(0) )/var1;
        SigCond = std::sqrt(var2- cov12*cov12/var1);
LogExpectedProbMisXij_mod1(meanCond,SigCond,alpha02, alpha12,logExpectMissProbUni,mean2,sd2);
//std::cout<<"ExpectMissProbUni due to mis xi2 at row "<<i<<" :"<<std::exp(logExpectMissProbUni)<<"\n";
// std::cout<<trans(Condi_musigbetav)<<"\n";
        loglikei=-0.5*std::log(var1)-0.5* ( (*pDat)(i,0)-mu(0) )\
         * ( (*pDat)(i,0)-mu(0) )/var1 \
        +logExpectMissProbUni  -std::log( 1+std::exp( alpha01 + alpha11*(xi(0)-mean1)/sd1 ) );
        }

        else
        //xi1 missing, xi2 missing
        {
            if(IfExpectedMissProbUnivGHRequired==0)
            {
LogExpectedProbMisXi1ProbMisXi2_mod1(mu, L,alpha01, alpha02, alpha11, alpha12,logExpectMissProbBiv,\
                                     mean1,sd1,mean2,sd2);
            IfExpectedMissProbUnivGHRequired = 1;
//if(printedBivExp==0) {std::cout<<" ExpectMissProbBiv used:\n"; printedBivExp=1;}
            }

        loglikei=logExpectMissProbBiv;
        }
      // std::cout<<"loglikei: "<<loglikei<<"\n";

        loglike =loglike+loglikei;
    }

  return (-1.0)*loglike;

}


dbcolvec mod_NoRi1_MNARd2_MLogLikGH_Grad(const dbcolvec& theta)
{
  return  gradcdif(mod_NoRi1_MNARd2_MLogLikGH,theta);
}

dbmat mod_NoRi1_MNARd2_MLogLikGH_Hess(const dbcolvec& theta)
{
  return  HESScdif(mod_NoRi1_MNARd2_MLogLikGH,theta);
  	}


class mod_NoRi1_MNARd2_MLogLikGH_class
{
public:
    typedef ::dbcolvec column_vector;
    typedef matrix<double> general_matrix;
    double operator() (
        const column_vector& x
    ) const { return mod_NoRi1_MNARd2_MLogLikGH(x); }

    void get_derivative_and_hessian (
        const column_vector& x,
        column_vector& der,
        general_matrix& hess
    ) const
    {
        der = gradcdif(mod_NoRi1_MNARd2_MLogLikGH,x);
        hess = HESScdif(mod_NoRi1_MNARd2_MLogLikGH,x);
    }
};


//####################
// model 2:
// logit{ prob xi1 mis }= alpha01 + alpha11*(xi1-mu1)/sd1;
// logit{ prob xi2 mis }= alpha02 + alpha12*(xi2-mu2)/sd2; + alpha22 * ri1;
//####################

// compute log{ \int f(ri1|xi1)f(xi1) dxi1 }
//==== Ensure pDat, PMisInd,  pQHPts are ready
void LogExpectedProbMisXi1_mod2(const double& mu,const double& sig,const double& alpha01,\
const double& alpha11,double& LogExpectedProbMisXi1,double&muXi1,double&sdXi1)
{
    int mQHpts = (*pQHPts).nr();   LogExpectedProbMisXi1= 0;

    for(int s=0; s<mQHpts; ++s)
    {
    LogExpectedProbMisXi1 = LogExpectedProbMisXi1 + ( (*pQHPts)(s,1) ) / ( 1+std::exp( -alpha01 \
    - alpha11* ( std::sqrt(2) *sig* (*pQHPts)(s,0) +mu -  muXi1)/sdXi1  ) );
    }
    LogExpectedProbMisXi1 = std::log( LogExpectedProbMisXi1 );
}

// compute log{ \int f(ri2|ri1,xi1)f(xi1) dxi1 }; logit{ xi2 }= alpha02 + alpha12*xi2 + alpha22 * ri1;  xi2 is missing xi1 is observed
//==== Ensure pDat, PMisInd,  pQHPts are ready
void LogExpectedProbMisXi2_mod2(const double& mu,const double& sig,const double& alpha02,\
const double& alpha12, double& LogExpectedProbMisXi2,\
double&muXi2,double&sdXi2)
{
    int mQHpts=(*pQHPts).nr();   LogExpectedProbMisXi2=0;

    for(int s=0; s<mQHpts; ++s)
    {
    LogExpectedProbMisXi2 = LogExpectedProbMisXi2 + ( (*pQHPts)(s,1) ) / ( 1+std::exp( -alpha02 \
    - alpha12 * ( std::sqrt(2) *sig* (*pQHPts)(s,0) +mu - muXi2  )/sdXi2   ) );
    }
    LogExpectedProbMisXi2 = std::log( LogExpectedProbMisXi2 );
}


// logit{ prob xi1 mis }= alpha01 + alpha11*xi1;  logit{ prob xi2 mis }= alpha02 + alpha12*xi2 + alpha22 * ri1;
// logExpectMissProbBiv= log{ \int f(ri1,ri2|xi1,xi2)f(xi1,xi2) dxi1 dxi2}
//==== Ensure pDat, PMisInd,  pQHPts are ready
void LogExpectedProbMisXi1MisXi2_mod2(const dbcolvec& mu, const dbmat& L, \
 const double& alpha01, const double& alpha02, const double& alpha11, const double& alpha12,\
const double& alpha22, double&out,\
double&mu1, double&sd1, double&mu2, double&sd2)
{
int mQHpts = (*pQHPts).nr();
// std::cout<<"mQHpts="<<mQHpts<<"\n";
   dbcolvec zj(2),qj(2); double wj1,wj2;
   out = 0;
   for(int j1=0;j1<mQHpts;++j1)
    {
        zj(0)=std::sqrt(2)*(*pQHPts)(j1,0);
        wj1=(*pQHPts)(j1,1);

        for(int j2=0;j2<mQHpts;++j2)
        {
        zj(1) = std::sqrt(2)*(*pQHPts)(j2,0);
        wj2 = (*pQHPts)(j2,1);
        qj = mu+L*zj;
        out = out + (wj1*wj2)/ (  ( 1+std::exp( -alpha01-alpha11*(qj(0)-mu1)/sd1 )  ) \
        *( 1+std::exp( -alpha02 -alpha12*(qj(1)-mu2)/sd2  - alpha22  ) )   );
        }
    }
    out = std::log(out) ; // - std::log(dlib::pi);
}


// log-likelihood for missing model 2:
// logit{ prob xi1 mis }= alpha01 + alpha11*xi1;  logit{ prob xi2 mis }= alpha02 + alpha12*xi2 + alpha22 * ri1;
//==== Ensure pDat, PMisInd,  pQHPts are ready
double mod2_MNARd2_MLogLikGH(const dbcolvec& theta)
{
    int Sz=(*pDat).nr();   int dim=2,dimEta=dim*(dim+1)/2;
    dbmat L(dim,dim),Sig(dim,dim),Linv(dim,dim); dbcolvec Eta(dimEta);
    Eta=subm(theta,range(dim,dim+dimEta-1),range(0,0));
    EtaToL(Eta,dim,L);  Sig=L*trans(L);  Linv=inv(L);  double detL=det(L);
    dbcolvec mu(dim); mu=subm(theta,range(0,dim-1),range(0,0) );

//if( std::abs( max(Sig))>1e8 || std::abs( min(Sig))<1e-8  )
//std::cout<<"reconstructed Sig: \n"<<Sig<<"\n";

     double var1 =   Sig(0,0), var2 =  Sig(1,1), cov12 = Sig(0,1);
     double mu1=mu(0), mu2=mu(1), sd1=std::sqrt(var1), sd2 = std::sqrt(var2);

    double alpha01,alpha11,alpha02,alpha12,alpha22 ;
    alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
    alpha22=theta(9);

    int i;  double logExpectMissProbBiv,logExpectMissProbUni;

    int IfExpectedMissProbUnivGHRequired = 0;
    //set IfExpectedMissProbUnivGHRequired=1
    // if there is an i such that xi1,xi2 are both missing
// int printedBivExp=0;

    double loglike = 0.0, loglikei; dbcolvec xi(2);
    double meanCond,SigCond;
    for ( i=0; i<Sz; ++i)
    {
        if( std::abs( (*PMisInd)(i,0) )<0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1,xi2 observed;  bivariate normal pdf use choleskey factor for Sig
        { xi(0) = (*pDat)(i,0);  xi(1) = (*pDat)(i,1);
        loglikei = -std::log(detL)-0.5*( trans(xi-mu)*trans(Linv)*Linv*(xi-mu) )(0,0)\
        -std::log( 1+std::exp( alpha01 + alpha11*(xi(0)-mu1)/sd1 ) )\
        -std::log( 1+std::exp( alpha02 + alpha12*(xi(1)-mu2)/sd2   ) );
        }

        else if( std::abs( (*PMisInd)(i,0) )>0.2 \
        && std::abs( (*PMisInd)(i,1) )<0.2 )
        //xi1 missing, xi2 observed
        {
            // xi1 conditional on xi2
         meanCond = mu(0) + cov12* ( (*pDat)(i,1)-mu(1) )/var2;
         SigCond = std::sqrt(var1 - cov12*cov12/var2);
LogExpectedProbMisXi1_mod2(meanCond,SigCond,alpha01, alpha11,logExpectMissProbUni,mu1,sd1);
        loglikei = -0.5*std::log(var2)-0.5* ( (*pDat)(i,1)-mu(1) )\
         * ( (*pDat)(i,1)-mu(1) )/var2 \
        +logExpectMissProbUni  - std::log( 1+std::exp( alpha02 + alpha12*(xi(1)-mu2)/sd2 + alpha22 ) );
        }

        else if( std::abs( (*PMisInd)(i,0) )<0.2\
            && std::abs( (*PMisInd)(i,1) )>0.2 )
        //xi1 observed, xi2 missing
        {   // xi2 conditional on xi1
        meanCond = mu(1) + cov12* ( (*pDat)(i,0)-mu(0) )/var1;
        SigCond = std::sqrt(var2- cov12*cov12/var1);

LogExpectedProbMisXi2_mod2(meanCond,SigCond,alpha02, alpha12,logExpectMissProbUni,mu2,sd2);
//std::cout<<"ExpectMissProbUni due to mis xi2 at row "<<i<<" :"<<std::exp(logExpectMissProbUni)<<"\n";
// std::cout<<trans(Condi_musigbetav)<<"\n";
        loglikei = -0.5*std::log(var1)-0.5* ( (*pDat)(i,0)-mu(0) )\
         * ( (*pDat)(i,0)-mu(0) )/var1\
        +logExpectMissProbUni  -std::log( 1+std::exp( alpha01 + alpha11*(xi(0)-mu1)/sd1 ) );
        }

        else
        //xi1 missing, xi2 missing
        {
            if(IfExpectedMissProbUnivGHRequired==0)
            {
          LogExpectedProbMisXi1MisXi2_mod2(mu, L,alpha01, alpha02, alpha11, \
                                            alpha12,alpha22,logExpectMissProbBiv,mu1,sd1,mu2,sd2);
            IfExpectedMissProbUnivGHRequired = 1;
//if(printedBivExp==0) {std::cout<<" ExpectMissProbBiv used:\n"; printedBivExp=1;}
            }

        loglikei = logExpectMissProbBiv;
        }
      // std::cout<<"loglikei: "<<loglikei<<"\n";

        loglike =loglike+loglikei;
    }

  return (-1.0)*loglike;

}


dbcolvec mod2_MNARd2_MLogLikGH_Grad(const dbcolvec& theta)
{
  return  gradcdif(mod2_MNARd2_MLogLikGH,theta);
}

dbmat mod2_MNARd2_MLogLikGH_Hess(const dbcolvec& theta)
{
  return  HESScdif(mod2_MNARd2_MLogLikGH,theta);
  	}



class mod2_MNARd2_MLogLikGH_class
{
public:
    typedef ::dbcolvec column_vector;
    typedef matrix<double> general_matrix;
    double operator() (
        const column_vector& x
    ) const { return mod2_MNARd2_MLogLikGH(x); }

    void get_derivative_and_hessian (
        const column_vector& x,
        column_vector& der,
        general_matrix& hess
    ) const
    {
        der = gradcdif(mod2_MNARd2_MLogLikGH,x);
        hess = HESScdif(mod2_MNARd2_MLogLikGH,x);
    }
};


  void  ZinExMod_GenerateMisInd(dlib::rand& rndMisInd,const dbmat& Dat,\
                                const dbcolvec& alphaVec,dbmat& MisInd,const dbcolvec& mu,\
                                const double& alphaMAR,const double&sd1,const double&sd2)
  {
    int Sz = Dat.nr();  double misProb,RandUniform;
    for(int i = 0; i < Sz; ++i)
    {
    misProb = 1/(1 + std::exp( -alphaVec(0) - alphaVec(1)*(Dat(i,0) -mu(0) )/sd1 ) );
    // xi1
    RandUniform = rndMisInd.get_random_double();
    if( RandUniform < misProb)  MisInd(i,0) = 1;
    else  MisInd(i,0) = 0;

    misProb = 1/(1 + std::exp( -alphaVec(2) - alphaVec(3) * (Dat(i,1) -mu(1) )/sd2 \
                              -alphaMAR*( Dat(i,0)-mu(0) )*(1-MisInd(i,0) ) /sd1 )  );
                              // alphaMAR*Standardized_x_{i1}*(1-r_{i1})
    // xi2
    RandUniform = rndMisInd.get_random_double();
    if( RandUniform < misProb)  MisInd(i,1) = 1;
    else  MisInd(i,1) = 0;
    }
  }


void GetMinMax(const dbcolvec& Input, double& Min, double& Max)
{
   int len = Input.nr();  Min = Input(0); Max = Input(0);
   for(int k=0; k<len; ++k)
   {
       if( Input(k) < Min) Min = Input(k);
       else if(Input(k) > Max) Max = Input(k);
   }
}

bool is_positive_definite_2by2(dbmat& M)
{
  int dim = M.nr();
  if(dim!=2) {std::cout<<"check_positive_definite_2by2 called by a matrix NOT of size 2*2\n"; return false;}
  else
  {
    if(M(0,0)>0 && M(1,1)>0 && M(0,0)*M(1,1)>M(0,1)*M(1,0) \
       && std::abs(M(0,1)-M(1,0))<1e-3*std::abs(M(1,0)) )
    return true;
    else return false;
  }
}


#endif
