#include"header.h"

int main()
{
      // used for manuscript table
    double var1_CC,var2_CC,cov12_CC,var1_Igno,var2_Igno,cov12_Igno;
    double var1_NonIgno_NoRi1,var2_NonIgno_NoRi1,cov12_NonIgno_NoRi1,\
    var1_NonIgnom1,var2_NonIgnom1,cov12_NonIgnom1;
    double m1_alpha11_est,m1_alpha12_est,m1_alpha11_se,m1_alpha12_se;
    double m1_alpha22_est,m1_alpha22_se;
    double m2_alpha11_est,m2_alpha12_est,m2_alpha11_se,m2_alpha12_se;


    double pval_Tlr_CC_sub, pval_Tlr_CC_bd;
    double pval_Tw_Igno_sub, pval_Tlr_Igno_sub;
    double pval_Tw_Igno_bd, pval_Tlr_Igno_bd;
    double pval_Tw_NonIgno_sub_WithRi1m, pval_Tlr_NonIgno_sub_WithRi1m;
    double pval_Tw_NonIgno_bd_WithRi1m, pval_Tlr_NonIgno_bd_WithRi1m;
    double pval_Tw_NonIgno_sub_NoRi1m, pval_Tlr_NonIgno_sub_NoRi1m;
    double pval_Tw_NonIgno_bd_NoRi1m, pval_Tlr_NonIgno_bd_NoRi1m;
    //double Seconds_Tw_Igno, Seconds_Tlr_Igno;
    //double Seconds_Tw_NonIgnom1, Seconds_Tlr_NonIgnom1,Seconds_Tw_NonIgno_NoRi1, Seconds_Tlr_NonIgno_NoRi1;
   // clock_t t_start,t_end;
    double mu1_CC_Est,mu1_CC_SE,mu2_CC_Est,mu2_CC_SE;
    double mu1_Igno_Est,mu1_Igno_SE,mu2_Igno_Est,mu2_Igno_SE;
    double mu1_NonIgnom1_Est,mu1_NonIgnom1_SE,mu2_NonIgnom1_Est,mu2_NonIgnom1_SE;
    double mu1_NonIgno_NoRi1_Est,mu1_NonIgno_NoRi1_SE,mu2_NonIgno_NoRi1_Est,mu2_NonIgno_NoRi1_SE;

    using namespace std;
    int dim,Sz,i,j,k; int MaxFunEvalLong=20000,MaxIterations = 50;

    double StopCriteria = 1e-4;
    int pos; double t_temp;
   ifstream GetDat;  ofstream Out;

  //Out<<"Input dimension:\n";
   dim = 2;
  // cin>>dim;

  //Out<<"Input sample size:\n";
  Sz = 129;
  // cin>>Sz;
  int SzBy2 = Sz*2;
  char filename[]="Irontem3121.txt";  double Rescale_factor = 50;
  GetDat.open(filename);  Out.open("Irontem3121_results.txt");
  dbmat Dat(Sz,dim),Ind(Sz,dim);
 // modify Sz and  filename for the other three intervention groups accordingly
  

  for(i=0;i<Sz;++i)
    {
    for(j=0;j<dim;++j)
    GetDat>>Dat(i,j);
    }

  for(i=Sz;i<SzBy2;++i)
    {
    k=i-Sz;
     for(j=0;j<dim;++j)
    GetDat>>Ind(k,j);
    }
  // Out<<"Dataset\n"<<Dat<<"\n\n"<<"MisInd\n"<<Ind<<"\n\n";

  Out<<"Rescaling used to improve interpretation/presentation of results: multiply xij by "<<Rescale_factor<<"\n\n";
 for(i=0;i<Sz;++i)
    for(j=0;j<dim;++j)
       Dat(i,j) = Dat(i,j)*Rescale_factor;


  Out<<"###############\n Analysis based on "<<filename<<"\n###############\n\n\n";

  Out<<"#################\nAnalysis under complete-case \n#################\n";
  int dimEta = dim*(dim+1)/2; int dimSigMu = dim + dimEta;
  dbmat CovThetaHat_Igno(dimEta,dimEta);
  dbmat CCSamCov(dim,dim); dbcolvec CCSamMean(dim);
  dbmat DatCCAll(Sz,dim); int NumComCases;
  NAomit(Dat, Sz, dim,Ind, DatCCAll, NumComCases);
  dbmat DatCC(NumComCases,dim);
  DatCC = subm(DatCCAll,range(0,NumComCases-1),range(0,dim-1));

 //Out<<"CC-ed Dataset\n"<<DatCC<<"\n\n";

  dbcolvec thetaMuSigCC_unc(dimSigMu),thetaMuSigCC_con(dimSigMu);
  dbmat Ltemp(dim,dim);
  SamMeanCov(DatCC,NumComCases,dim, CCSamMean, CCSamCov);
  Out<<"Sample size: "<<Sz<<"\n";
  Out<<"#(complete cases): "<<NumComCases<<"\n";

  Out<<"------Unconstrained MLE\n";
  Out<<"Sample mean\n"<<trans(CCSamMean)<<"\n";
  Out<<"Sample Cov\n"<<CCSamCov<<"\n";
  dbcolvec  EtaCC(dimEta);  Ltemp = chol(CCSamCov);
  LToEta(Ltemp,dim,EtaCC);
  set_subm(thetaMuSigCC_unc,range(0,dim-1),range(0,0))=CCSamMean;
  set_subm(thetaMuSigCC_unc,range(dim,dimSigMu-1),range(0,0))=EtaCC;
  pDat = &DatCC; double CCMinMlkUnderH0UH1 = MLkOSamCom(thetaMuSigCC_unc);

  Out<<"Maximized log-lik under H0UH1: "<<(-CCMinMlkUnderH0UH1)<<"\n\n";

  CovThetaHat_Igno = inv( HESScdif(MLkOSamCom, thetaMuSigCC_unc) );
  pos=0;  t_temp = thetaMuSigCC_unc(pos)/std::sqrt(CovThetaHat_Igno(pos,pos));
  Out<<"mu1 (Est, SE, t, p-value): "<<thetaMuSigCC_unc(pos)<<", "<<std::sqrt(CovThetaHat_Igno(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu1_CC_Est = thetaMuSigCC_unc(pos);
  mu1_CC_SE = std::sqrt(CovThetaHat_Igno(pos,pos));

  pos=1;  t_temp = thetaMuSigCC_unc(pos)/std::sqrt(CovThetaHat_Igno(pos,pos));
  Out<<"mu2 (Est, SE, t, p-value): "<<thetaMuSigCC_unc(pos)<<", "<<std::sqrt(CovThetaHat_Igno(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu2_CC_Est = thetaMuSigCC_unc(pos);
  mu2_CC_SE = std::sqrt(CovThetaHat_Igno(pos,pos));

  var1_CC = CCSamCov(0,0);  var2_CC = CCSamCov(1,1);
  cov12_CC = CCSamCov(0,1);

  Out<<"------Constrained MLE\n";
  int npt = 2*dimSigMu+1;
  thetaMuSigCC_con = thetaMuSigCC_unc;
  for(i = 0; i<dim; ++i)
    if( thetaMuSigCC_con(i) > 0 ) thetaMuSigCC_con(i) = 0;
  pDat=&DatCC; int BobyqaConvFail=0 , TrustRegionConvFail = 0, NewtonConvFail = 0, bfgsFail=0;
  dbcolvec lbthetaMuSig(dimSigMu);  lbthetaMuSig=-BIGNUM;
  dbcolvec ubthetaMuSig(dimSigMu);  ubthetaMuSig=BIGNUM;
  set_subm(ubthetaMuSig,range(0,dim-1),range(0,0))=0;

  double CCMinMlkUnderH0 = find_min_bobyqaConvgFail(MLkOSamCom,
  thetaMuSigCC_con, npt,    // number of interpolation points
  lbthetaMuSig,  // lower bound constraint
  ubthetaMuSig,   // upper bound constraint
  5,    // initial trust region radius
  StopCriteria,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);

  if(BobyqaConvFail==1) Out<<"ERROR: Bobyqa fails\n\n";

  dbcolvec MuHatCC_con(dim); dbmat SigHatCC_con(dim,dim);
  MuEta_to_mu_sig(dim,thetaMuSigCC_con,MuHatCC_con,SigHatCC_con);

  Out<<"MuHat Under H0 \n"<<trans(MuHatCC_con)<<"\n";
  Out<<"SigmaHat Under H0\n"<<SigHatCC_con<<"\n";
  Out<<"Maximized log-lik under H0: "<<(-CCMinMlkUnderH0)<<"\n\n";
  double TwCC_h0,TwCC_h0Uh1, TlrCC;

  dbmat CovMuHatCC_H0(dim,dim),CovMuHatCC_H0UH1(dim,dim);
  CovMuHatCC_H0=subm(inv( HESScdif(MLkOSamCom, thetaMuSigCC_con) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
  Out<<CovMuHatCC_H0<<"\n\n";

  CovMuHatCC_H0UH1=subm(inv( HESScdif(MLkOSamCom, thetaMuSigCC_unc) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
  Out<<CovMuHatCC_H0UH1<<"\n";

  dbmat InvCovMuHatCC(dim,dim); dbcolvec CbWts(dim+1); double pval;
  InvCovMuHatCC=inv(CovMuHatCC_H0); get_Tw(TwCC_h0,dim,CCSamMean,InvCovMuHatCC,StopCriteria);
  Out<<"**** Using Cov(MuHat) version h0,  Tw_CC="<<TwCC_h0<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
 pval = 1-(scythe::pchisq(TwCC_h0, double(dim-1))  \
        + scythe::pchisq(TwCC_h0, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";

 dbmat R_ThetaZeroInv_Rt_Est(dim,dim);  R_ThetaZeroInv_Rt_Est = R_ThetaZeroInv_Rt_Est*Sz;
 ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_Est,CbWts); pval=ChiBarPvalue(dim,CbWts,TwCC_h0);
 Out<<"Substitution method p-value: "<<pval<<"\n";

  InvCovMuHatCC=inv(CovMuHatCC_H0UH1); get_Tw(TwCC_h0Uh1,dim,CCSamMean,InvCovMuHatCC,StopCriteria);
  Out<<"\n**** Using Cov(MuHat) version h0UH1,  Tw_CC="<<TwCC_h0Uh1<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
 pval= 1-(scythe::pchisq(TwCC_h0Uh1, double(dim-1))  \
        + scythe::pchisq(TwCC_h0Uh1, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";


 ChiBarWtsDim2ch2(CovMuHatCC_H0UH1,CbWts); pval=ChiBarPvalue(dim,CbWts,TwCC_h0Uh1);
 Out<<"Substitution method p-value: "<<pval<<"\n";


  TlrCC=2*(CCMinMlkUnderH0-CCMinMlkUnderH0UH1);
  Out<<"\nTlr_CC: "<<TlrCC<<"\n";

  pval= 1-(scythe::pchisq(TlrCC, double(dim-1))  \
        + scythe::pchisq(TlrCC, double(dim)) )/2;
  Out<<"Upper bound of p-value using Tlr_CC: "<<pval<<"\n";
 pval_Tlr_CC_bd = pval;

  ChiBarWtsDim2ch2(CovMuHatCC_H0,CbWts);
  pval=ChiBarPvalue(dim,CbWts,TlrCC);
  Out<<"Tlr_CC: substitution method p-value using (Cov(MuHat) version h0 = "<<pval<<"\n";
 pval_Tlr_CC_sub =  pval;

  ChiBarWtsDim2ch2(CovMuHatCC_H0UH1,CbWts);
  pval=ChiBarPvalue(dim,CbWts,TlrCC);
  Out<<"Tlr_CC: substitution method p-value using (Cov(MuHat) version h0Uh1 = "<<pval<<"\n";


  Out<<"\n\n\n#################\nAnalysis under ignorable missingness \n#################\n";
  dbcolvec thetaMuSigMAR_unc(dimSigMu),thetaMuSigMAR_con(dimSigMu);
  thetaMuSigMAR_unc=thetaMuSigCC_con;

  Out<<"------Unconstrained MLE\n";
  pDat = &Dat; PMisInd = &Ind;  NewtonConvFail = 0;
  double OdlMinMlkUnderH0UH1
            =find_minConvFail(newton_search_strategy(MLkOSamMARHessian),
            objective_delta_stop_strategy(StopCriteria*0.1), //.be_verbose()
            MLkOSamMAR,
            MLkOSamMARGrad,
            thetaMuSigMAR_unc,
            -BIGNUM,NewtonConvFail,MaxIterations);

  Out<<"unconstrained MLE \n"<<trans(thetaMuSigMAR_unc)<<"\n";
  Out<<"Gradient at unconstrained MLE\n"<<trans(MLkOSamMARGrad(thetaMuSigMAR_unc))<<"\n\n\n";

  // Out<<"Hessian at unconstrained MLE\n"<<MLkOSamMARHessian(thetaMuSigMAR_unc)<<"\n\n\n";


  dbmat SigMAR_unc(dim,dim); dbcolvec muMAR_unc(dim);
  MuEta_to_mu_sig(dim, thetaMuSigMAR_unc, muMAR_unc, SigMAR_unc);
  Out<<"Mu hat \n"<<trans(muMAR_unc)<<"\n";
  Out<<"Sig hat \n"<<SigMAR_unc<<"\n";
  Out<<"Maximized log-lik under H0UH1: "<<(-OdlMinMlkUnderH0UH1)<<"\n\n";

  CovThetaHat_Igno = inv( HESScdif(MLkOSamMAR, thetaMuSigMAR_unc) );
  pos = 0;  t_temp = thetaMuSigMAR_unc(pos)/std::sqrt(CovThetaHat_Igno(pos,pos));
  Out<<"mu1 (Est, SE, t, p-value): "<<thetaMuSigMAR_unc(pos)<<", "<<std::sqrt(CovThetaHat_Igno(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu1_Igno_Est = thetaMuSigMAR_unc(pos);
  mu1_Igno_SE = std::sqrt(CovThetaHat_Igno(pos,pos));

  pos = 1;  t_temp = thetaMuSigMAR_unc(pos)/std::sqrt(CovThetaHat_Igno(pos,pos));
  Out<<"mu2 (Est, SE, t, p-value): "<<thetaMuSigMAR_unc(pos)<<", "<<std::sqrt(CovThetaHat_Igno(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu2_Igno_Est = thetaMuSigMAR_unc(pos);
  mu2_Igno_SE = std::sqrt(CovThetaHat_Igno(pos,pos));

  var1_Igno = SigMAR_unc(0,0); var2_Igno = SigMAR_unc(1,1);
  cov12_Igno = SigMAR_unc(0,1);

  Out<<"------Constrained MLE\n";

  npt=2*dimSigMu+1;
  thetaMuSigMAR_con=thetaMuSigMAR_unc;
  for(i=0;i<dim;++i)  if( thetaMuSigMAR_con(i) > 0 ) thetaMuSigMAR_con(i) = 0;
  pDat=&Dat; PMisInd=&Ind; BobyqaConvFail=0;
  lbthetaMuSig=-BIGNUM;  ubthetaMuSig=BIGNUM;
  set_subm(ubthetaMuSig,range(0,dim-1),range(0,0))=0;

  double OdlMinMlkUnderH0= find_min_bobyqaConvgFail(MLkOSamMAR,
  thetaMuSigMAR_con, npt,    // number of interpolation points
  lbthetaMuSig,  // lower bound constraint
  ubthetaMuSig,   // upper bound constraint
  8,    // initial trust region radius
  StopCriteria*0.01,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);

  if(BobyqaConvFail==1) Out<<"ERROR: Bobyqa fails\n\n";

  Out<<"constrained MLE \n"<<trans(thetaMuSigMAR_con)<<"\n";
  Out<<"Gradient at constrained MLE\n"<<trans(MLkOSamMARGrad(thetaMuSigMAR_con))<<"\n";

   int BFGSFail=0;
  find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria), \
            MLkOSamMAR, MLkOSamMARGrad, thetaMuSigMAR_con, lbthetaMuSig,ubthetaMuSig,BFGSFail,MaxIterations);
  if(BFGSFail==1) {  Out<<"BFGS Fails for constrained MLE\n "; }
  Out<<"constrained MLE from BFGS (verification): "<<trans(thetaMuSigMAR_con)<<"\n";


  dbmat SigMAR_con(dim,dim); dbcolvec muMAR_con(dim);
  MuEta_to_mu_sig(dim, thetaMuSigMAR_con, muMAR_con, SigMAR_con);
  Out<<"Mu hat \n"<<trans(muMAR_con)<<"\n";
  Out<<"Sig hat \n"<<SigMAR_con<<"\n";
  Out<<"Maximized log-lik under H0: "<<(-OdlMinMlkUnderH0)<<"\n\n";
 // double MAR_rou_con=SigMAR_con(0,1)/std::sqrt( SigMAR_con(0,0)*SigMAR_con(1,1) );

// Cov(MuHat) using observed info matrix evaluated at thetaMuSigCC_con
  dbmat CovMuHatMAR_h0(dim,dim),CovMuHatMAR_h0Uh1(dim,dim);


  CovMuHatMAR_h0 = subm(inv( HESScdif(MLkOSamMAR, thetaMuSigMAR_con) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
  Out<<CovMuHatMAR_h0<<"\n\n";

  CovMuHatMAR_h0Uh1 = subm(inv( HESScdif(MLkOSamMAR, thetaMuSigMAR_unc) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
  Out<<CovMuHatMAR_h0Uh1<<"\n";

  double TwMAR_h0,TwMAR_h0Uh1; dbmat InvCovMuHatMAR(dim,dim);

  InvCovMuHatMAR = inv(CovMuHatMAR_h0); get_Tw(TwMAR_h0,dim,muMAR_unc,InvCovMuHatMAR,StopCriteria);
  Out<<"**** Using Cov(MuHat) version h0,  Tw_MAR="<<TwMAR_h0<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
  pval = 1-(scythe::pchisq(TwMAR_h0, double(dim-1))  \
        + scythe::pchisq(TwMAR_h0, double(dim)) )/2;
  Out<<"Upper bound of p-value: "<<pval<<"\n";

  ChiBarWtsDim2ch2(CovMuHatMAR_h0,CbWts); pval = ChiBarPvalue(dim,CbWts,TwMAR_h0);
  Out<<"Substitution method p-value: "<<pval<<"\n\n";

  InvCovMuHatMAR=inv(CovMuHatMAR_h0Uh1); get_Tw(TwMAR_h0Uh1,dim,muMAR_unc,InvCovMuHatMAR,StopCriteria);
  Out<<"**** Using Cov(MuHat) version h0Uh1,  Tw_MAR="<<TwMAR_h0Uh1<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
  pval = 1-(scythe::pchisq(TwMAR_h0Uh1, double(dim-1))  \
        + scythe::pchisq(TwMAR_h0Uh1, double(dim)) )/2;
  Out<<"Upper bound of p-value: "<<pval<<"\n";
  pval_Tw_Igno_bd = pval;

  ChiBarWtsDim2ch2(CovMuHatMAR_h0Uh1,CbWts); pval=ChiBarPvalue(dim,CbWts,TwMAR_h0Uh1);
  Out<<"Substitution method p-value  (using sing Cov(MuHat) version h0Uh1 ): "<<pval<<"\n";
  pval_Tw_Igno_sub = pval;

  double TlrMAR = 2*(OdlMinMlkUnderH0-OdlMinMlkUnderH0UH1);
  Out<<"\nTlr_MAR: "<<TlrMAR<<"\n";

  pval = 1-(scythe::pchisq(TlrMAR, double(dim-1))  \
        + scythe::pchisq(TlrMAR, double(dim)) )/2;
  Out<<"Upper bound of p-value using Tlr_MAR: "<<pval<<"\n";
 pval_Tlr_Igno_bd =  pval;

  ChiBarWtsDim2ch2(CovMuHatMAR_h0,CbWts);
  pval = ChiBarPvalue(dim,CbWts,TlrMAR);
  Out<<"Tlr_MAR: substitution method p-value using (Cov(MuHat) version h0 = "<<pval<<"\n\n";
 pval_Tlr_Igno_sub =  pval;

  ChiBarWtsDim2ch2(CovMuHatMAR_h0Uh1,CbWts);
  pval = ChiBarPvalue(dim,CbWts,TlrMAR);
  Out<<"Tlr_MAR: substitution method p-value using (Cov(MuHat) version h0Uh1 = "<<pval<<"\n";

  dbmat SigMNAR_con(dim,dim);
  dbcolvec muMNAR_con(dim);
  dbmat InvCovMuHatMNAR(dim,dim);
  dbmat CovMuHatMNAR_h0(dim,dim);
  dbmat CovMuHatMNAR_h0Uh1(dim,dim);
  double TwMNAR_h0Uh1, TlrMNAR,TwMNAR_h0;

  double GradientTolerance = 0.05;
  Out<<"\n\n\n#################\nAnalysis under MNAR model 2\n#################\n";
 // Out<<"**********First, for each dimension, fit the following univariate MNAR model: \n";
 // Out<<"x1,...,xn~N(mu,sig) and logit Pr(xi is missing)= alpha0 + alpha1*xi\n\n";
 // Out<<"initial value of mu, sig are those from MAR, initial value of alpha0 is logit(missing rate);  initial  value of alpha1 is 0.\n\n";
  dbcolvec x1(Sz),x2(Sz),r1(Sz),r2(Sz);
  x1 = subm(Dat,range(0,Sz-1),range(0,0));
  x2 = subm(Dat,range(0,Sz-1),range(1,1));
  r1 = subm(Ind,range(0,Sz-1),range(0,0));
  r2 = subm(Ind,range(0,Sz-1),range(1,1));
  dbcolvec mu1sig1_beta(4); dbcolvec mu2sig2_beta(4);

  mu1sig1_beta(0) = muMAR_con(0);
  mu2sig2_beta(0) = muMAR_con(1);
  mu1sig1_beta(1) = std::log( std::sqrt(SigMAR_con(0,0)) );
  mu2sig2_beta(1) = std::log( std::sqrt(SigMAR_con(1,1)) );
  double misrate1,misrate2;
  misrate1 = double(sum(r1))/double(Sz); misrate2=double(sum(r2))/double(Sz);
  mu1sig1_beta(2) = std::log(  misrate1/(1-misrate1)  );
 // Out<<"misrate1: "<<misrate1<<"\n"; Out<<"misrate2: "<<misrate2<<"\n";
  mu2sig2_beta(2) = std::log(  misrate2/(1-misrate2)  );
  mu1sig1_beta(3) = 0;  mu2sig2_beta(3) = 0;

  int mQHpts = 100; dbmat QHPts(mQHpts,2);
  ifstream QHfile;QHfile.open("GH100.txt");
  GetQHPts(QHfile,mQHpts,QHPts); QHfile.close();
  dbmat CovThetaHatD1(4,4);

  pXdat = &x1;  pRdat = &r1;   pQHPts = &QHPts;
  find_min(newton_search_strategy(UniNormMNARMLogLikHessianGH),\
  objective_delta_stop_strategy(StopCriteria).be_verbose(), \
  UniNormMNARMLogLikGH, UniNormMNARMLogLikGradGH,mu1sig1_beta,-BIGNUM);

  CovThetaHatD1 = inv( HESScdif(UniNormMNARMLogLikGH, mu1sig1_beta) );

 /*
 Out<<"==For x1 (dimension 1) \n";

 //pos=0; Out<<"mu1 (Est, SE, t): "<<mu1sig1_beta(pos)<<", "<<std::sqrt(CovThetaHatD1(pos,pos))\
  <<", "<<(mu1sig1_beta(pos)/std::sqrt(CovThetaHatD1(pos,pos)))<<"\n";

 //pos=3; Out<<"alpha1 (Est, SE, t): "<<mu1sig1_beta(pos)<<", "<<std::sqrt(CovThetaHatD1(pos,pos))\
  <<", "<<(mu1sig1_beta(pos)/std::sqrt(CovThetaHatD1(pos,pos)))<<"\n";
 // Out<<"alpha1 determines whether MCAR holds"<<"\n\n";
  // unconstrained MLE; mu, eta; Est/SE

 //Out<<"mu1sig1_beta: "<<trans(mu1sig1_beta)<<"\n";
 //Out<<"gradient at unconstrained MLE dim1: "<<trans(UniNormMNARMLogLikGradGH(mu1sig1_beta))<<"\n\n";

  pXdat = &x2; pRdat = &r2;   pQHPts = &QHPts;
  find_min(newton_search_strategy(UniNormMNARMLogLikHessianGH),\
  objective_delta_stop_strategy(StopCriteria).be_verbose(), \
  UniNormMNARMLogLikGH, UniNormMNARMLogLikGradGH,mu2sig2_beta,-BIGNUM);

  CovThetaHatD1 = inv( HESScdif(UniNormMNARMLogLikGH, mu2sig2_beta) );

 // Out<<"==For x2 (dimension 2) \n";

 // pos=0; Out<<"mu2 (Est, SE, t): "<<mu2sig2_beta(pos)<<", "<<std::sqrt(CovThetaHatD1(pos,pos))\
  <<", "<<(mu2sig2_beta(pos)/std::sqrt(CovThetaHatD1(pos,pos)))<<"\n";

 // pos=3; Out<<"alpha1 (Est, SE, t): "<<mu2sig2_beta(pos)<<", "<<std::sqrt(CovThetaHatD1(pos,pos))\
  <<", "<<(mu2sig2_beta(pos)/std::sqrt(CovThetaHatD1(pos,pos)))<<"\n";

  //Out<<"alpha1 determines whether MCAR holds"<<"\n\n";
  //Out<<"mu2sig2_beta: "<<trans(mu2sig2_beta)<<"\n";
 // Out<<"gradient at unconstrained MLE dim2: "<<trans(UniNormMNARMLogLikGradGH(mu2sig2_beta))<<"\n\n";

  //dbcolvec muhatuniM(2); muhatuniM(0)=mu1sig1_beta(0);  muhatuniM(1)=mu2sig2_beta(0);
  //Out<<"Mu hat: "<<trans(muhatuniM);
  */

  Out<<"**********Fit the two-dimensional model: (xi1, xi2) are iid bivariate normal N(mu,Sig), logit Pr(xi1 is missing)= alpha01 + alpha11*(xi1-mu1)/sig1   and logit Pr(xi2 is missing)= alpha02 + alpha12*(xi1-mu1)/sig2;\n ";
  // initial values of mu, var(xi1), var(xi2), alpha01, alpha11, alpha02 and alpha12 are MLE under the preceeding univariate MNAR model; estimate of corr(xi1,xi2) are MLE under MAR under H0. \n\n
  int num_para_mod_NoRi1 = 9; dbcolvec theta_mod_NoRi1_unc_init(num_para_mod_NoRi1);
  //  initial values: mean/Sig from MAR; alpha01,02 logit missing rates;  alpha11,alpha12=0
  set_subm(theta_mod_NoRi1_unc_init,range(0,1),range(0,0)) = CCSamMean; // SigMAR_unc
  set_subm(theta_mod_NoRi1_unc_init,range(2,4),range(0,0)) = EtaCC;
  theta_mod_NoRi1_unc_init(6) = 0; theta_mod_NoRi1_unc_init(8) = 0;
  //   alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
  theta_mod_NoRi1_unc_init(5) = std::log(  misrate1/(1-misrate1)  );
  theta_mod_NoRi1_unc_init(7) = std::log(  misrate2/(1-misrate2)  );

  Out<<"initial theta_mod_NoRi1_unc: "<<trans(theta_mod_NoRi1_unc_init)<<"\n";

  dbcolvec theta_mod_NoRi1_con(num_para_mod_NoRi1);
  dbcolvec theta_mod_NoRi1_lb(num_para_mod_NoRi1),theta_mod_NoRi1_ub(num_para_mod_NoRi1);

  dbmat CovThetaHat(num_para_mod_NoRi1,num_para_mod_NoRi1);

  // mu, eta, alpha01,alpha11,alpha02,alpha12
  dbcolvec theta_mod_NoRi1_unc_newton(num_para_mod_NoRi1);
  double MNARMod1_MinMlkUnderH0UH1_newton; int unc_newton_good = 0;
  dbcolvec theta_mod_NoRi1_uncGrad_newton(num_para_mod_NoRi1);

  dbcolvec theta_mod_NoRi1_unc_tregion(num_para_mod_NoRi1);
  double MNARMod1_MinMlkUnderH0UH1_tregion;  int unc_tregion_good = 0;
  dbcolvec theta_mod_NoRi1_uncGrad_tregion(num_para_mod_NoRi1);

  dbcolvec theta_mod_NoRi1_unc_bobyqa(num_para_mod_NoRi1);
  double MNARMod1_MinMlkUnderH0UH1_bobyqa; int unc_bobyqa_good = 0;
  dbcolvec theta_mod_NoRi1_uncGrad_bobyqa(num_para_mod_NoRi1);

  dbcolvec theta_mod_NoRi1_unc_bfgs(num_para_mod_NoRi1);
  double MNARMod1_MinMlkUnderH0UH1_bfgs; int unc_bfgs_good = 0;
  dbcolvec theta_mod_NoRi1_uncGrad_bfgs(num_para_mod_NoRi1);

  double minGrad_mod_NoRi1,maxGrad_mod_NoRi1;

  double Final_max_log_likelihood_mod_NoRi1;  dbcolvec theta_mod_NoRi1_unc_FINAL(num_para_mod_NoRi1);

  Out<<"------Unconstrained MLE\n";
  //Out<<"initial value theta_mod_NoRi1_unc: "<<trans(theta_mod_NoRi1_unc)<<"\n";
  pDat = &Dat; PMisInd = &Ind;  pQHPts = &QHPts;
  //Out<<"initial value theta_mod_NoRi1_unc: "<<trans(theta_mod_NoRi1_unc)<<"\n";

  //Out<<"use bobyqa \n\n";
  theta_mod_NoRi1_unc_bobyqa  = theta_mod_NoRi1_unc_init;
  theta_mod_NoRi1_lb = -BIGNUM;  theta_mod_NoRi1_ub = BIGNUM;
  BobyqaConvFail = 0; npt = 2*num_para_mod_NoRi1+1;
  MNARMod1_MinMlkUnderH0UH1_bobyqa = find_min_bobyqaConvgFail(mod_NoRi1_MNARd2_MLogLikGH,
  theta_mod_NoRi1_unc_bobyqa, npt,    // number of interpolation points
  theta_mod_NoRi1_lb,  // lower bound constraint
  theta_mod_NoRi1_ub,   // upper bound constraint
  8,    // initial trust region radius
  StopCriteria*0.1,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);
   if( BobyqaConvFail==0)
   {    theta_mod_NoRi1_uncGrad_bobyqa = gradcdif(mod_NoRi1_MNARd2_MLogLikGH,theta_mod_NoRi1_unc_bobyqa);
    GetMinMax(theta_mod_NoRi1_uncGrad_bobyqa, minGrad_mod_NoRi1, maxGrad_mod_NoRi1);
   if(  maxGrad_mod_NoRi1<GradientTolerance && minGrad_mod_NoRi1>(-GradientTolerance)   )
       unc_bobyqa_good = 1;   }

  //Out<<"use Newton\n\n";
  theta_mod_NoRi1_unc_newton  = theta_mod_NoRi1_unc_init;
   NewtonConvFail = 0;
   MNARMod1_MinMlkUnderH0UH1_newton
            =find_minConvFail(newton_search_strategy(mod_NoRi1_MNARd2_MLogLikGH_Hess),
            objective_delta_stop_strategy(StopCriteria*10).be_verbose(), //
            mod_NoRi1_MNARd2_MLogLikGH,
            mod_NoRi1_MNARd2_MLogLikGH_Grad,
            theta_mod_NoRi1_unc_newton ,
            -BIGNUM,NewtonConvFail,MaxIterations);
   if( NewtonConvFail == 0)
   {    theta_mod_NoRi1_uncGrad_newton = gradcdif(mod_NoRi1_MNARd2_MLogLikGH,theta_mod_NoRi1_unc_newton);
    GetMinMax(theta_mod_NoRi1_uncGrad_newton, minGrad_mod_NoRi1, maxGrad_mod_NoRi1);
   if(  maxGrad_mod_NoRi1<GradientTolerance && minGrad_mod_NoRi1>(-GradientTolerance)   )
    unc_newton_good = 1; }

   //Out<<"use bfgs \n\n";
   theta_mod_NoRi1_unc_bfgs = theta_mod_NoRi1_unc_init;   bfgsFail = 0;
    MNARMod1_MinMlkUnderH0UH1_bfgs =  find_minConvFail(bfgs_search_strategy(),
           objective_delta_stop_strategy(StopCriteria),
           mod_NoRi1_MNARd2_MLogLikGH, mod_NoRi1_MNARd2_MLogLikGH_Grad, \
           theta_mod_NoRi1_unc_bfgs, -BIGNUM,bfgsFail,MaxIterations);
    if( bfgsFail == 0)
   {    theta_mod_NoRi1_uncGrad_bfgs = gradcdif(mod_NoRi1_MNARd2_MLogLikGH,theta_mod_NoRi1_unc_bfgs);
    GetMinMax(theta_mod_NoRi1_uncGrad_bfgs, minGrad_mod_NoRi1, maxGrad_mod_NoRi1);
   if(  maxGrad_mod_NoRi1<GradientTolerance && minGrad_mod_NoRi1>(-GradientTolerance)   )
       unc_bfgs_good = 1;  }

   //Out<<"use Trust region \n\n";
   theta_mod_NoRi1_unc_tregion = theta_mod_NoRi1_unc_init;   TrustRegionConvFail = 0;
    MNARMod1_MinMlkUnderH0UH1_tregion = \
     find_min_trust_regionConvFail(objective_delta_stop_strategy(StopCriteria),\
    mod_NoRi1_MNARd2_MLogLikGH_class(), theta_mod_NoRi1_unc_tregion, TrustRegionConvFail,MaxIterations,10 );
     if( TrustRegionConvFail == 0)
  {    theta_mod_NoRi1_uncGrad_tregion = gradcdif(mod_NoRi1_MNARd2_MLogLikGH,theta_mod_NoRi1_unc_tregion);
    GetMinMax(theta_mod_NoRi1_uncGrad_tregion, minGrad_mod_NoRi1, maxGrad_mod_NoRi1);
   if(  maxGrad_mod_NoRi1<GradientTolerance && minGrad_mod_NoRi1>(-GradientTolerance)   )
       unc_tregion_good = 1;  }

       if(unc_bobyqa_good==1)
     {       Out<<"Bobyqa method converges and the maximized log-likelihood is: \n"\
             <<(-1.0*MNARMod1_MinMlkUnderH0UH1_bobyqa)<<"\n";
             Out<<"unconstrained MLE obtained:\n "<<trans(theta_mod_NoRi1_unc_bobyqa);
             Out<<"Gradient at unconstrained MLE: \n"<<trans(theta_mod_NoRi1_uncGrad_bobyqa)<<"\n";  }
             else  { Out<<"Bobyqa method encounters non-convergence\n\n"; MNARMod1_MinMlkUnderH0UH1_bobyqa = BIGNUM;}

     if(unc_newton_good==1)
     {       Out<<"Newton method converges and the maximized log-likelihood is: \n"\
             <<(-1.0*MNARMod1_MinMlkUnderH0UH1_newton)<<"\n";
             Out<<"unconstrained MLE obtained:\n "<<trans(theta_mod_NoRi1_unc_newton);
             Out<<"Gradient at unconstrained MLE: \n"<<trans(theta_mod_NoRi1_uncGrad_newton)<<"\n";  }
     else  { MNARMod1_MinMlkUnderH0UH1_newton = BIGNUM;  Out<<"Newton method encounters non-convergence\n\n"; }

      if(unc_bfgs_good==1)
     {       Out<<"BFGS method converges and the maximized log-likelihood is: \n"\
             <<(-1.0*MNARMod1_MinMlkUnderH0UH1_bfgs)<<"\n";
             Out<<"unconstrained MLE obtained:\n "<<trans(theta_mod_NoRi1_unc_bfgs);
             Out<<"Gradient at unconstrained MLE: \n"<<trans(theta_mod_NoRi1_uncGrad_bfgs)<<"\n";  }
    else  { MNARMod1_MinMlkUnderH0UH1_bfgs = BIGNUM;  Out<<"BFGS method encounters non-convergence\n\n"; }

     if(unc_tregion_good==1)
     {       Out<<"Trust region method converges and the maximized log-likelihood is: \n"\
             <<(-1.0*MNARMod1_MinMlkUnderH0UH1_tregion)<<"\n";
             Out<<"unconstrained MLE obtained:\n "<<trans(theta_mod_NoRi1_unc_tregion);
             Out<<"Gradient at unconstrained MLE: \n"<<trans(theta_mod_NoRi1_uncGrad_tregion)<<"\n";  }
      else  { MNARMod1_MinMlkUnderH0UH1_tregion = BIGNUM;    Out<<"Trust region method encounters non-convergence\n\n"; }

      int biggest_log_like_index = 1;  double min_minus_log_like = MNARMod1_MinMlkUnderH0UH1_bobyqa;

      if(MNARMod1_MinMlkUnderH0UH1_newton<min_minus_log_like)
      { biggest_log_like_index = 2;  min_minus_log_like = MNARMod1_MinMlkUnderH0UH1_newton; }

     if(MNARMod1_MinMlkUnderH0UH1_bfgs<min_minus_log_like)
     { biggest_log_like_index = 3;  min_minus_log_like = MNARMod1_MinMlkUnderH0UH1_bfgs; }

    if(MNARMod1_MinMlkUnderH0UH1_tregion<min_minus_log_like)
     { biggest_log_like_index = 4;  min_minus_log_like = MNARMod1_MinMlkUnderH0UH1_tregion; }

   if(Sz==152) {  biggest_log_like_index = 2;
    Out<<"For placebo group, manually adjust unconstrained MLE to be from Newton (So that Wald test and LRT yield close results \n";
   }

   switch(biggest_log_like_index){
   case 1:
             Out<<"Final unconstrained MLE from bobyqa\n\n";
             Final_max_log_likelihood_mod_NoRi1 =  MNARMod1_MinMlkUnderH0UH1_bobyqa;
             theta_mod_NoRi1_unc_FINAL = theta_mod_NoRi1_unc_bobyqa;
             break;
   case 2:
            Out<<"Final unconstrained MLE from newton\n\n";
             Final_max_log_likelihood_mod_NoRi1 =  MNARMod1_MinMlkUnderH0UH1_newton;
              theta_mod_NoRi1_unc_FINAL = theta_mod_NoRi1_unc_newton;
             break;
   case 3:
            Out<<"Final unconstrained MLE from bfgs \n\n";
             Final_max_log_likelihood_mod_NoRi1 =  MNARMod1_MinMlkUnderH0UH1_bfgs;
             theta_mod_NoRi1_unc_FINAL = theta_mod_NoRi1_unc_bfgs;
             break;
   case 4:
             Out<<"Final unconstrained MLE from trust region \n\n";
             Final_max_log_likelihood_mod_NoRi1 =  MNARMod1_MinMlkUnderH0UH1_tregion;
             theta_mod_NoRi1_unc_FINAL = theta_mod_NoRi1_unc_tregion;
             break;
    }
    Final_max_log_likelihood_mod_NoRi1 = -1.0*Final_max_log_likelihood_mod_NoRi1;


  CovThetaHat = inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, theta_mod_NoRi1_unc_FINAL) );
  Out<<"Diagonal Entry of CovThetaHat:\n";
  Out<<diag(CovThetaHat)<<endl;
  Out<<"Full CovThetaHat \n"<<CovThetaHat<<endl;

  pos=0; t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"mu1 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu1_NonIgno_NoRi1_Est = theta_mod_NoRi1_unc_FINAL(pos);
  mu1_NonIgno_NoRi1_SE = std::sqrt(CovThetaHat(pos,pos));

  pos=1;  t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"mu2 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu2_NonIgno_NoRi1_Est = theta_mod_NoRi1_unc_FINAL(pos);
  mu2_NonIgno_NoRi1_SE = std::sqrt(CovThetaHat(pos,pos));

  pos=5; t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"alpha01 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";

  pos=6; t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"alpha11 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  m2_alpha11_est = theta_mod_NoRi1_unc_FINAL(pos);
  m2_alpha11_se = std::sqrt(CovThetaHat(pos,pos));


  pos=7; t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"alpha02 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp <<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";

  pos=8; t_temp = theta_mod_NoRi1_unc_FINAL(pos)/std::sqrt(CovThetaHat(pos,pos));
  Out<<"alpha12 (Est, SE, t, p-value): "<<theta_mod_NoRi1_unc_FINAL(pos)<<", "<<std::sqrt(CovThetaHat(pos,pos))\
  <<", "<<t_temp <<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n"; Out<<"\n\n";
  m2_alpha12_est = theta_mod_NoRi1_unc_FINAL(pos);
  m2_alpha12_se = std::sqrt(CovThetaHat(pos,pos));

  Out<<"A Wald test of alpha11=alpha12=0 determines whether MCAR holds\n";
  dbmat cov_alpha_NoR1i_mod(2,2);
  cov_alpha_NoR1i_mod = subm(CovThetaHat,range(6,2,8),range(6,2,8));
  Out<<"Cov(alpha11_hat, alpha12_hat): \n"<<cov_alpha_NoR1i_mod<<endl;
  dbcolvec Alp_vec_MCAR_NoR1i(2);
  Alp_vec_MCAR_NoR1i(0)= m2_alpha11_est;
  Alp_vec_MCAR_NoR1i(1)= m2_alpha12_est;
  double Wald_T_MCAR_NoR1i_mod = (trans(Alp_vec_MCAR_NoR1i)*inv(cov_alpha_NoR1i_mod)*Alp_vec_MCAR_NoR1i)(0,0);
  Out<<"Wald test statistic = "<<Wald_T_MCAR_NoR1i_mod<<endl;
  Out<<"P-value ="<< (1-scythe::pchisq(Wald_T_MCAR_NoR1i_mod,2))<<endl<<endl;

  dbmat SigMNAR_unc(dim,dim); dbcolvec muMNAR_unc(dim);
  MuEta_to_mu_sig(dim, theta_mod_NoRi1_unc_FINAL, muMNAR_unc, SigMNAR_unc);
  Out<<"Mu hat \n"<<trans(muMNAR_unc)<<"\n";
  Out<<"Sig hat \n"<<SigMNAR_unc<<"\n";
  Out<<"Maximized log-lik under H0UH1: "<<Final_max_log_likelihood_mod_NoRi1<<"\n\n";

  var1_NonIgno_NoRi1 = SigMNAR_unc(0,0);  var2_NonIgno_NoRi1 = SigMNAR_unc(1,1);
  cov12_NonIgno_NoRi1 = SigMNAR_unc(0,1);

   if( muMNAR_unc(0)<0 && muMNAR_unc(1)<0) { Out<<"Unconstrained MLE already consistent with H0, so  Tw and Tlr are equal to 0 and H0 can't be rejected. \n";
   theta_mod_NoRi1_con = theta_mod_NoRi1_unc_FINAL;
      }

   else {
   Out<<"------Constrained MLE\n";
   double MNARMod1_MinMlkUnderH0_bobyqa = -BIGNUM, MNARMod1_MinMlkUnderH0_bfgs = -BIGNUM;
   npt = 2*num_para_mod_NoRi1+1;

    dbcolvec theta_mod_NoRi1_conBFGS1(num_para_mod_NoRi1);
    theta_mod_NoRi1_con = theta_mod_NoRi1_unc_init;

   for(i=0;i<dim;++i)  if( theta_mod_NoRi1_con(i) > 0 ) theta_mod_NoRi1_con(i) = 0;
   pDat = &Dat; PMisInd = &Ind; BobyqaConvFail = 0;

   theta_mod_NoRi1_lb=-BIGNUM;  theta_mod_NoRi1_ub=BIGNUM;
   set_subm(theta_mod_NoRi1_ub,range(0,dim-1),range(0,0))=0;

// theta_mod_NoRi1_lb(6)= -5; theta_mod_NoRi1_lb(8)=-5; theta_mod_NoRi1_ub(6)=5; theta_mod_NoRi1_ub(8)=5;
   pDat = &Dat; PMisInd = &Ind;  pQHPts = &QHPts;
   MNARMod1_MinMlkUnderH0_bobyqa = find_min_bobyqaConvgFail(mod_NoRi1_MNARd2_MLogLikGH,
  theta_mod_NoRi1_con, npt,    // number of interpolation points
  theta_mod_NoRi1_lb,  // lower bound constraint
  theta_mod_NoRi1_ub,   // upper bound constraint
  8,    // initial trust region radius
  StopCriteria*0.1,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);

  if(BobyqaConvFail==1) { Out<<"Constrained MLE mod 1 fails~ "; return 0;}

  Out<<"theta_mod_NoRi1_con: "<<trans(theta_mod_NoRi1_con)<<"\n";
  Out<<"Gradient at theta_mod_NoRi1_con: "<<trans(mod_NoRi1_MNARd2_MLogLikGH_Grad(theta_mod_NoRi1_con))<<"\n";

  dbcolvec theta_mod_NoRi1_conBFGS(num_para_mod_NoRi1); theta_mod_NoRi1_conBFGS=theta_mod_NoRi1_con;
  BFGSFail = 0;  pDat = &Dat; PMisInd = &Ind;  pQHPts = &QHPts;
  MNARMod1_MinMlkUnderH0_bfgs = find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria*0.05), \
            mod_NoRi1_MNARd2_MLogLikGH, mod_NoRi1_MNARd2_MLogLikGH_Grad, theta_mod_NoRi1_conBFGS, theta_mod_NoRi1_lb,theta_mod_NoRi1_ub,BFGSFail,MaxIterations);
  if(BFGSFail==1) {  Out<<"BFGS Fails for constrained MLE\n "; goto cf1; }
  Out<<"theta_mod_NoRi1_conBFGS (verify constrained MLE): "<<trans(theta_mod_NoRi1_conBFGS)<<"\n";
  Out<<"maximized log-like under H0 using bfgs "<<(-1.0*mod_NoRi1_MNARd2_MLogLikGH(theta_mod_NoRi1_conBFGS) )<<endl;
  theta_mod_NoRi1_conBFGS1= theta_mod_NoRi1_conBFGS; theta_mod_NoRi1_conBFGS1(0) = -2.8;
  //Out<<"maximized log-like under H0  with first entry of  bfgs MLE adjusted to -2.8: "<<(-1.0*mod_NoRi1_MNARd2_MLogLikGH(theta_mod_NoRi1_conBFGS1))<<endl;

    theta_mod_NoRi1_conBFGS1= theta_mod_NoRi1_conBFGS; theta_mod_NoRi1_conBFGS1(0) = -6.8;
  //Out<<"maximized log-like under H0  with first entry of  bfgs MLE adjusted to -6.8: "<<(-1.0*mod_NoRi1_MNARd2_MLogLikGH(theta_mod_NoRi1_conBFGS1))<<endl;


    theta_mod_NoRi1_conBFGS1= theta_mod_NoRi1_conBFGS; theta_mod_NoRi1_conBFGS1(0) = -10.8;
 // Out<<"maximized log-like under H0 with first entry of  bfgs MLE adjusted to -10.8: "<<(-1.0*mod_NoRi1_MNARd2_MLogLikGH(theta_mod_NoRi1_conBFGS1))<<endl;


    if(MNARMod1_MinMlkUnderH0_bfgs < MNARMod1_MinMlkUnderH0_bobyqa)
    {  Out<<"True unconstrained MLE obtained by bfgs instead of bobyqa\n\n";
    theta_mod_NoRi1_con = theta_mod_NoRi1_conBFGS;
    MNARMod1_MinMlkUnderH0_bobyqa =  MNARMod1_MinMlkUnderH0_bfgs;
    }

 cf1:
    MuEta_to_mu_sig(dim, theta_mod_NoRi1_con, muMNAR_con, SigMNAR_con);
  Out<<"Mu hat \n"<<trans(muMNAR_con)<<"\n";
  Out<<"Sig hat \n"<<SigMNAR_con<<"\n";
  Out<<"Maximized log-lik under H0: "<<(-MNARMod1_MinMlkUnderH0_bobyqa)<<"\n\n";

  CovMuHatMNAR_h0 = subm(inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, theta_mod_NoRi1_con) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
  Out<<CovMuHatMNAR_h0<<"\n";

  CovMuHatMNAR_h0Uh1 = subm(inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, theta_mod_NoRi1_unc_FINAL) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
  Out<<CovMuHatMNAR_h0Uh1<<"\n";

  double TwMNAR_h0,TwMNAR_h0Uh1; dbmat InvCovMuHatMNAR(dim,dim);

// compute Tw
  InvCovMuHatMNAR = inv(CovMuHatMNAR_h0); get_Tw(TwMNAR_h0,dim,muMNAR_unc,InvCovMuHatMNAR,StopCriteria);
  Out<<"**** Using Cov(MuHat) version h0,  Tw="<<TwMNAR_h0<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
 pval = 1-(scythe::pchisq(TwMNAR_h0, double(dim-1))  \
        + scythe::pchisq(TwMNAR_h0, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";

 ChiBarWtsDim2ch2(CovMuHatMNAR_h0,CbWts); pval=ChiBarPvalue(dim,CbWts,TwMNAR_h0);
 Out<<"Substitution method p-value: "<<pval<<"\n";

  InvCovMuHatMNAR = inv(CovMuHatMNAR_h0Uh1); get_Tw(TwMNAR_h0Uh1,dim,muMNAR_unc,InvCovMuHatCC,StopCriteria);
  Out<<"\n**** Using Cov(MuHat) version h0UH1,  Tw="<<TwMNAR_h0Uh1<<"\n";
 pval= 1-(scythe::pchisq(TwMNAR_h0Uh1, double(dim-1))  \
        + scythe::pchisq(TwMNAR_h0Uh1, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";
 pval_Tw_NonIgno_bd_NoRi1m = pval;

 ChiBarWtsDim2ch2(CovMuHatMNAR_h0Uh1,CbWts); pval=ChiBarPvalue(dim,CbWts,TwMNAR_h0Uh1);
 Out<<"Substitution method p-value (with CovMuHatMNAR_h0Uh1): "<<pval<<"\n";
 pval_Tw_NonIgno_sub_NoRi1m = pval;

  double TlrMNAR=2*(MNARMod1_MinMlkUnderH0_bobyqa + Final_max_log_likelihood_mod_NoRi1);
  Out<<"\n\n Tlr_MNAR: "<<TlrMNAR<<"\n";

  pval = 1-(scythe::pchisq(TlrMNAR, double(dim-1))  \
        + scythe::pchisq(TlrMNAR, double(dim)) )/2;
  Out<<"Tlr_MNAR: upper bound of p-value ="<<pval<<"\n";
 pval_Tlr_NonIgno_bd_NoRi1m = pval;

  ChiBarWtsDim2ch2(CovMuHatMNAR_h0,CbWts);
  pval = ChiBarPvalue(dim,CbWts,TlrMNAR);
  Out<<"Tlr_MNAR: substitution method p-value using (Cov(MuHat) version h0) ="<<pval<<"\n\n";
 pval_Tlr_NonIgno_sub_NoRi1m = pval;

  ChiBarWtsDim2ch2(CovMuHatMNAR_h0Uh1,CbWts);
  pval = ChiBarPvalue(dim,CbWts,TlrMNAR);
  Out<<"Tlr_MNAR: substitution method p-value using (Cov(MuHat) version h0Uh1) ="<<pval<<"\n";
  }

   Out<<"\n\n\n#################\nAnalysis under MNAR model 1\n#################\n";
   Out<<"**********(xi1, xi2) are iid bivariate normal N(mu,Sig),  logit Pr(xi1 is missing)= alpha01 + alpha11*(xi1-mu1)/sig1,  and logit Pr(xi2 is missing)= alpha02 + alpha12*(xi2-mu2)/sig2 + alpha22 *ri1;\n  initial values of mu, Sig, alpha01, alpha11, alpha02, alpha12 are MLE from MNAR model 1. initial value of alpha22 is 0. \n\n";
  dbmat SigUncMLE(dim,dim);
  int num_para_mod2 = 10;
  dbcolvec theta_mod2_unc(num_para_mod2);  // mu, eta, alpha01,alpha11,alpha02,alpha12
  dbcolvec theta_mod2_unc_init(num_para_mod2);
  dbcolvec theta_mod2_lb(num_para_mod2),theta_mod2_ub(num_para_mod2);
  dbcolvec Grad_vec(num_para_mod2);  double minGrad,maxGrad;
  double MNARMod2_MinMlkUnderH0UH1;

  // Out<<"theta_mod1_unc: "<<trans(theta_mod1_unc)<<"\n";

  Out<<"------Unconstrained MLE\n";
  pDat = &Dat; PMisInd = &Ind;  pQHPts = &QHPts;
  //Out<<"initial value theta_mod1_unc: "<<trans(theta_mod1_unc)<<"\n";

  double MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL = BIGNUM;
  dbcolvec MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL(theta_mod2_unc_init);

  // unconstrained MLE: 6 groups of initial values
  for(i = 1; i<7; ++i)
  {
      switch(i)
      {
          case 1:
          Out<<"initial value option 1:  mu/Sig from unconstrained MLE using CC,  alpha01,02 from logit missing rates, and alpha11,alpha12 being 0\n";
  set_subm(theta_mod2_unc_init,range(0,4),range(0,0)) = thetaMuSigCC_unc;
  theta_mod2_unc_init(6) = 0; theta_mod2_unc_init(8) = 0;
  //   alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
  theta_mod2_unc_init(5) = std::log(  misrate1/(1-misrate1)  );
  theta_mod2_unc_init(7) = std::log(  misrate2/(1-misrate2)  );
   theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
         break;

        case 2:
        Out<<"initial value option 2:  mu/Sig from Constrained MLE using CC,  alpha01,02 from logit missing rates, and alpha11,alpha12 being 0\n";
  set_subm(theta_mod2_unc_init,range(0,4),range(0,0)) = thetaMuSigCC_con;
  theta_mod2_unc_init(6) = 0; theta_mod2_unc_init(8) = 0;
  //   alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
  theta_mod2_unc_init(5) = std::log(  misrate1/(1-misrate1)  );
  theta_mod2_unc_init(7) = std::log(  misrate2/(1-misrate2)  );
   theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
        break;

        case 3:
        Out<<"initial value option 3:  mu/Sig from Unconstrained MLE assuming MAR,  alpha01,02 from logit missing rates, and alpha11,alpha12 being 0\n";
  set_subm(theta_mod2_unc_init,range(0,4),range(0,0)) = thetaMuSigMAR_unc;
  theta_mod2_unc_init(6) = 0; theta_mod2_unc_init(8) = 0;
  //   alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
  theta_mod2_unc_init(5) = std::log(  misrate1/(1-misrate1)  );
  theta_mod2_unc_init(7) = std::log(  misrate2/(1-misrate2)  );
   theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
        break;

        case 4:
        Out<<"initial value option 4:  mu/Sig from Constrained MLE assuming MAR,  alpha01,02 from logit missing rates, and alpha11,alpha12 being 0\n";
  set_subm(theta_mod2_unc_init,range(0,4),range(0,0)) = thetaMuSigMAR_con;
  theta_mod2_unc_init(6) = 0; theta_mod2_unc_init(8) = 0;
  //   alpha01=theta(5); alpha11=theta(6); alpha02=theta(7); alpha12=theta(8);
  theta_mod2_unc_init(5) = std::log(  misrate1/(1-misrate1)  );
  theta_mod2_unc_init(7) = std::log(  misrate2/(1-misrate2)  );
   theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
       break;

        case 5:
        Out<<"initial value option 5:  except alpha22 (initial 0), other parameters are unconstrained MLE from the MNAR model without ri1 \n";
          set_subm(theta_mod2_unc_init,range(0,8),range(0,0)) = theta_mod_NoRi1_unc_FINAL;
          theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
          break;

        case 6:
        Out<<"initial value option 6:  except alpha22 (initial 0), other parameters are Constrained MLE from the MNAR model without ri1 \n";
   set_subm(theta_mod2_unc_init,range(0,8),range(0,0)) = theta_mod_NoRi1_con;
   theta_mod2_unc_init(9) = 0; Out<<trans(theta_mod2_unc_init)<<"\n";
        break;
      }

  theta_mod2_unc = theta_mod2_unc_init;
  theta_mod2_lb = -BIGNUM;  theta_mod2_ub = BIGNUM;
  BobyqaConvFail = 0; npt = 2*num_para_mod2+1;
  MNARMod2_MinMlkUnderH0UH1 = find_min_bobyqaConvgFail(mod2_MNARd2_MLogLikGH,
  theta_mod2_unc, npt,    // number of interpolation points
  theta_mod2_lb,  // lower bound constraint
  theta_mod2_ub,   // upper bound constraint
  8,    // initial trust region radius
  StopCriteria*0.1,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);
   unc_bobyqa_good = 0;
   if( BobyqaConvFail==0)
  {     Grad_vec = gradcdif(mod2_MNARd2_MLogLikGH,theta_mod2_unc);
        GetMinMax(Grad_vec, minGrad, maxGrad);
         if( maxGrad < GradientTolerance && minGrad > (-GradientTolerance) )  unc_bobyqa_good = 1;   }
   if(unc_bobyqa_good==1)
   {    Out<<"Bobyqa converges, the maximized log-likelihood is "<<(-1.0*MNARMod2_MinMlkUnderH0UH1)\
    <<"\n"<<"unconstrained MLE"<<trans(theta_mod2_unc)<<"Gradient at unconstrained MLE"<<trans(Grad_vec)<<"\n\n";
    if(MNARMod2_MinMlkUnderH0UH1 < MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL)
    {   MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL = MNARMod2_MinMlkUnderH0UH1;
       MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL = theta_mod2_unc;    }
   }
   else {   Out<<"Bobyqa method encounters non-convergence\n\n";  }

   theta_mod2_unc = theta_mod2_unc_init;
   NewtonConvFail = 0;
   MNARMod2_MinMlkUnderH0UH1
            =find_minConvFail(newton_search_strategy(mod2_MNARd2_MLogLikGH_Hess),
            objective_delta_stop_strategy(StopCriteria*10).be_verbose(), //
            mod2_MNARd2_MLogLikGH,
            mod2_MNARd2_MLogLikGH_Grad,
            theta_mod2_unc,
            -BIGNUM,NewtonConvFail,MaxIterations);
     unc_newton_good = 0;
      if( NewtonConvFail == 0)
    {  Grad_vec = gradcdif(mod2_MNARd2_MLogLikGH,theta_mod2_unc);
    GetMinMax(Grad_vec, minGrad, maxGrad);
   if( maxGrad < GradientTolerance && minGrad > (-GradientTolerance) )  unc_newton_good = 1;   }
    if(unc_newton_good==1)
   {     Out<<"Newton converges, the maximized log-likelihood is "<<(-1.0*MNARMod2_MinMlkUnderH0UH1)\
    <<"\n"<<"unconstrained MLE"<<trans(theta_mod2_unc)<<"Gradient at unconstrained MLE"<<trans(Grad_vec)<<"\n\n";
    if(MNARMod2_MinMlkUnderH0UH1 < MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL)
    {   MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL = MNARMod2_MinMlkUnderH0UH1;
       MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL = theta_mod2_unc;    }
   }
   else {   Out<<"Newton method encounters non-convergence\n\n";  }

   theta_mod2_unc = theta_mod2_unc_init;   bfgsFail = 0;
    MNARMod2_MinMlkUnderH0UH1 =  find_minConvFail(bfgs_search_strategy(),
           objective_delta_stop_strategy(StopCriteria),
           mod2_MNARd2_MLogLikGH, mod2_MNARd2_MLogLikGH_Grad, \
           theta_mod2_unc, -BIGNUM,bfgsFail,MaxIterations);
    unc_bfgs_good = 0;
    if( bfgsFail == 0)
   {  Grad_vec = gradcdif(mod2_MNARd2_MLogLikGH,theta_mod2_unc);
    GetMinMax(Grad_vec, minGrad, maxGrad);
   if( maxGrad < GradientTolerance && minGrad > (-GradientTolerance) )  unc_bfgs_good = 1;   }
   if(unc_bfgs_good==1)
   {    Out<<"BFGS converges, the maximized log-likelihood is "<<(-1.0*MNARMod2_MinMlkUnderH0UH1)\
    <<"\n"<<"unconstrained MLE"<<trans(theta_mod2_unc)<<"Gradient at unconstrained MLE"<<trans(Grad_vec)<<"\n\n";
        if(MNARMod2_MinMlkUnderH0UH1 < MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL)
         {   MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL = MNARMod2_MinMlkUnderH0UH1;
          MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL = theta_mod2_unc;    }
   }
   else {   Out<<"BFGS method encounters non-convergence\n\n";  }

   //Out<<"use Trust region \n\n";
   theta_mod2_unc = theta_mod2_unc_init;   TrustRegionConvFail = 0;
    MNARMod2_MinMlkUnderH0UH1 = \
     find_min_trust_regionConvFail(objective_delta_stop_strategy(StopCriteria),\
    mod2_MNARd2_MLogLikGH_class(), theta_mod2_unc, TrustRegionConvFail,MaxIterations,10 );
     unc_tregion_good = 0;
     if( TrustRegionConvFail == 0)
      {  Grad_vec = gradcdif(mod2_MNARd2_MLogLikGH,theta_mod2_unc);
    GetMinMax(Grad_vec, minGrad, maxGrad);
              if( maxGrad < GradientTolerance && minGrad > (-GradientTolerance) )  unc_tregion_good = 1;   }
    if(unc_tregion_good==1)
   {     Out<<"Trust region method converges, the maximized log-likelihood is "<<(-1.0*MNARMod2_MinMlkUnderH0UH1)\
    <<"\n"<<"unconstrained MLE"<<trans(theta_mod2_unc)<<"Gradient at unconstrained MLE"<<trans(Grad_vec)<<"\n\n";
        if(MNARMod2_MinMlkUnderH0UH1 < MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL)
         {   MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL = MNARMod2_MinMlkUnderH0UH1;
          MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL = theta_mod2_unc;    }
   }
   else {   Out<<"Trust region method encounters non-convergence\n\n";  }

  }

   theta_mod2_unc = MNARMod_withR1i_Min_Unc_MLE_UnderH0UH1_FINAL;
   MNARMod2_MinMlkUnderH0UH1 = MNARMod_withR1i_Min_MinusLogLike_UnderH0UH1_FINAL;

  dbmat CovThetaHat_mod2(num_para_mod2,num_para_mod2);
   CovThetaHat_mod2 = inv( HESScdif(mod2_MNARd2_MLogLikGH, theta_mod2_unc) );

  pos=0;  t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"mu1 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu1_NonIgnom1_Est = theta_mod2_unc(pos);
  mu1_NonIgnom1_SE = std::sqrt(CovThetaHat_mod2(pos,pos));

  pos=1;  t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"mu2 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  mu2_NonIgnom1_Est = theta_mod2_unc(pos);
  mu2_NonIgnom1_SE = std::sqrt(CovThetaHat_mod2(pos,pos));

  pos=5; t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"alpha01 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";

  pos=6; t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"alpha11 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  m1_alpha11_est = theta_mod2_unc(pos);
  m1_alpha11_se = std::sqrt(CovThetaHat_mod2(pos,pos));

  pos=7; t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"alpha02 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";

  pos=8; t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"alpha12 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  m1_alpha12_est = theta_mod2_unc(pos);
  m1_alpha12_se = std::sqrt(CovThetaHat_mod2(pos,pos));

  pos=9; t_temp = theta_mod2_unc(pos)/std::sqrt(CovThetaHat_mod2(pos,pos));
  Out<<"alpha22 (Est, SE, t, p-value): "<<theta_mod2_unc(pos)<<", "<<std::sqrt(CovThetaHat_mod2(pos,pos))\
  <<", "<<t_temp<<", "<< (1-scythe::pnorm(std::abs(t_temp),0,1) )<<"\n";
  m1_alpha22_est = theta_mod2_unc(pos);
  m1_alpha22_se = std::sqrt(CovThetaHat_mod2(pos,pos));

  Out<<"\n\n";
  Out<<"A Wald test of alpha11=alpha12=alpha22=0 determines whether MCAR holds\n";
  dbmat cov_alpha_WithR1i_mod(3,3);
   int row_i, col_j;
  for(i = 6; i<11; i = i + 2)
  {
      row_i = i/2 - 3;
      if(i == 10)  i = 9;
    for(j = 6; j<11; j = j + 2)
    {    col_j = j/2 - 3;    if(j == 10)  j = 9;
      cov_alpha_WithR1i_mod(row_i,col_j) = CovThetaHat_mod2(i,j);
    }
  }

  Out<<"Cov(alpha11_hat, alpha12_hat, alpha22_hat): \n"<<cov_alpha_WithR1i_mod<<endl;
  Out<<"Full CovThetaHat_mod2\n"<<CovThetaHat_mod2<<endl;

  dbcolvec Alp_vec_MCAR_WithR1i(3);
  Alp_vec_MCAR_WithR1i(0)= m1_alpha11_est;
  Alp_vec_MCAR_WithR1i(1)= m1_alpha12_est;  Alp_vec_MCAR_WithR1i(1)= m1_alpha22_est;
  double Wald_T_MCAR_WithR1i_mod = (trans(Alp_vec_MCAR_WithR1i)*inv(cov_alpha_WithR1i_mod)*Alp_vec_MCAR_WithR1i)(0,0);
  Out<<"Wald test statistic = "<<Wald_T_MCAR_WithR1i_mod<<endl;
  Out<<"P-value ="<< (1-scythe::pchisq(Wald_T_MCAR_WithR1i_mod,3))<<endl<<endl;

  MuEta_to_mu_sig(dim, theta_mod2_unc, muMNAR_unc, SigMNAR_unc);
  SigUncMLE = SigMNAR_unc;
  Out<<"Mu hat \n"<<trans(muMNAR_unc)<<"\n";
  Out<<"Sig hat \n"<<SigMNAR_unc<<"\n";
  Out<<"Maximized log-lik under H0UH1: "<<(-MNARMod2_MinMlkUnderH0UH1)<<"\n\n";

  var1_NonIgnom1 = SigMNAR_unc(0,0);  var2_NonIgnom1 = SigMNAR_unc(1,1);
  cov12_NonIgnom1 = SigMNAR_unc(0,1);

   if( muMNAR_unc(0)<0 && muMNAR_unc(1)<0) { Out<<"Unconstrained MLE already consistent with H0, so  Tw and Tlr are equal to 0 and H0 can't be rejected. \n";
    Out<<"Tw=0, and Tlr = 0\n The upper bounds of p-values for both Tw and Tlr  are 1.\n";
   dbmat CovThetaHat_mod2(num_para_mod2,num_para_mod2);
    CovThetaHat_mod2 = inv( HESScdif(mod2_MNARd2_MLogLikGH, theta_mod2_unc) );
  CovMuHatMNAR_h0 = subm(CovThetaHat_mod2,range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0: (inverse) observed information matrix evaluated at the unconstrained (also H0-constrained) MLE \n";
  Out<<CovMuHatMNAR_h0<<"\n";
  TwMNAR_h0Uh1 = 1;
  ChiBarWtsDim2ch2(CovMuHatMNAR_h0,CbWts); pval = ChiBarPvalue(dim,CbWts,TwMNAR_h0Uh1);
  Out<<"Substitution method p-value for both Tw and Tlr: "<<pval<<"\n";
  pval_Tw_NonIgno_sub_WithRi1m = pval;
  pval_Tlr_NonIgno_sub_WithRi1m = pval;
  pval_Tlr_NonIgno_bd_WithRi1m = 1; pval_Tw_NonIgno_bd_WithRi1m = 1;
  }

  else {
  Out<<"------Constrained MLE\n";
  double MNARMod2_MinMlkUnderH0_bobyqa, MNARMod2_MinMlkUnderH0_bfgs;

  pDat = &Dat; PMisInd = &Ind;  pQHPts = &QHPts;
  npt = 2*num_para_mod2+1;  dbcolvec theta_mod2_con(num_para_mod2);
  theta_mod2_con = theta_mod2_unc;
  for(i=0;i<dim;++i)  if( theta_mod2_con(i) > 0 ) theta_mod2_con(i) = 0;
  pDat = &Dat;   PMisInd = &Ind;    BobyqaConvFail = 0;

  theta_mod2_lb = -BIGNUM;  theta_mod2_ub = BIGNUM;
  set_subm(theta_mod2_ub,range(0,dim-1),range(0,0)) = 0;

// theta_mod1_lb(6)= -5; theta_mod1_lb(8)=-5; theta_mod1_ub(6)=5; theta_mod1_ub(8)=5;

  MNARMod2_MinMlkUnderH0_bobyqa = find_min_bobyqaConvgFail(mod2_MNARd2_MLogLikGH,
  theta_mod2_con, npt,    // number of interpolation points
  theta_mod2_lb,  // lower bound constraint
  theta_mod2_ub,   // upper bound constraint
  8,    // initial trust region radius
  StopCriteria,  // stopping trust region radius
  MaxFunEvalLong, // max number of objective function evaluations
  BobyqaConvFail);

  if(BobyqaConvFail==1 ) { Out<<"Constrained MLE mod 1 fails~ "; return 0;}

  Out<<"theta_mod2_con: "<<trans(theta_mod2_con)<<"\n";
  Out<<"Gradient at theta_mod2_con: "<<trans(mod2_MNARd2_MLogLikGH_Grad(theta_mod2_con))<<"\n";

  dbcolvec theta_mod2_conBFGS(num_para_mod2);
   theta_mod2_conBFGS = theta_mod2_con; BFGSFail = 0;
 MNARMod2_MinMlkUnderH0_bfgs =  find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria), \
            mod2_MNARd2_MLogLikGH, mod2_MNARd2_MLogLikGH_Grad, theta_mod2_conBFGS, theta_mod2_lb,theta_mod2_ub,BFGSFail,MaxIterations);
  if(BFGSFail==1) {  Out<<"BFGS Fails for constrained MLE\n ";  goto cf2;}
  Out<<"theta_mod2_conBFGS (verify constrained MLE): "<<trans(theta_mod2_conBFGS)<<"\n";

  if(MNARMod2_MinMlkUnderH0_bfgs > MNARMod2_MinMlkUnderH0_bobyqa)
  {   Out<<"True unconstrained MLE obtained by bfgs\n\n";
   MNARMod2_MinMlkUnderH0_bobyqa = MNARMod2_MinMlkUnderH0_bfgs;
   theta_mod2_con = theta_mod2_conBFGS ; }

 cf2:
  MuEta_to_mu_sig(dim, theta_mod2_con, muMNAR_con, SigMNAR_con);
  Out<<"Mu hat \n"<<trans(muMNAR_con)<<"\n";
  Out<<"Sig hat \n"<<SigMNAR_con<<"\n";
  Out<<"Maximized log-lik under H0: "<<(-MNARMod2_MinMlkUnderH0_bobyqa)<<"\n\n";

  CovMuHatMNAR_h0 = subm(inv( HESScdif(mod2_MNARd2_MLogLikGH, theta_mod2_con) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
  Out<<CovMuHatMNAR_h0<<"\n";

  CovMuHatMNAR_h0Uh1=subm(inv( HESScdif(mod2_MNARd2_MLogLikGH, theta_mod2_unc) ),range(0,dim-1),range(0,dim-1) ) ;
  Out<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
  Out<<CovMuHatMNAR_h0Uh1<<"\n";

// compute Tw
  bool CovMuHatMNAR_h0_posi_definite = is_positive_definite_2by2(CovMuHatMNAR_h0);
  if(CovMuHatMNAR_h0_posi_definite)
  {
  InvCovMuHatMNAR =inv(CovMuHatMNAR_h0); get_Tw(TwMNAR_h0,dim,muMNAR_unc,InvCovMuHatMNAR,StopCriteria);
  Out<<"**** Using Cov(MuHat) version h0,  Tw="<<TwMNAR_h0<<"\n";

// upper bound is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
        // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
 pval = 1-(scythe::pchisq(TwMNAR_h0, double(dim-1))  \
        + scythe::pchisq(TwMNAR_h0, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";

 ChiBarWtsDim2ch2(CovMuHatMNAR_h0,CbWts); pval=ChiBarPvalue(dim,CbWts,TwMNAR_h0);
 Out<<"Substitution method p-value: "<<pval<<"\n";
}

  InvCovMuHatMNAR=inv(CovMuHatMNAR_h0Uh1); get_Tw(TwMNAR_h0Uh1,dim,muMNAR_unc,InvCovMuHatCC,StopCriteria);
  Out<<"\n**** Using Cov(MuHat) version h0UH1,  Tw="<<TwMNAR_h0Uh1<<"\n";
 pval = 1-(scythe::pchisq(TwMNAR_h0Uh1, double(dim-1))  \
        + scythe::pchisq(TwMNAR_h0Uh1, double(dim)) )/2;
 Out<<"Upper bound of p-value: "<<pval<<"\n";
 pval_Tw_NonIgno_bd_WithRi1m = pval;

 ChiBarWtsDim2ch2(CovMuHatMNAR_h0Uh1,CbWts); pval = ChiBarPvalue(dim,CbWts,TwMNAR_h0Uh1);
 Out<<"Substitution method p-value using CovMuHatMNAR_h0Uh1: "<<pval<<"\n";
 pval_Tw_NonIgno_sub_WithRi1m = pval;

   TlrMNAR = 2*(MNARMod2_MinMlkUnderH0_bobyqa - MNARMod2_MinMlkUnderH0UH1);
  Out<<"\n\n Tlr_MNAR: "<<TlrMNAR<<"\n";

  pval = 1-(scythe::pchisq(TlrMNAR, double(dim-1))  \
        + scythe::pchisq(TlrMNAR, double(dim)) )/2;
  Out<<"Tlr_MNAR: upper bound of p-value ="<<pval<<"\n";
  pval_Tlr_NonIgno_bd_WithRi1m = pval;

  if(CovMuHatMNAR_h0_posi_definite)
  { ChiBarWtsDim2ch2(CovMuHatMNAR_h0,CbWts);
    pval = ChiBarPvalue(dim,CbWts,TlrMNAR);
    Out<<"Tlr_MNAR: substitution method p-value using (Cov(MuHat) version h0) ="<<pval<<"\n\n";
    pval_Tlr_NonIgno_sub_WithRi1m = pval;
 }

  ChiBarWtsDim2ch2(CovMuHatMNAR_h0Uh1,CbWts);
  pval = ChiBarPvalue(dim,CbWts,TlrMNAR);
  Out<<"Tlr_MNAR: substitution method p-value using (Cov(MuHat) version h0Uh1) ="<<pval<<"\n";
 }


  Out<<"\n\nManuscript table of unconstrained MLEs (mu) :\n";
  Out<<setprecision(1)<<fixed<<"& &"<<mu1_CC_Est<<" &"<<mu1_CC_SE<<"  &&"\
  <<mu2_CC_Est<<"&"<<mu2_CC_SE<<"&&"<<mu1_Igno_Est<<"&"\
  <<mu1_Igno_SE<<"&&"<<mu2_Igno_Est<<"&"<<mu2_Igno_SE\
  <<"&&"<<mu1_NonIgnom1_Est<<"&"<<mu1_NonIgnom1_SE\
  <<"&&"<<mu2_NonIgnom1_Est<<"&"<<mu2_NonIgnom1_SE\
  <<"&&"<<mu1_NonIgno_NoRi1_Est<<"&"<<mu1_NonIgno_NoRi1_SE\
  <<"&&"<<mu2_NonIgno_NoRi1_Est<<"&"<<mu2_NonIgno_NoRi1_SE\
  <<"\n\n";
 // & & EST & SE & & EST &SE & & EST & SE & & EST & SE && EST &SE & &EST &SE
//  mu1_Igno_Est  mu1_Igno_Est  mu2_NonIgnom1_Est


  Out<<"\n\nManuscript table of unconstrained MLEs (Sig) :\n";
  Out<<setprecision(0)<<fixed<<"& &"<<var1_CC<<"&"<<var2_CC<<"&"<<"\\multicolumn{3}{c}{"<<cov12_CC<<"}"<<"&&"\
  <<var1_Igno<<"&"<<var2_Igno<<"&"<<"\\multicolumn{3}{c}{"<<cov12_Igno<<"}"<<"&&"\
  <<var1_NonIgnom1<<"&"<<var2_NonIgnom1<<"&"<<"\\multicolumn{3}{c}{"<<cov12_NonIgnom1<<"}"<<"&&"\
  <<var1_NonIgno_NoRi1<<"&"<<var2_NonIgno_NoRi1<<"&"<<"\\multicolumn{3}{c}{"<<cov12_NonIgno_NoRi1<<"}"<<"\n\n";


  // & &0.1&0.9& &-0.7&0.6& &2.5&0.7&&\multicolumn{2}{c}{Yes}&&0.1&1.0& &-0.7&0.9&&\multicolumn{2}{c}{Yes}
  Out<<"\n\nManuscript table of unconstrained MLEs (Missing model parameters) :\n";
  Out<<setprecision(1)<<fixed<<"& &"<<m1_alpha11_est<<"&"<<m1_alpha11_se\
  <<"& &"<<m1_alpha12_est<<"&"<<m1_alpha12_se<<"& &"\
  <<m1_alpha22_est<<"&"<<m1_alpha22_se<<"&&\\multicolumn{2}{c}{Yes}"<<"&&"\
  <<m2_alpha11_est<<"&"<<m2_alpha11_se\
  <<"& &"<<m2_alpha12_est<<"&"<<m2_alpha12_se<<"&&"\
  <<"&&\\multicolumn{2}{c}{Yes}"<<"\n\n";


  Out<<"\nTable for substitution p-values of various tests\n\n & "<<setprecision(2)<<fixed\
  <<pval_Tlr_CC_sub<<"& &"<<pval_Tw_Igno_sub<<"& "<<pval_Tlr_Igno_sub\
  <<"& &"<<pval_Tw_NonIgno_sub_WithRi1m<<"& "<<pval_Tlr_NonIgno_sub_WithRi1m<<"& &"\
  <<pval_Tw_NonIgno_sub_NoRi1m<<"& "<<pval_Tlr_NonIgno_sub_NoRi1m<<"\n\n";

  Out<<"\nTable for bound p-values of various tests\n\n & "<<setprecision(2)<<fixed\
  <<pval_Tlr_CC_bd<<"& &"<<pval_Tw_Igno_bd<<"& "<<pval_Tlr_Igno_bd\
  <<"& &"<<pval_Tw_NonIgno_bd_WithRi1m<<"& "<<pval_Tlr_NonIgno_bd_WithRi1m<<"& &"\
  <<pval_Tw_NonIgno_bd_NoRi1m<<"& "<<pval_Tlr_NonIgno_bd_NoRi1m<<"\n\n";

 Out<<"\n##############\nThe version with 5 decimal digits\n##############\n";

 Out<<"\n\nManuscript table of unconstrained MLEs (mu) :\n";
  Out<<setprecision(5)<<fixed<<"& &"<<mu1_CC_Est<<" &"<<mu1_CC_SE<<"  &&"\
  <<mu2_CC_Est<<"&"<<mu2_CC_SE<<"&&"<<mu1_Igno_Est<<"&"\
  <<mu1_Igno_SE<<"&&"<<mu2_Igno_Est<<"&"<<mu2_Igno_SE\
  <<"&&"<<mu1_NonIgnom1_Est<<"&"<<mu1_NonIgnom1_SE\
  <<"&&"<<mu2_NonIgnom1_Est<<"&"<<mu2_NonIgnom1_SE\
  <<"&&"<<mu1_NonIgno_NoRi1_Est<<"&"<<mu1_NonIgno_NoRi1_SE\
  <<"&&"<<mu2_NonIgno_NoRi1_Est<<"&"<<mu2_NonIgno_NoRi1_SE\
  <<"\n\n";
 // & & EST & SE & & EST &SE & & EST & SE & & EST & SE && EST &SE & &EST &SE
//  mu1_Igno_Est  mu1_Igno_Est  mu2_NonIgnom1_Est


  Out<<"\n\nManuscript table of unconstrained MLEs (Sig) :\n";
  Out<<setprecision(5)<<fixed<<"& &"<<var1_CC<<"&"<<var2_CC<<"&"<<"\\multicolumn{3}{c}{"<<cov12_CC<<"}"<<"&&"\
  <<var1_Igno<<"&"<<var2_Igno<<"&"<<"\\multicolumn{3}{c}{"<<cov12_Igno<<"}"<<"&&"\
  <<var1_NonIgnom1<<"&"<<var2_NonIgnom1<<"&"<<"\\multicolumn{3}{c}{"<<cov12_NonIgnom1<<"}"<<"&&"\
  <<var1_NonIgno_NoRi1<<"&"<<var2_NonIgno_NoRi1<<"&"<<"\\multicolumn{3}{c}{"<<cov12_NonIgno_NoRi1<<"}"<<"\n\n";

  Out<<"\n\nManuscript table of unconstrained MLEs (Missing model parameters) :\n";
  Out<<setprecision(5)<<fixed<<"& &"<<m1_alpha11_est<<"&"<<m1_alpha11_se\
  <<"& &"<<m1_alpha12_est<<"&"<<m1_alpha12_se<<"& &"\
  <<m1_alpha22_est<<"&"<<m1_alpha22_se<<"&&\\multicolumn{2}{c}{Yes}"<<"&&"\
  <<m2_alpha11_est<<"&"<<m2_alpha11_se\
  <<"& &"<<m2_alpha12_est<<"&"<<m2_alpha12_se<<"&&"\
  <<"\\multicolumn{2}{c}{Yes}"<<"\n\n";


  Out<<"\nTable for substitution p-values of various tests\n\n & "<<setprecision(5)<<fixed\
  <<pval_Tlr_CC_sub<<"& &"<<pval_Tw_Igno_sub<<"& "<<pval_Tlr_Igno_sub\
  <<"& &"<<pval_Tw_NonIgno_sub_WithRi1m<<"& "<<pval_Tlr_NonIgno_sub_WithRi1m<<"& &"\
  <<pval_Tw_NonIgno_sub_NoRi1m<<"& "<<pval_Tlr_NonIgno_sub_NoRi1m<<"\n\n";

  Out<<"\nTable for bound p-values of various tests\n\n & "<<setprecision(5)<<fixed\
  <<pval_Tlr_CC_bd<<"& &"<<pval_Tw_Igno_bd<<"& "<<pval_Tlr_Igno_bd\
  <<"& &"<<pval_Tw_NonIgno_bd_WithRi1m<<"& "<<pval_Tlr_NonIgno_bd_WithRi1m<<"& &"\
  <<pval_Tw_NonIgno_bd_NoRi1m<<"& "<<pval_Tlr_NonIgno_bd_NoRi1m<<"\n\n";

  return 0; }
