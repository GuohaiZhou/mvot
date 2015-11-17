#include"header.h"

int main()
{
    using namespace std;
    int mQHpts = 100; dbmat QHPts(mQHpts,2);
    ifstream QHfile; QHfile.open("GH100.txt");
    GetQHPts(QHfile,mQHpts,QHPts); QHfile.close();

   // double initial_trust_region = 10;

    ifstream in1;
    in1.open("simulation_configurations.txt");
    string str,seedO,MisType;
    stringstream ss,ss2;
    ofstream osum,otLog,osumDts;
    osum.open("Summary.txt");
    otLog.open("SimuLog.txt");
    osumDts.open("SummaryDetails.txt");
    string ManualConfirm;
    int MaxFunEval,MaxIterations;
    double StopCriteria;
    long int MaxFunEvalLong;
    ofstream NumerErrorInstance;
    NumerErrorInstance.open("ExampleOfNumericalFailure.txt");
// int NumerErrorDatExported=0;

    int NewtonConvFail, bfgsConvFail,BobyqaConvFail, TrustRegionConvFail;
    // read simulation setting
    int i,j,NumSimu,Sz;
    int dim = 2;
    int dimEta = dim*(dim+1)/2;
    int dimAlpha = 2*dim;  int dimAlphaNew = dimAlpha + 1;

    NumSimu = atoi(str.c_str());
    Sz = atoi(str.c_str());

    dbcolvec mu(dim),Var1Var2Cov12(dimEta),alphaVec(dimAlpha), alphaVecNew(dimAlphaNew);


    for(i = 1; i < 32; ++i)
    {
        switch (i)
        {
        case 4: // integers
            getline(in1,str);
            NumSimu = atoi(str.c_str());
            break;

        case 7: //  integers
            getline(in1,str);
            Sz = atoi(str.c_str());
            break;

        case 10:  // vectors
            for(j=0; j < (dim-1); ++j)
            {
                getline(in1,str,',');
                mu(j) = strtod(str.c_str(),NULL);
            }
            getline(in1,str);
            mu(dim-1) = strtod(str.c_str(),NULL);
            break;

        case 13:  // vectors
            for(j = 0; j < (dimEta-1); ++j)
            {
                getline(in1,str,',');
                Var1Var2Cov12(j)=strtod(str.c_str(),NULL);
            }
            getline(in1,str);
            Var1Var2Cov12(dimEta-1)=strtod(str.c_str(),NULL);
            break;

        case 16:  // vectors
            for(j = 0; j < (dimAlphaNew-1); ++j)
            {
                getline(in1,str,',');
                alphaVecNew(j)=strtod(str.c_str(),NULL);
            }
            getline(in1,str);
            alphaVecNew(dimAlphaNew-1) = strtod(str.c_str(),NULL);
            break;

        case 19: // string
            getline(in1,seedO);
            break;

        case 22: // double
            getline(in1,str);
            StopCriteria = strtod(str.c_str(),NULL);
            break;

        case 25: // for integers
            getline(in1,str);
            MaxFunEval = atoi(str.c_str());
            break;

        case 28:
            getline(in1,str);
            MaxIterations = atoi(str.c_str());
            break;

        case 31: // for strings
            getline(in1, ManualConfirm);
            break;

        default:
            getline(in1,str);
        }
    }
    alphaVec = subm(alphaVecNew,range(0,dimAlpha-1),range(0,0));
    double alphaMAR = alphaVecNew(dimAlphaNew-1);
    dbmat Sig(dim,dim);
    Sig(0,0) = Var1Var2Cov12(0);
    Sig(1,1) = Var1Var2Cov12(1);
    Sig(0,1) = Var1Var2Cov12(2);
    Sig(1,0) = Sig(0,1);
    double sd1 = std::sqrt(Var1Var2Cov12(0)), sd2 = std::sqrt(Var1Var2Cov12(1));

    double StopCriteriaTrustRegion = StopCriteria*10;
    // trust region method can have a looser stop criterion
    MaxFunEvalLong = long(MaxFunEval);

 // double BigPositiveEstimate=200,BigNegativeEstimate=-200;

    cout<<"Simulation Settings: \n"<<endl;
    cout<<"Number of simulation replications: "<<NumSimu<<endl<<endl;
    cout<<"Sample size: "<<Sz<<endl<<endl;
    cout<<"true value of mu: "<<trans(mu)<<endl<<endl;
    cout<<"true value of Sig: \n"<<Sig<<endl<<endl;
    cout<<"true value of alpha: \n"<<trans(alphaVecNew)<<endl<<endl;
    cout<<"seed: "<<seedO<<endl<<endl;
    cout<<"Stopping criteria: "<<StopCriteria<<"\n";
    cout<<"Maximum number of likelihood evaluations for bobyqa: "<<MaxFunEval<<"\n";
    cout<<"Maximum number of iterations evaluations for newton, trust-region and bfgs: "<<MaxIterations<<"\n";

    string confirmed;
    if(ManualConfirm.compare("N") == 0 || ManualConfirm.compare("n") == 0)
        confirmed = "y";
    else if (ManualConfirm.compare("Y") == 0 || ManualConfirm.compare("y") == 0)
    {
        cout<<"Please confirm the above settings:\n";
        cout<<"type Y or y to start the simulation:\n "<<endl;
        cout<<"type N or n to exit if the settings need changes:\n "<<endl;
        cin>>confirmed;
    }
    else   confirmed ="else";

    if(confirmed == "Y" || confirmed == "y")
    {   // set up random num generator
        dlib::rand rndSample,rndMisInd;
        rndSample.set_seed(seedO);

        double SeedForMisInd = rndSample.get_random_gaussian();
        ss2<<SeedForMisInd;
        rndMisInd.set_seed( ss2.str() );

        // variables used to record running time
        clock_t tStartTw, tStartTlr;
    //double Seconds_Tw_CC_average = 0,Seconds_Tlr_CC_average = 0;
    //double Seconds_Tw_Igno_average = 0,Seconds_Tlr_Igno_average = 0;
    // double Seconds_Tw_NonIgno_average = 0,Seconds_Tlr_NonIgno_average = 0;

        // temp var to record running time
        double Seconds_unconstrained_MLE,Seconds_constrained_MLE,Seconds_Step2_Tw;
        double Seconds_Tw_Igno, Seconds_Tlr_Igno,Seconds_Tw_NonIgno,Seconds_Tlr_NonIgno;

        // allocate space ( no re-allocation in each simu loop!)
        int numPara = dim + dimEta + dimAlpha;
        int ParaMuSig = dim + dimEta;
        dbmat NonIgno_Avg_CovThetaHat_vh0Uh1(numPara,numPara);
        dbmat NonIgno_Avg_CovThetaHat_vh0(numPara,numPara);
        dbmat Igno_Avg_CovThetaHat_vh0Uh1(ParaMuSig,ParaMuSig);
        dbmat Igno_Avg_CovThetaHat_vh0(ParaMuSig,ParaMuSig);
        dbmat CC_Avg_CovThetaHat_vh0Uh1(ParaMuSig,ParaMuSig);
        dbmat CC_Avg_CovThetaHat_vh0(ParaMuSig,ParaMuSig);
        // estimated cov(ThetaHat) over NumSimu instances, to be compared with actual sampling cov(ThetaHat)

      // temporary matrices/vars
        dbmat NonIgno_CovThetaHat_this(numPara,numPara),Igno_CovThetaHat_this(ParaMuSig,ParaMuSig);
        dbmat CC_CovThetaHat_this(ParaMuSig,ParaMuSig);
        dbmat CovMuHat_this_vh0(dim,dim),InvCovMuHat_this_vh0(dim,dim);
        dbmat CovMuHat_this_vh0Uh1(dim,dim),InvCovMuHat_this_vh0Uh1(dim,dim);
        dbmat R_ThetaZeroInv_Rt_vh0(dim,dim),R_ThetaZeroInv_Rt_vh0Uh1(dim,dim); // used for substitution p-val
        dbmat SamCovInit(dim,dim);  SamCovInit = 0;
        dbcolvec SamMeanInit(dim);  SamMeanInit = 0;
        double CC_MinNegLogLike_UnderH0, CC_MinNegLogLike_UnderH0UH1;
        double Igno_MinNegLogLike_UnderH0, Igno_MinNegLogLike_UnderH0UH1;
        double NonIgno_MinNegLogLike_UnderH0, NonIgno_MinNegLogLike_UnderH0UH1;
        dbmat Dat(Sz,dim),MisInd(Sz,dim),DatCCtemp(Sz,dim);
        dbmat Dat_BeforeRescal(Sz,dim);
        int NumComCases;
        dbcolvec MissingRates(dim);
        dbcolvec lb_Igno(ParaMuSig,ParaMuSig), ub_Igno(ParaMuSig,ParaMuSig);
        dbcolvec lb_NonIgno(numPara,numPara), ub_NonIgno(numPara,numPara);
        dbcolvec CbWts_vh0Uh1(dim+1),CbWts_vh0(dim+1);
        long int npt_Igno = 2*ParaMuSig + 1, npt_NonIgno = 2*numPara + 1;

        dbcolvec CC_theta_unc(ParaMuSig), CC_theta_con(ParaMuSig), CC_Eta_unc(dimEta);
        dbmat CC_L_unc(dim,dim);
        dbcolvec CC_mu_con(dim), CC_Eta_con(dimEta); dbmat CC_L_con(dim,dim);

        dbcolvec Igno_theta_unc(ParaMuSig), Igno_theta_con(ParaMuSig), Igno_Eta_unc(dimEta);
        dbmat Igno_L_unc(dim,dim);  dbcolvec Igno_mu_unc(dim);
        dbcolvec Igno_mu_con(dim), Igno_Eta_con(dimEta); dbmat Igno_L_con(dim,dim);

        dbcolvec NonIgno_theta_unc(numPara), NonIgno_theta_con(numPara), NonIgno_Eta_unc(dimEta);
        dbmat NonIgno_L_unc(dim,dim);  dbcolvec NonIgno_mu_unc(dim);
        dbcolvec NonIgno_mu_con(dim), NonIgno_Eta_con(dimEta); dbmat NonIgno_L_con(dim,dim);

        dbcolvec Igno_Gradient(ParaMuSig), NonIgno_Gradient(numPara); double MinGrad,MaxGrad;
        dbcolvec Igno_theta_initial_con(ParaMuSig), Igno_theta_initial_unc(numPara);
        // reserved for initial value
        dbcolvec NonIgno_theta_initial_con(numPara), NonIgno_theta_initial_unc(numPara);
        dbcolvec NonIgno_theta_initial_unc_bfgs(numPara);
        // reserved for initial value

     //   dbcolvec NonIgno_theta_true(numPara);
        // true parameter vector, if non-convergence occurs, use this value as initial value

        dbmat NonIgno_Sig_unc(dim,dim), NonIgno_Sig_con(dim,dim);

        double pval_Tlr_CC_bd,pval_Tw_CC_bd, pval_Tw_Igno_bd;
        double pval_Tlr_Igno_bd, pval_Tw_NonIgno_bd, pval_Tlr_NonIgno_bd;
        double pval_Tlr_CC_sub_vh0Uh1,pval_Tw_CC_sub_vh0Uh1, pval_Tw_Igno_sub_vh0Uh1;
        double pval_Tlr_Igno_sub_vh0Uh1, pval_Tw_NonIgno_sub_vh0Uh1, pval_Tlr_NonIgno_sub_vh0Uh1;
        double pval_Tlr_CC_sub_vh0,pval_Tw_CC_sub_vh0, pval_Tw_Igno_sub_vh0;
        double pval_Tlr_Igno_sub_vh0, pval_Tw_NonIgno_sub_vh0, pval_Tlr_NonIgno_sub_vh0;
        double misratesXi1, misratesXi2, Tw_CC,Tlr_CC, Tw_Igno, Tlr_Igno, Tw_NonIgno, Tlr_NonIgno;
        double mu1_unc, mu2_unc, var1_unc, var2_unc, cov12_unc;
        double alpha01_unc, alpha11_unc, alpha02_unc, alpha12_unc;
        double mu1_con, mu2_con, var1_con, var2_con, cov12_con;
        double alpha01_con, alpha11_con, alpha02_con, alpha12_con;

        int SimuDone = 0, SimuSkipped = 0;
        int Rej_H0_counts_completeData = 0;
        int Rej_H0_counts_completeData_skippedInstances = 0;

        NonIgno_Avg_CovThetaHat_vh0Uh1 = 0; NonIgno_Avg_CovThetaHat_vh0 = 0;
        Igno_Avg_CovThetaHat_vh0Uh1 = 0; Igno_Avg_CovThetaHat_vh0 = 0;
        CC_Avg_CovThetaHat_vh0Uh1 = 0;  CC_Avg_CovThetaHat_vh0 = 0;


        otLog<<"pval_Tlr_CC_bd,pval_Tw_CC_bd,pval_Tw_Igno_bd,pval_Tlr_Igno_bd,\
        pval_Tw_NonIgno_bd,pval_Tlr_NonIgno_bd,pval_Tlr_CC_sub_vh0Uh1,pval_Tw_CC_sub_vh0Uh1,\
        pval_Tw_Igno_sub_vh0Uh1,pval_Tlr_Igno_sub_vh0Uh1, pval_Tw_NonIgno_sub_vh0Uh1,\
        pval_Tlr_NonIgno_sub_vh0Uh1,pval_Tlr_CC_sub_vh0,pval_Tw_CC_sub_vh0, pval_Tw_Igno_sub_vh0,\
        pval_Tlr_Igno_sub_vh0, pval_Tw_NonIgno_sub_vh0, pval_Tlr_NonIgno_sub_vh0,\
        misratesXi1, misratesXi2, Tw_CC,Tlr_CC, Tw_Igno, Tlr_Igno, Tw_NonIgno, Tlr_NonIgno,\
        mu1_unc, mu2_unc, var1_unc, var2_unc, cov12_unc,\
        alpha01_unc, alpha11_unc, alpha02_unc, alpha12_unc,\
        mu1_con, mu2_con, var1_con, var2_con, cov12_con,\
        alpha01_con, alpha11_con, alpha02_con, alpha12_con,\
        Seconds_Tw_Igno, Seconds_Tlr_Igno,Seconds_Tw_NonIgno,Seconds_Tlr_NonIgno"<<endl;

         dbmat L_scale(dim,dim);  dbcolvec mu_scale(dim);
         // standardized data by CCed sample mean/covariance matrix

        dbmat full_data_SamCov(dim,dim),full_data_CovMuHat_inv(dim,dim);
        dbcolvec full_data_SamMean(dim);
        double full_data_Tw, full_data_Tw_pval;

        double gradient_tolerance = 0.01;
        int num_alpha = 4; dbcolvec alpha_Vec_MNAR(num_alpha);
        double min_alpha,max_alpha;  // simulation instances with unusually large/small alpha (outside -8,8) would be skipped

        // SIMULATION LOOP
         while(SimuDone<NumSimu)
        {
            if( SimuSkipped>(0.4*NumSimu) && SimuDone>5 )
            {   osum<<"Simulation fails as more than 40% instances were dropped due to numerical errors \n";
                osumDts<<"Simulation fails as more than 40% instances were dropped due to numerical errors \n";
                return 0;    }

            //########### Simulate data and create missing data
            cout<<"Simulation Run No."<<(SimuDone+1)<<endl;
            osumDts<<"\n\n###################\n"\
                   <<"Simulation Run No."<<(SimuDone+1)<<"\n"<<"###################\n\n";

            // generate full sample
            mvrnorm(rndSample,Sz,dim,mu,Sig,Dat);

            // full data Tw
            SamMeanCov(Dat,Sz,dim, full_data_SamMean, full_data_SamCov);
            full_data_CovMuHat_inv = inv(full_data_SamCov)*Sz;
            get_Tw(full_data_Tw,dim,full_data_SamMean,full_data_CovMuHat_inv,StopCriteria);
            ChiBarWtsDim2ch2(full_data_SamCov,CbWts_vh0);
            full_data_Tw_pval = ChiBarPvalue(dim,CbWts_vh0,full_data_Tw);
            if(full_data_Tw_pval<0.05) Rej_H0_counts_completeData++;

            // Generate missing data according to logit{ prob xi1 mis }= alpha01 + alpha11*xi1;  logit{ prob xi2 mis }= alpha02 + alpha12*xi2 + alphaMAR*x_{i1}*(1-r_{i1})
            ZinExMod_GenerateMisInd(rndMisInd,Dat,alphaVec,MisInd,mu,alphaMAR,sd1,sd2);


            // ###############
            // complete-case  method:  unconstrained MLE/Cov are CCed sample mean/cov
            //#################
            NAomit_MisRates(Dat, MisInd, DatCCtemp, NumComCases, MissingRates);

            dbmat DatCCsam(NumComCases,dim);
            osumDts<<"Sample size after CC:\n"<<NumComCases<<"\n\n";
            DatCCsam = subm( DatCCtemp,range(0,NumComCases-1),range(0,dim-1) );
            osumDts<<"Missing rates (in %):\n"<<trans(MissingRates*100)<<"\n\n";

            osumDts<<"\n\n-----------Results from complete-case \n\n";
            SamMeanCov(DatCCsam,NumComCases,dim, SamMeanInit, SamCovInit);
            CC_L_unc = chol(SamCovInit);
            LToEta(CC_L_unc,dim,CC_Eta_unc);
            set_subm(CC_theta_unc,range(0,dim-1),range(0,0) )=SamMeanInit;
            set_subm(CC_theta_unc,range(dim,dim+dimEta-1),range(0,0) )=CC_Eta_unc;

            tStartTw = clock();
            pDat = &DatCCsam;  CC_MinNegLogLike_UnderH0UH1 = MLkOSamCom(CC_theta_unc);
            if (!is_finite(CC_MinNegLogLike_UnderH0UH1) )
            {   osumDts<<"Error in computing unconstrained MLE, simu instance skipped\n\n";
                SimuSkipped++;  continue;  }

            osumDts<<"Maximized log-likelihood under H0UH1: "<<(-1.0*CC_MinNegLogLike_UnderH0UH1)<<"\n";
            Seconds_unconstrained_MLE =  (double)(clock() - tStartTw)/CLOCKS_PER_SEC;


            osumDts<<"unconstrained MLE of theta (mu,Eta): "<<trans(CC_theta_unc)<<endl;
            osumDts<<"\nCC, unconstrained MLE of mean (i.e sample mean) \n"<<trans(SamMeanInit)<<endl;
            osumDts<<"CC, unconstrained MLE of Cov (i.e sample cov)\n"<<SamCovInit<<endl;

       // osumDts<<"Unconstrained MLE: "<<trans(CC_theta_unc)<<"\n";
       // osumDts<<"Max entry in unconstrained MLE: "<<dlib::max(CC_theta_unc)<<"\n";

            if( dlib::max(SamMeanInit) < 0.0001   ) // Constrained MLE = unconstrained MLE
            {
            osumDts<<"Since the unconstrained MLE already satisfies constraints (up to a numerical tolerance of 1e-3), it is also the constrained MLE. H0 won't be rejected.  \n\n";
      // Seconds_Tw = Seconds_unconstrained_MLE;
    //  osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw<<"\n";
     //    tStartTlr = clock(); CC_theta_con = CC_theta_unc;
    //    Seconds_constrained_MLE = (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
     //    Seconds_Tlr = Seconds_constrained_MLE + Seconds_unconstrained_MLE;
     //   osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr<<"\n";
     // osumDts<<"Tlr requires how much additional time (in \%): "<<( (Seconds_Tlr - Seconds_Tw)/Seconds_Tw )*100<<"% \n";

      // Seconds_Tw_CC_average  = Seconds_Tw_CC_average + Seconds_Tw;
     //  Seconds_Tlr_CC_average = Seconds_Tlr_CC_average + Seconds_Tlr;

            Tlr_CC = 0; pval_Tw_CC_bd = 1;  pval_Tlr_CC_bd = 1;
            // compute substitution p-value
            CC_CovThetaHat_this = inv( HESScdif(MLkOSamCom, CC_theta_unc) );
            CovMuHat_this_vh0Uh1 = subm(CC_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            CC_Avg_CovThetaHat_vh0 = CC_Avg_CovThetaHat_vh0 + CC_CovThetaHat_this;
            CC_Avg_CovThetaHat_vh0Uh1 = CC_Avg_CovThetaHat_vh0Uh1 + CC_CovThetaHat_this;
            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_CC_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_CC);
            pval_Tw_CC_sub_vh0 = pval_Tw_CC_sub_vh0Uh1;
            pval_Tlr_CC_sub_vh0 = pval_Tw_CC_sub_vh0Uh1;  pval_Tlr_CC_sub_vh0Uh1 = pval_Tw_CC_sub_vh0Uh1;
            }
         else
            {
            // step 2 of Tw
     //   tStartTw = clock();
            CC_CovThetaHat_this = inv( HESScdif(MLkOSamCom, CC_theta_unc) );
            CC_Avg_CovThetaHat_vh0Uh1 = CC_Avg_CovThetaHat_vh0Uh1 + CC_CovThetaHat_this;
            CovMuHat_this_vh0Uh1 = subm(CC_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0Uh1 = inv(CovMuHat_this_vh0Uh1);
            get_Tw(Tw_CC,dim,SamMeanInit,InvCovMuHat_this_vh0Uh1,StopCriteria);
        //    Seconds_Step2_Tw =  (double)(clock() - tStartTw)/CLOCKS_PER_SEC;
      //    Seconds_Tw = Seconds_unconstrained_MLE + Seconds_Step2_Tw;

        // for bd-constraint MLE: if initial value do not satisfy constraints, adjust to satisfy constraints
            CC_theta_con = CC_theta_unc;
              for(i=0; i<dim; ++i)  { if(CC_theta_con(i) > 0 ) CC_theta_con(i) = 0; }
                tStartTlr = clock();  MLkOSamComCtr = 0; pDat = &DatCCsam;
                BobyqaConvFail=0;  lb_Igno = -BIGNUM;
                ub_Igno = BIGNUM;  set_subm(ub_Igno,range(0,dim-1),range(0,0)) = 0;
  osumDts<<"\nbegin constrained MLE: \n";
  osumDts<<"CC_theta_con: "<<trans(CC_theta_con);
  osumDts<<"lb_Igno: "<<trans(lb_Igno);
  osumDts<<"ub_Igno: "<<trans(ub_Igno)<<endl;

         CC_MinNegLogLike_UnderH0 = find_min_bobyqaConvgFail(MLkOSamCom,
                                 CC_theta_con, npt_Igno,    // number of interpolation points
                                 lb_Igno,  // lower bound constraint
                                 ub_Igno,   // upper bound constraint
                                 10,    // initial trust region radius
                                 StopCriteria,  // stopping trust region radius
                                 MaxFunEvalLong, // max number of objective function evaluations
                                 BobyqaConvFail);

                if (!is_finite(CC_MinNegLogLike_UnderH0) || BobyqaConvFail==1)
                { osumDts<<"Error in computing constrained MLE, simu instance skipped\n\n";
                    SimuSkipped++;  continue;  }

            osumDts<<"Maximized log-likelihood under H0: "<<(-CC_MinNegLogLike_UnderH0)<<"\n";
     //   Seconds_constrained_MLE =   (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
    //   tStartTlr = clock();
            Tlr_CC = 2*(CC_MinNegLogLike_UnderH0 - CC_MinNegLogLike_UnderH0UH1);
            if(Tlr_CC < 0)
             {  osumDts<<"!warning: maximized log-like under H0 is unexpectedly greater than that under H0UH1, Tlr adjusted to 0.\n";   Tlr_CC = 0; }
     //   Seconds_Tlr = Seconds_constrained_MLE + Seconds_unconstrained_MLE  + (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;

            CC_mu_con = subm(CC_theta_con,range(0,dim-1),range(0,0) );
            CC_Eta_con = subm(CC_theta_con,range(dim,ParaMuSig-1),range(0,0) );
            EtaToL(CC_Eta_con,dim,CC_L_con);

            osumDts<<"constrained MLE of theta (mu,Eta): "<<trans(CC_theta_con)<<endl;
            osumDts<<"\nCC, bd-constrained MLE of mean\n"<<trans(CC_mu_con)<<endl;
            osumDts<<"CC, bd-constrained MLE of Cov\n"<<CC_L_con*trans(CC_L_con)<<endl;
           // osumDts<<"\n#(loglike evaluation during numerical optimization): "<<MLkOSamComCtr<<endl<<endl;

       //     osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw<<"\n";
       //     osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr<<"\n";
      //      osumDts<<"Tlr requires how much additional time (in \%): "<<( (Seconds_Tlr - Seconds_Tw)/Seconds_Tw )*100<<"% \n";

            CC_CovThetaHat_this = inv( HESScdif(MLkOSamCom, CC_theta_con) );
            CovMuHat_this_vh0 = subm(CC_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0 = inv(CovMuHat_this_vh0);
            osumDts<<"\nCov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
            osumDts<<CovMuHat_this_vh0<<"\n";

            osumDts<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
            osumDts<<CovMuHat_this_vh0Uh1<<"\n";

            osumDts<<"Cov(MuHat) reference: Sigma/SampleSize; it is exact when Muhat is the complete-sample mean \n";
            osumDts<<(Sig/Sz)<<"\n\n";

            osumDts<<"Test results based on Tw:\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            osumDts<<"Tw = "<<Tw_CC<<"\n";
            pval_Tw_CC_bd = 1-(scythe::pchisq(Tw_CC, double(dim-1)) \
            + scythe::pchisq(Tw_CC, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tw: "<<pval_Tw_CC_bd<<"\n";

            // substitution p-value
            R_ThetaZeroInv_Rt_vh0 =  Sz*CovMuHat_this_vh0;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0,CbWts_vh0);
            pval_Tw_CC_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tw_CC);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tw_CC_sub_vh0<<"\n";

            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_CC_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_CC);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tw_CC_sub_vh0Uh1<<"\n";

            osumDts<<"\n\n";
            osumDts<<"Test results based on Tlr:\n";
            osumDts<<"Tlr = "<<Tlr_CC<<"\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            pval_Tlr_CC_bd = 1-(scythe::pchisq(Tlr_CC, double(dim-1)) \
            + scythe::pchisq(Tlr_CC, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tlr: "<<pval_Tlr_CC_bd<<"\n";

            // substitution p-value
            pval_Tlr_CC_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tlr_CC);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tlr_CC_sub_vh0<<"\n";

            pval_Tlr_CC_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tlr_CC);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tlr_CC_sub_vh0Uh1<<"\n";
            }


            // ###############
            // Assuming ignorable missingness
            //#################
            osumDts<<"\n\n-----------Results from assuming ignorable missingness \n\n";
            Igno_theta_initial_unc =  CC_theta_unc;

            tStartTw = clock();
            pDat=&Dat; PMisInd=&MisInd;
            BobyqaConvFail=0; lb_Igno = -BIGNUM; ub_Igno = BIGNUM;
            Igno_theta_unc = Igno_theta_initial_unc;
            Igno_MinNegLogLike_UnderH0UH1 = find_min_bobyqaConvgFail(MLkOSamMAR,
            Igno_theta_unc, npt_Igno,    // number of interpolation points
            lb_Igno,  // lower bound constraint
            ub_Igno,   // upper bound constraint
            8,    // initial trust region radius
            StopCriteria*0.1,  // stopping trust region radius
            MaxFunEvalLong, // max number of objective function evaluations
BobyqaConvFail
            );

          Igno_Gradient = MLkOSamMARGrad(Igno_theta_unc);
          GetMinMax(Igno_Gradient, MinGrad, MaxGrad);
         // bobyqa converges  if gradient at its solution is not greater than 0.5 in absolute value
        if(is_finite(Igno_MinNegLogLike_UnderH0UH1) && BobyqaConvFail==0 \
           && MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )
        goto Unconstrained_MLE_Success;

        osumDts<<"non-convergence occurs in bobyqa, backup option NEWTON activated:\n";
        NewtonConvFail = 0; Igno_theta_unc = Igno_theta_initial_unc;
            Igno_MinNegLogLike_UnderH0UH1
            =find_minConvFail(newton_search_strategy(MLkOSamMARHessian),
            objective_delta_stop_strategy(StopCriteria), //.be_verbose()
            MLkOSamMAR,
            MLkOSamMARGrad,
            Igno_theta_unc,
            -BIGNUM,NewtonConvFail,MaxIterations);
        if(is_finite(Igno_MinNegLogLike_UnderH0UH1) && NewtonConvFail==0)
        {
             Igno_Gradient = MLkOSamMARGrad(Igno_theta_unc);
            GetMinMax(Igno_Gradient, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success;
        }

        osumDts<<"non-convergence occurs in NEWTON, backup option TRUST-REGION activated:\n";
        TrustRegionConvFail = 0; Igno_theta_unc = Igno_theta_initial_unc;
        Igno_MinNegLogLike_UnderH0UH1 = \
        find_min_trust_regionConvFail(objective_delta_stop_strategy(StopCriteriaTrustRegion),\
D2MAR_model(), Igno_theta_unc, TrustRegionConvFail,MaxIterations,10 );
        if(is_finite(Igno_MinNegLogLike_UnderH0UH1) && TrustRegionConvFail==0)
          {
             Igno_Gradient = MLkOSamMARGrad(Igno_theta_unc);
            GetMinMax(Igno_Gradient, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success;
        }

        osumDts<<"non-convergence occurs in TRUST-REGION, backup option BFGS activated:\n";
        bfgsConvFail = 0; Igno_theta_unc = Igno_theta_initial_unc;
        Igno_MinNegLogLike_UnderH0UH1 = \
        find_minConvFail(bfgs_search_strategy(),
           objective_delta_stop_strategy(StopCriteria),
           MLkOSamMAR, MLkOSamMARGrad, Igno_theta_unc, -BIGNUM,bfgsConvFail,MaxIterations);
        if(is_finite(Igno_MinNegLogLike_UnderH0UH1) && bfgsConvFail==0)
          {
             Igno_Gradient = MLkOSamMARGrad(Igno_theta_unc);
            GetMinMax(Igno_Gradient, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success;
        }

        osumDts<<"calculation of unconstrained MLE encounters non-convergence \
        after all options tried, skip this simu instance\n";
        SimuSkipped++; continue;

  Unconstrained_MLE_Success:
            Seconds_unconstrained_MLE = (double)(clock() - tStartTw)/CLOCKS_PER_SEC;
      //osumDts<<"Seconds_unconstrained_MLE: "<<Seconds_unconstrained_MLE<<endl;
            Igno_Gradient = MLkOSamMARGrad(Igno_theta_unc);
            GetMinMax(Igno_Gradient,MinGrad,MaxGrad);
            osumDts<<"\nmin/max value of the gradient at unconstrained MLE: "<<MinGrad<<"/"<<MaxGrad<<"\n";
            Igno_mu_unc = subm(Igno_theta_unc, range(0,dim-1), range(0,0) );
            Igno_Eta_unc = subm(Igno_theta_unc, range(dim,ParaMuSig-1), range(0,0) );
            EtaToL(Igno_Eta_unc,dim,Igno_L_unc);
            osumDts<<"Maximized log-likelihood under H0UH1: "\
            <<(-1.0)*Igno_MinNegLogLike_UnderH0UH1<<"\n";
            osumDts<<"unconstrained MLE of theta (mu,Eta): "<<trans(Igno_theta_unc);
            osumDts<<"Assuming ignorable missingness, unconstrained MLE of mean \n"\
            <<trans(Igno_mu_unc)<<endl;
            osumDts<<"Assuming ignorable missingness, unconstrained MLE of Cov \n"\
            <<( Igno_L_unc*trans(Igno_L_unc) )<<endl;

      if( dlib::max(Igno_mu_unc) < 0.0001   ) // Constrained MLE = unconstrained MLE
            {
            osumDts<<"Since the unconstrained MLE already satisfies constraints (up to a numerical tolerance of 1e-3), it is also the constrained MLE. H0 won't be rejected.  \n\n";

 // osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw<<"\n";
           Seconds_Tw_Igno = Seconds_unconstrained_MLE;
            tStartTlr = clock(); Igno_theta_con = Igno_theta_unc;
            Seconds_constrained_MLE = (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
     //osumDts<<"Seconds_constrained_MLE: "<<Seconds_constrained_MLE<<endl;
            Seconds_Tlr_Igno = Seconds_constrained_MLE + Seconds_unconstrained_MLE;
            osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw_Igno<<"\n";
            osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr_Igno<<"\n";
            osumDts<<"Tlr requires how much additional time (in \%): "\
            <<( (Seconds_Tlr_Igno - Seconds_Tw_Igno)/Seconds_Tw_Igno )*100<<"% \n";

     //  Seconds_Tw_Igno_average  = Seconds_Tw_Igno_average + Seconds_Tw;
    //  Seconds_Tlr_Igno_average = Seconds_Tlr_Igno_average + Seconds_Tlr;

            Tw_Igno = 0; Tlr_Igno = 0; pval_Tw_Igno_bd = 1;  pval_Tlr_Igno_bd = 1;
            // compute substitution p-value
            Igno_CovThetaHat_this = inv( HESScdif(MLkOSamMAR, Igno_theta_unc) );
            Igno_Avg_CovThetaHat_vh0 = Igno_Avg_CovThetaHat_vh0 + Igno_CovThetaHat_this;
            Igno_Avg_CovThetaHat_vh0Uh1 = Igno_Avg_CovThetaHat_vh0Uh1 + Igno_CovThetaHat_this;

            CovMuHat_this_vh0Uh1 = subm(Igno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0Uh1 = inv(CovMuHat_this_vh0Uh1);
            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_Igno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_Igno);
            pval_Tw_Igno_sub_vh0 = pval_Tw_Igno_sub_vh0Uh1;
            pval_Tlr_Igno_sub_vh0 = pval_Tw_Igno_sub_vh0Uh1;  pval_Tlr_Igno_sub_vh0Uh1 = pval_Tw_Igno_sub_vh0Uh1;
            }

      else
            {
            // step 2 of Tw
            tStartTw = clock();
            Igno_CovThetaHat_this = inv( HESScdif(MLkOSamMAR, Igno_theta_unc) );
            Igno_Avg_CovThetaHat_vh0Uh1 = Igno_Avg_CovThetaHat_vh0Uh1 + Igno_CovThetaHat_this;

            CovMuHat_this_vh0Uh1 = subm(Igno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0Uh1 = inv(CovMuHat_this_vh0Uh1);
            Igno_mu_unc = subm(Igno_theta_unc, range(0,dim-1), range(0,0) );
            get_Tw(Tw_Igno,dim,Igno_mu_unc,InvCovMuHat_this_vh0Uh1,StopCriteria);
            Seconds_Step2_Tw =  (double)(clock() - tStartTw)/CLOCKS_PER_SEC;
     // osumDts<<"Seconds_Step2_Tw: "<<Seconds_Step2_Tw<<endl;
            Seconds_Tw_Igno = Seconds_unconstrained_MLE + Seconds_Step2_Tw;

            Igno_theta_initial_con = Igno_theta_unc;
            for(i = 0; i< dim; ++i)
            { if(Igno_theta_initial_con(i)>0)  Igno_theta_initial_con(i) = 0; }

            tStartTlr = clock();   pDat = &Dat; PMisInd = &MisInd;
            BobyqaConvFail = 0;  lb_Igno = -BIGNUM;
            ub_Igno = BIGNUM;  set_subm(ub_Igno,range(0,dim-1),range(0,0)) = 0;
            Igno_theta_con = Igno_theta_initial_con;

            osumDts<<"\nbegin constrained MLE: \n";
            osumDts<<"Igno_theta_initial_con: "<<trans(Igno_theta_initial_con);
            osumDts<<"lb_Igno: "<<trans(lb_Igno);
            osumDts<<"ub_Igno: "<<trans(ub_Igno)<<endl;

           if( !(dlib::min(Igno_theta_con - lb_Igno) >= 0 \
               && dlib::max(ub_Igno - Igno_theta_con) >= 0) )
            {    osumDts<<"calculation of constrained MLE fails: skip\n\n";
               SimuSkipped++; continue;   }

                Igno_MinNegLogLike_UnderH0 = find_min_bobyqaConvgFail(MLkOSamMAR,
                 Igno_theta_con, npt_Igno,    // number of interpolation points
                  lb_Igno,  // lower bound constraint
                  ub_Igno,   // upper bound constraint
                   10,    // initial trust region radius
                   StopCriteria,  // stopping trust region radius
                   MaxFunEvalLong, // max number of objective function evaluations
                   BobyqaConvFail);

                if (is_finite(Igno_MinNegLogLike_UnderH0) && BobyqaConvFail==0)
                    goto ConstrainedMLE_success;

        osumDts<<"non-convergence occurs in BOBYQA, backup option BFGS activated:\n";
        bfgsConvFail = 0;  pDat=&Dat; PMisInd=&MisInd;
        Igno_theta_con = Igno_theta_initial_con;
        Igno_MinNegLogLike_UnderH0 = \
        find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria), \
            MLkOSamMAR, MLkOSamMARGrad, Igno_theta_con, lb_Igno,ub_Igno,bfgsConvFail,MaxIterations);
        if(is_finite(Igno_MinNegLogLike_UnderH0) && bfgsConvFail==0)
        goto ConstrainedMLE_success;

                 osumDts<<"Error in computing constrained MLE, simu instance skipped\n\n";
                    SimuSkipped++; continue;


     ConstrainedMLE_success:
            osumDts<<"Maximized log-likelihood under H0: "<<(-Igno_MinNegLogLike_UnderH0)<<"\n";
            Seconds_constrained_MLE =   (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
 //osumDts<<"Seconds_constrained_MLE: "<<Seconds_constrained_MLE<<endl;
            tStartTlr = clock();
            Tlr_Igno = 2*(Igno_MinNegLogLike_UnderH0 - Igno_MinNegLogLike_UnderH0UH1);
            if(Tlr_Igno < 0)
             {  osumDts<<"!warning: maximized log-like under H0 is unexpectedly greater than that under H0UH1, Tlr adjusted to 0.\n";   Tlr_Igno = 0; }
            Seconds_Tlr_Igno = Seconds_constrained_MLE + Seconds_unconstrained_MLE \
            + (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;

            Igno_mu_con = subm(Igno_theta_con,range(0,dim-1),range(0,0) );
            Igno_Eta_con = subm(Igno_theta_con,range(dim,ParaMuSig-1),range(0,0) );
            EtaToL(Igno_Eta_con,dim,Igno_L_con);
            osumDts<<"constrained MLE of theta (mu,Eta): "<<trans(Igno_theta_con);
            osumDts<<"Assuming ignorable missingness, bd-constrained MLE of mean\n"<<trans(Igno_mu_con)<<endl;
            osumDts<<"Assuming ignorable missingness, bd-constrained MLE of Cov\n"<<Igno_L_con*trans(Igno_L_con)<<endl;
           // osumDts<<"\n#(loglike evaluation during numerical optimization): "<<MLkOSamComCtr<<endl<<endl;

            osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw_Igno<<"\n";
            osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr_Igno<<"\n";
            osumDts<<"Tlr requires how much additional time (in \%): "\
            <<( (Seconds_Tlr_Igno - Seconds_Tw_Igno)/Seconds_Tw_Igno )*100<<"% \n";

            Igno_CovThetaHat_this = inv( HESScdif(MLkOSamMAR, Igno_theta_con) );
            Igno_Avg_CovThetaHat_vh0 = Igno_Avg_CovThetaHat_vh0 + Igno_CovThetaHat_this;
            CovMuHat_this_vh0 = subm(Igno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0 = inv(CovMuHat_this_vh0);
            osumDts<<"\nCov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
            osumDts<<CovMuHat_this_vh0<<"\n";

            osumDts<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
            osumDts<<CovMuHat_this_vh0Uh1<<"\n";

            osumDts<<"Cov(MuHat) reference: Sigma/SampleSize; it is exact when Muhat is the complete-sample mean \n";
            osumDts<<(Sig/Sz)<<"\n\n";

            osumDts<<"Test results based on Tw:\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            osumDts<<"Tw = "<<Tw_Igno<<"\n";
            pval_Tw_Igno_bd = 1-(scythe::pchisq(Tw_Igno, double(dim-1)) \
            + scythe::pchisq(Tw_Igno, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tw: "<<pval_Tw_Igno_bd<<"\n";

            // substitution p-value
            R_ThetaZeroInv_Rt_vh0 =  Sz*CovMuHat_this_vh0;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0,CbWts_vh0);
            pval_Tw_Igno_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tw_Igno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tw_Igno_sub_vh0<<"\n";

            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_Igno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_Igno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tw_Igno_sub_vh0Uh1<<"\n";

            osumDts<<"\n\n";
            osumDts<<"Test results based on Tlr:\n";
            osumDts<<"Tlr = "<<Tlr_Igno<<"\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            pval_Tlr_Igno_bd = 1-(scythe::pchisq(Tlr_Igno, double(dim-1)) \
            + scythe::pchisq(Tlr_Igno, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tlr: "<<pval_Tlr_Igno_bd<<"\n";

            // substitution p-value
            pval_Tlr_Igno_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tlr_Igno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tlr_Igno_sub_vh0<<"\n";

            pval_Tlr_Igno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tlr_Igno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tlr_Igno_sub_vh0Uh1<<"\n";
            }

             // ###############
            // Assuming non-ignorable missingness:
            // logit{Pr(xi1 is missing)} = alpha01 + alpha11*xi1;
            // logit{Pr(xi2 is missing)} = alpha02 + alpha12*xi2;
            //#################


            osumDts<<"\n\n-----------Results from assuming a non-ignorable model:\n \
logit{Pr(xi1 is missing)} = alpha01 + alpha11*(xi1-mu1)/sig1;\n \
logit{Pr(xi2 is missing)} = alpha02 + alpha12*(xi2-mu2)/sig2; \n\n";
            set_subm(NonIgno_theta_initial_unc,range(0,ParaMuSig-1),range(0,0)) =  Igno_theta_unc;

            NonIgno_theta_initial_unc (6) = 0; NonIgno_theta_initial_unc (8) = 0;
            // initial values for alpha11, alpha12 are 0

            NonIgno_theta_initial_unc (5) = std::log( MissingRates(0)/ ( 1-MissingRates(0) ) );
            NonIgno_theta_initial_unc (7) = std::log( MissingRates(1)/ ( 1-MissingRates(1) ) );
            // initial values for alpha01, alpha02 are logit missing rates

            /*
            int Rej_H0_counts_completeData = 0;
            try BFGS first,
            if  {  BFGS converges }
            if( gradient within (-0.01,0.01) goto unconstrained_MLE_obtained;
             else  { use the MLE from BFSG as initial value}

              }
             else { initial value as usual}

              try bobyqa, newton and trust region

            if (all the three fails_
            {
             calculate the complete-data based Tw, carry test,
             if(H0_rejecetd)  Rej_H0_counts_completeData++;
              skip_counts++;
            }

              unconstrained_MLE_obtained:

              after simu loop,   cout<<Rej_H0_counts_completeData;  cout<<skip_counts;
             */

            tStartTw = clock();


            pDat = &Dat; PMisInd = &MisInd;  pQHPts = &QHPts;
            lb_NonIgno = -BIGNUM; ub_NonIgno = BIGNUM;

            NonIgno_theta_unc = NonIgno_theta_initial_unc;
            int BFGSFail = 0;
            NonIgno_MinNegLogLike_UnderH0UH1 = find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria), \
            mod_NoRi1_MNARd2_MLogLikGH, mod_NoRi1_MNARd2_MLogLikGH_Grad, NonIgno_theta_unc, lb_NonIgno,ub_NonIgno,BFGSFail,MaxIterations);

            if( BFGSFail == 0 )
            {
                NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
          GetMinMax(NonIgno_Gradient, MinGrad, MaxGrad);
                if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )
                      {
                        osumDts<<"unconstrained MLE successfully obtained from BFGS\n";
                        goto Unconstrained_MLE_Success_NonIgno;}
                else  {
                        osumDts<<"BFGS converges but gradients are too large, BFGS solutions are re-used as new initial values for newton, bobyqa and trust_region methods \n";
                        NonIgno_theta_initial_unc_bfgs = NonIgno_theta_unc;  }
            }

      BobyqaConvFail = 0;
            NonIgno_theta_unc = NonIgno_theta_initial_unc_bfgs;
            NonIgno_MinNegLogLike_UnderH0UH1 = find_min_bobyqaConvgFail(mod_NoRi1_MNARd2_MLogLikGH,
            NonIgno_theta_unc, npt_NonIgno,    // number of interpolation points
            lb_NonIgno,  // lower bound constraint
            ub_NonIgno,   // upper bound constraint
            8,    // initial trust region radius
            StopCriteria*0.1,  // stopping trust region radius
            MaxFunEvalLong, // max number of objective function evaluations
            BobyqaConvFail
            );

          NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
          GetMinMax(NonIgno_Gradient, MinGrad, MaxGrad);
         // bobyqa converges  if gradient at its solution is not greater than 0.5 in absolute value
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && BobyqaConvFail==0 \
           && MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )
        goto Unconstrained_MLE_Success_NonIgno;

        osumDts<<"non-convergence occurs in bobyqa, backup option NEWTON activated:\n";
        NewtonConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc_bfgs;
            NonIgno_MinNegLogLike_UnderH0UH1
            =find_minConvFail(newton_search_strategy(mod_NoRi1_MNARd2_MLogLikGH_Hess),
            objective_delta_stop_strategy(StopCriteria), //.be_verbose()
            mod_NoRi1_MNARd2_MLogLikGH,
            mod_NoRi1_MNARd2_MLogLikGH_Grad,
            NonIgno_theta_unc,
            -BIGNUM,NewtonConvFail,MaxIterations);
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && NewtonConvFail==0)
               {
             NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }


        osumDts<<"non-convergence occurs in NEWTON, backup option TRUST-REGION activated:\n";
        TrustRegionConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc_bfgs;
        NonIgno_MinNegLogLike_UnderH0UH1 = \
        find_min_trust_regionConvFail(objective_delta_stop_strategy(StopCriteriaTrustRegion),\
mod_NoRi1_MNARd2_MLogLikGH_class(), NonIgno_theta_unc, TrustRegionConvFail,MaxIterations,10 );
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && TrustRegionConvFail==0)
        {
            NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }

        osumDts<<"non-convergence occurs in TRUST-REGION, backup option BFGS activated:\n";
        bfgsConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc_bfgs;
        NonIgno_MinNegLogLike_UnderH0UH1 = \
        find_minConvFail(bfgs_search_strategy(),
           objective_delta_stop_strategy(StopCriteria),
           mod_NoRi1_MNARd2_MLogLikGH, mod_NoRi1_MNARd2_MLogLikGH_Grad, NonIgno_theta_unc, -BIGNUM,bfgsConvFail,MaxIterations);
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && bfgsConvFail==0)
        {
            NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }

            osumDts<<"Switch back to original initial values, try bfgs, newton, trust_region method again\n";
            BobyqaConvFail = 0;
            NonIgno_theta_unc = NonIgno_theta_initial_unc;
            NonIgno_MinNegLogLike_UnderH0UH1 = find_min_bobyqaConvgFail(mod_NoRi1_MNARd2_MLogLikGH,
            NonIgno_theta_unc, npt_NonIgno,    // number of interpolation points
            lb_NonIgno,  // lower bound constraint
            ub_NonIgno,   // upper bound constraint
            8,    // initial trust region radius
            StopCriteria*0.1,  // stopping trust region radius
            MaxFunEvalLong, // max number of objective function evaluations
            BobyqaConvFail
            );

          NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
          GetMinMax(NonIgno_Gradient, MinGrad, MaxGrad);
         // bobyqa converges  if gradient at its solution is not greater than 0.5 in absolute value
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && BobyqaConvFail==0 \
           && MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )
        goto Unconstrained_MLE_Success_NonIgno;

        osumDts<<"non-convergence occurs in bobyqa, backup option NEWTON activated:\n";
        NewtonConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc;
            NonIgno_MinNegLogLike_UnderH0UH1
            =find_minConvFail(newton_search_strategy(mod_NoRi1_MNARd2_MLogLikGH_Hess),
            objective_delta_stop_strategy(StopCriteria), //.be_verbose()
            mod_NoRi1_MNARd2_MLogLikGH,
            mod_NoRi1_MNARd2_MLogLikGH_Grad,
            NonIgno_theta_unc,
            -BIGNUM,NewtonConvFail,MaxIterations);
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && NewtonConvFail==0)
               {
             NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }


        osumDts<<"non-convergence occurs in NEWTON, backup option TRUST-REGION activated:\n";
        TrustRegionConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc;
        NonIgno_MinNegLogLike_UnderH0UH1 = \
        find_min_trust_regionConvFail(objective_delta_stop_strategy(StopCriteriaTrustRegion),\
mod_NoRi1_MNARd2_MLogLikGH_class(), NonIgno_theta_unc, TrustRegionConvFail,MaxIterations,10 );
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && TrustRegionConvFail==0)
        {
            NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }

        osumDts<<"non-convergence occurs in TRUST-REGION, backup option BFGS activated:\n";
        bfgsConvFail = 0; NonIgno_theta_unc = NonIgno_theta_initial_unc;
        NonIgno_MinNegLogLike_UnderH0UH1 = \
        find_minConvFail(bfgs_search_strategy(),
           objective_delta_stop_strategy(StopCriteria),
           mod_NoRi1_MNARd2_MLogLikGH, mod_NoRi1_MNARd2_MLogLikGH_Grad, NonIgno_theta_unc, -BIGNUM,bfgsConvFail,MaxIterations);
        if(is_finite(NonIgno_MinNegLogLike_UnderH0UH1) && bfgsConvFail==0)
        {
            NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_theta_unc, MinGrad, MaxGrad);
          if( MaxGrad < gradient_tolerance && MinGrad > (-1.0*gradient_tolerance) )  goto Unconstrained_MLE_Success_NonIgno;
        }

        osumDts<<"calculation of unconstrained MLE encounters non-convergence \
        after all options tried, skip this simu instance\n";
        SimuSkipped++; continue;

  Unconstrained_MLE_Success_NonIgno:
      alpha_Vec_MNAR = subm(NonIgno_theta_unc,range(5,8),range(0,0));
             GetMinMax(alpha_Vec_MNAR,min_alpha,max_alpha);
      if(min_alpha <  (-8) || max_alpha > 8 )
      {
       osumDts<<"Min/Max of MLE for missing data model parameters: "<<min_alpha <<"/"<<max_alpha<<endl;
       osumDts<<"MLE of missing data model parameters are unusual, simulation instance skipped\n";
       SimuSkipped++; continue;
      }

            Seconds_unconstrained_MLE = (double)(clock() - tStartTw)/CLOCKS_PER_SEC;
      //osumDts<<"Seconds_unconstrained_MLE: "<<Seconds_unconstrained_MLE<<endl;
            NonIgno_Gradient = mod_NoRi1_MNARd2_MLogLikGH_Grad(NonIgno_theta_unc);
            GetMinMax(NonIgno_Gradient,MinGrad,MaxGrad);
            osumDts<<"\nmin/max value of the gradient at unconstrained MLE: "<<MinGrad<<"/"<<MaxGrad<<"\n";
            NonIgno_mu_unc = subm(NonIgno_theta_unc, range(0,dim-1), range(0,0) );
            NonIgno_Eta_unc = subm(NonIgno_theta_unc, range(dim,ParaMuSig-1), range(0,0) );
            EtaToL(NonIgno_Eta_unc,dim,NonIgno_L_unc);
            osumDts<<"Maximized log-likelihood under H0UH1: "\
            <<(-1.0)*NonIgno_MinNegLogLike_UnderH0UH1<<"\n\n";
            osumDts<<"unconstrained MLE of theta (mu,Eta,alpha): "<<trans(NonIgno_theta_unc);
            osumDts<<"Assuming the non-ignorable model, unconstrained MLE of mean \n"\
            <<trans(NonIgno_mu_unc)<<endl;
            NonIgno_Sig_unc = NonIgno_L_unc*trans(NonIgno_L_unc);
            osumDts<<"Assuming the non-ignorable model, unconstrained MLE of Cov \n"\
            <<(NonIgno_Sig_unc)<<endl;

      if( dlib::max(NonIgno_mu_unc) < 0.0001   ) // Constrained MLE = unconstrained MLE
            {
            osumDts<<"Since the unconstrained MLE already satisfies constraints (up to a numerical tolerance of 1e-3), it is also the constrained MLE. H0 won't be rejected.  \n\n";

 // osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw<<"\n";
            Seconds_Tw_NonIgno = Seconds_unconstrained_MLE;
            tStartTlr = clock(); NonIgno_theta_con = NonIgno_theta_unc;
            NonIgno_Sig_con = NonIgno_Sig_unc;
            Seconds_constrained_MLE = (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
     //osumDts<<"Seconds_constrained_MLE: "<<Seconds_constrained_MLE<<endl;
            Seconds_Tlr_NonIgno = Seconds_constrained_MLE + Seconds_unconstrained_MLE;
            osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr_NonIgno<<"\n";
            osumDts<<"Tlr requires how much additional time (in \%): "\
            <<( (Seconds_Tlr_NonIgno - Seconds_Tw_NonIgno)/Seconds_Tw_NonIgno )*100<<"% \n";

  //  Seconds_Tw_NonIgno_average  = Seconds_Tw_NonIgno_average + Seconds_Tw;
  // Seconds_Tlr_NonIgno_average = Seconds_Tlr_NonIgno_average + Seconds_Tlr;

            Tw_NonIgno = 0; Tlr_NonIgno = 0; pval_Tw_NonIgno_bd = 1;  pval_Tlr_NonIgno_bd = 1;
            // compute substitution p-value
            NonIgno_CovThetaHat_this = inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, NonIgno_theta_unc) );
            NonIgno_Avg_CovThetaHat_vh0 = NonIgno_Avg_CovThetaHat_vh0 + NonIgno_CovThetaHat_this;
            NonIgno_Avg_CovThetaHat_vh0Uh1 = NonIgno_Avg_CovThetaHat_vh0Uh1 \
            + NonIgno_CovThetaHat_this;

            CovMuHat_this_vh0Uh1 = subm(NonIgno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0Uh1 = inv(CovMuHat_this_vh0Uh1);
            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_NonIgno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_NonIgno);
            pval_Tw_NonIgno_sub_vh0 = pval_Tw_NonIgno_sub_vh0Uh1;
            pval_Tlr_NonIgno_sub_vh0 = pval_Tw_NonIgno_sub_vh0Uh1;
            pval_Tlr_NonIgno_sub_vh0Uh1 = pval_Tw_NonIgno_sub_vh0Uh1;
            }

      else
            {
            // step 2 of Tw
            tStartTw = clock();
            NonIgno_CovThetaHat_this = inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, NonIgno_theta_unc) );
            NonIgno_Avg_CovThetaHat_vh0Uh1 = NonIgno_Avg_CovThetaHat_vh0Uh1 \
            + NonIgno_CovThetaHat_this;

            CovMuHat_this_vh0Uh1 = subm(NonIgno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0Uh1 = inv(CovMuHat_this_vh0Uh1);
            NonIgno_mu_unc = subm(NonIgno_theta_unc, range(0,dim-1), range(0,0) );
            get_Tw(Tw_NonIgno,dim,NonIgno_mu_unc,InvCovMuHat_this_vh0Uh1,StopCriteria);
            Seconds_Step2_Tw =  (double)(clock() - tStartTw)/CLOCKS_PER_SEC;
     // osumDts<<"Seconds_Step2_Tw: "<<Seconds_Step2_Tw<<endl;
            Seconds_Tw_NonIgno = Seconds_unconstrained_MLE + Seconds_Step2_Tw;

            NonIgno_theta_initial_con = NonIgno_theta_initial_unc;
            for(i=0; i<dim; ++i)
            if(NonIgno_theta_initial_con(i) > 0 ) NonIgno_theta_initial_con(i)=0;

            tStartTlr = clock();   pDat = &Dat; PMisInd = &MisInd; pQHPts=&QHPts;
            BobyqaConvFail = 0;  lb_NonIgno = -BIGNUM;
            ub_NonIgno = BIGNUM;  set_subm(ub_NonIgno,range(0,dim-1),range(0,0)) = 0;
            NonIgno_theta_con = NonIgno_theta_initial_con;


           osumDts<<"\nbegin constrained MLE: \n";
           osumDts<<"NonIgno_theta_initial_con: "<<trans(NonIgno_theta_initial_con);
           osumDts<<"lb_NonIgno: "<<trans(lb_NonIgno);
           osumDts<<"ub_NonIgno: "<<trans(ub_NonIgno)<<endl;

              if( !(dlib::min(NonIgno_theta_initial_con - lb_NonIgno) >= 0 \
               && dlib::max(ub_NonIgno - NonIgno_theta_initial_con) >= 0) )
         {     osumDts<<"calculation of constrained MLE fails: skip\n\n";
               SimuSkipped++; continue;
         }

            NonIgno_MinNegLogLike_UnderH0 = find_min_bobyqaConvgFail(mod_NoRi1_MNARd2_MLogLikGH,
            NonIgno_theta_con, npt_NonIgno,    // number of interpolation points
            lb_NonIgno,  // lower bound constraint
            ub_NonIgno,   // upper bound constraint
            10,    // initial trust region radius
            StopCriteria,  // stopping trust region radius
            MaxFunEvalLong, // max number of objective function evaluations
            BobyqaConvFail);

                if (is_finite(NonIgno_MinNegLogLike_UnderH0) && BobyqaConvFail==0)
                    goto ConstrainedMLE_success_NonIgno;

        osumDts<<"non-convergence occurs in BOBYQA, backup option BFGS activated:\n";
        bfgsConvFail = 0;  pDat=&Dat; PMisInd=&MisInd; pQHPts=&QHPts;
        NonIgno_theta_con = NonIgno_theta_initial_con;
        NonIgno_MinNegLogLike_UnderH0 = \
        find_min_box_constrainedConvFail(bfgs_search_strategy(),\
            objective_delta_stop_strategy(StopCriteria), \
            MLkOSamMAR, MLkOSamMARGrad, NonIgno_theta_con, \
            lb_NonIgno,ub_NonIgno,bfgsConvFail,MaxIterations);
        if(is_finite(NonIgno_MinNegLogLike_UnderH0) && bfgsConvFail==0)
        goto ConstrainedMLE_success_NonIgno;

                 osumDts<<"Error in computing constrained MLE, simu instance skipped\n\n";
                    SimuSkipped++; continue;

     ConstrainedMLE_success_NonIgno:
            osumDts<<"Maximized log-likelihood under H0: "<<(-NonIgno_MinNegLogLike_UnderH0)<<"\n";
            Seconds_constrained_MLE =  (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;
 //osumDts<<"Seconds_constrained_MLE: "<<Seconds_constrained_MLE<<endl;
            tStartTlr = clock();
            Tlr_NonIgno = 2*(NonIgno_MinNegLogLike_UnderH0 - NonIgno_MinNegLogLike_UnderH0UH1);
            if(Tlr_NonIgno < 0)
             {  osumDts<<"!warning: maximized log-like under H0 is unexpectedly greater than that under H0UH1, Tlr adjusted to 0.\n";   Tlr_NonIgno = 0; }
            Seconds_Tlr_NonIgno = Seconds_constrained_MLE + Seconds_unconstrained_MLE \
            + (double)(clock() - tStartTlr)/CLOCKS_PER_SEC;

            NonIgno_mu_con = subm(NonIgno_theta_con,range(0,dim-1),range(0,0) );
            NonIgno_Eta_con = subm(NonIgno_theta_con,range(dim,ParaMuSig-1),range(0,0) );
            EtaToL(NonIgno_Eta_con,dim,NonIgno_L_con);

            osumDts<<"constrained MLE of theta (mu,Eta,alpha): "<<trans(NonIgno_theta_con);

            osumDts<<"Assuming ignorable missingness, bd-constrained MLE of mean\n"\
            <<trans(NonIgno_mu_con)<<endl;
            NonIgno_Sig_con = NonIgno_L_con*trans(NonIgno_L_con);

            osumDts<<"Assuming ignorable missingness, bd-constrained MLE of Cov\n"\
            <<NonIgno_Sig_con<<endl;
           // osumDts<<"\n#(loglike evaluation during numerical optimization): "<<MLkOSamComCtr<<endl<<endl;

            osumDts<<"Seconds used in the computation of Tw: "<<Seconds_Tw_NonIgno<<"\n";
            osumDts<<"Seconds used in the computation of Tlr: "<<Seconds_Tlr_NonIgno<<"\n";
            osumDts<<"Tlr requires how much additional time (in \%): "\
            <<( (Seconds_Tlr_NonIgno - Seconds_Tw_NonIgno)/Seconds_Tw_NonIgno )*100<<"% \n";
   // Seconds_Tw_NonIgno_average  = Seconds_Tw_NonIgno_average + Seconds_Tw;
   // Seconds_Tlr_NonIgno_average = Seconds_Tlr_NonIgno_average + Seconds_Tlr;


            NonIgno_CovThetaHat_this = inv( HESScdif(mod_NoRi1_MNARd2_MLogLikGH, NonIgno_theta_con) );
            NonIgno_Avg_CovThetaHat_vh0 = NonIgno_Avg_CovThetaHat_vh0 + NonIgno_CovThetaHat_this;
            CovMuHat_this_vh0 = subm(NonIgno_CovThetaHat_this,range(0,dim-1),range(0,dim-1));
            InvCovMuHat_this_vh0 = inv(CovMuHat_this_vh0);
            osumDts<<"\nCov(MuHat) version h0: (inverse) observed information matrix evaluated at MLE under H0\n";
            osumDts<<CovMuHat_this_vh0<<"\n";

            osumDts<<"Cov(MuHat) version h0Uh1: (inverse) observed information matrix evaluated at MLE under H0UH1\n";
            osumDts<<CovMuHat_this_vh0Uh1<<"\n";

            osumDts<<"Cov(MuHat) reference: Sigma/SampleSize; it is exact when Muhat is the complete-sample mean \n";
            osumDts<<(Sig/Sz)<<"\n\n";

            osumDts<<"Test results based on Tw:\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            osumDts<<"Tw = "<<Tw_NonIgno<<"\n";
            pval_Tw_NonIgno_bd = 1-(scythe::pchisq(Tw_NonIgno, double(dim-1)) \
            + scythe::pchisq(Tw_NonIgno, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tw: "<<pval_Tw_NonIgno_bd<<"\n";

            // substitution p-value
            R_ThetaZeroInv_Rt_vh0 =  Sz*CovMuHat_this_vh0;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0,CbWts_vh0);
            pval_Tw_NonIgno_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tw_NonIgno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tw_NonIgno_sub_vh0<<"\n";

            R_ThetaZeroInv_Rt_vh0Uh1 =  Sz*CovMuHat_this_vh0Uh1;
            ChiBarWtsDim2ch2(R_ThetaZeroInv_Rt_vh0Uh1,CbWts_vh0Uh1);
            pval_Tw_NonIgno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tw_NonIgno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tw_NonIgno_sub_vh0Uh1<<"\n";

            osumDts<<"\n\n";
            osumDts<<"Test results based on Tlr:\n";
            osumDts<<"Tlr = "<<Tlr_NonIgno<<"\n";
            // upper bound p-value  is ( Pr(Xi_{dim-1}^2>=c)+ Pr(Xi_{dim}^2>=c) )/2
            // pchisq(x, double(dim-1) ) returns Pr(Xi_{dim-1}^2<=x)
            pval_Tlr_NonIgno_bd = 1-(scythe::pchisq(Tlr_NonIgno, double(dim-1)) \
            + scythe::pchisq(Tlr_NonIgno, double(dim)) )/2;
            osumDts<<"Upper bound of p-value using Tlr: "<<pval_Tlr_NonIgno_bd<<"\n";

            // substitution p-value
            pval_Tlr_NonIgno_sub_vh0 = ChiBarPvalue(dim,CbWts_vh0,Tlr_NonIgno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0: "\
            <<pval_Tlr_NonIgno_sub_vh0<<"\n";

            pval_Tlr_NonIgno_sub_vh0Uh1 = ChiBarPvalue(dim,CbWts_vh0Uh1,Tlr_NonIgno);
            osumDts<<"Substitution method p-value, using Cov(MuHat) version h0Uh1: "\
            <<pval_Tlr_NonIgno_sub_vh0Uh1<<"\n";
            }


            misratesXi1 = MissingRates(0);  misratesXi2 = MissingRates(1);
            var1_unc = NonIgno_Sig_unc(0,0);  var2_unc = NonIgno_Sig_unc(1,1);
            cov12_unc = NonIgno_Sig_unc(1,0);
            var1_con = NonIgno_Sig_con(0,0);  var2_con = NonIgno_Sig_con(1,1);
            cov12_con = NonIgno_Sig_con(1,0);

            mu1_unc = NonIgno_theta_unc(0); mu2_unc = NonIgno_theta_unc(1);
            alpha01_unc = NonIgno_theta_unc(5); alpha11_unc = NonIgno_theta_unc(6);
            alpha02_unc = NonIgno_theta_unc(7); alpha12_unc = NonIgno_theta_unc(8);

            mu1_con = NonIgno_theta_con(0); mu2_con = NonIgno_theta_con(1);
            alpha01_con = NonIgno_theta_con(5); alpha11_con = NonIgno_theta_con(6);
            alpha02_con = NonIgno_theta_con(7); alpha12_con = NonIgno_theta_con(8);

            // write log-file
        otLog<<setprecision (6)<< fixed\
        <<pval_Tlr_CC_bd<<","<<pval_Tw_CC_bd<<","<<pval_Tw_Igno_bd<<","<<pval_Tlr_Igno_bd<<","<<\
        pval_Tw_NonIgno_bd<<","<<pval_Tlr_NonIgno_bd<<","<<pval_Tlr_CC_sub_vh0Uh1<<","<<pval_Tw_CC_sub_vh0Uh1<<","<<\
        pval_Tw_Igno_sub_vh0Uh1<<","<<pval_Tlr_Igno_sub_vh0Uh1<<","<< pval_Tw_NonIgno_sub_vh0Uh1<<","<<\
        pval_Tlr_NonIgno_sub_vh0Uh1<<","<<pval_Tlr_CC_sub_vh0<<","<<pval_Tw_CC_sub_vh0<<","<< pval_Tw_Igno_sub_vh0<<","<<\
        pval_Tlr_Igno_sub_vh0<<","<< pval_Tw_NonIgno_sub_vh0<<","<< pval_Tlr_NonIgno_sub_vh0<<","<<\
        misratesXi1<<","<< misratesXi2<<","<< Tw_CC<<","<<Tlr_CC<<","<< Tw_Igno<<","<< Tlr_Igno<<","<< Tw_NonIgno<<","<< Tlr_NonIgno<<","<<\
        mu1_unc<<","<< mu2_unc<<","<< var1_unc<<","<< var2_unc<<","<< cov12_unc<<","<<\
        alpha01_unc<<","<< alpha11_unc<<","<< alpha02_unc<<","<< alpha12_unc<<","<<\
        mu1_con<<","<< mu2_con<<","<< var1_con<<","<< var2_con<<","<< cov12_con<<","<<\
        alpha01_con<<","<< alpha11_con<<","<< alpha02_con<<","<< alpha12_con<<","\
        <<Seconds_Tw_Igno<<","<<Seconds_Tlr_Igno<<","\
        <<Seconds_Tw_NonIgno<<","<<Seconds_Tlr_NonIgno<<"\n";
 // Seconds_Tw_Igno, Seconds_Tlr_Igno,Seconds_Tw_NonIgno,Seconds_Tlr_NonIgno

        SimuDone++; }

        osum<<"Number of simulation instances skipped due to non-convergence: "<<SimuSkipped<<endl;
        osum<<"Number of H_0 rejections among these skipped instances (based on full data  (missing+observed) Wald test: "\
        <<Rej_H0_counts_completeData_skippedInstances<<endl;
        osum<<"Proportion of H_0 rejections among these skipped instances (based on full data  (missing+observed) Wald test: "\
        <<double(Rej_H0_counts_completeData_skippedInstances)/double(SimuSkipped)<<endl;


        osum<<"NonIgno_Avg_CovThetaHat_vh0Uh1:\n"<<NonIgno_Avg_CovThetaHat_vh0Uh1<<"\n\n";
        osum<<"NonIgno_Avg_CovThetaHat_vh0:\n"<<NonIgno_Avg_CovThetaHat_vh0<<"\n\n";

        osum<<"Igno_Avg_CovThetaHat_vh0Uh1:\n"<<Igno_Avg_CovThetaHat_vh0Uh1<<"\n\n";
        osum<<"Igno_Avg_CovThetaHat_vh0:\n"<<Igno_Avg_CovThetaHat_vh0<<"\n\n";

        osum<<"CC_Avg_CovThetaHat_vh0Uh1:\n"<<CC_Avg_CovThetaHat_vh0Uh1<<"\n\n";
       // osum<<"CC_Avg_CovThetaHat_vh0:\n"<<CC_Avg_CovThetaHat_vh0<<"\n\n";

        osum<<"Reference value for power/type I error rate: number of H_0 rejections base on the complete data Wald/LRT test="\
        <<Rej_H0_counts_completeData<<endl;
        osum<<"Proportion of H_0 rejections base on the complete data Wald/LRT test="\
        <<double(Rej_H0_counts_completeData)/double(SimuDone)<<endl;

    }
    else
    {   cout<<"Goodbye\n";
    }

    return 0;

}
