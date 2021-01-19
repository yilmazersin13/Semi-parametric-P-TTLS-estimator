%SIMULATION STUDY FOR SEMIPARAMETRIC ESTIMATOR P-TTLS (RATIONAL
%APROXIMATION)
%--------------------------------------------------------------
%Simulation designed for 18 different combinations that based on two
%nonparametric functýons, three semaple sizes (n=50,150,300) and three
%variance functýons. 
%This simulation study needs m.file for smoothing spline function for
%semiparametric regression model and Rational approximation. 
%-----------------------------------------------------------------
clear
clc
sim=20;
for nf=1:2          %for loop for two nonparametric reg. functions.
    for vf=1:3      %for loop for three variance functions.
        for ss=1:3  %for loop for three sample sizes.
            if ss==1 
                n=50;
            end
            if ss==2
                n=150;
            end
            if ss==3
                n=300;
            end
         for sj=1:sim
%%
%DATA GENERATION-------------------------------------------------
%Parametric Component---------------------------------------------
beta=[-2.5;3];           %Regression coeffs.
x=rand(n,2);             %Parametric covariates from Uniform dist.
%Nonparametric component-----------------------------------------
z=zeros(n,1); z1=zeros(n/2,1); z2=zeros(n/2,1);
for i=1:n/2
    z1(i)=(i-0.5)/(n/2);      %Nonparametric covariate z_i
    z2(i)=(i-0.5)/(n/2);      %Nonparametric covariate z_i
end
for i=1:n
    z(i)=(i-0.5)/(n);         %Nonparametric covariate z_i
end
if nf==1                      
    f1=exp(z1)+z1;
    f2=sin(2*pi*z2);
    func1=[f1;f2];            %Nonparametric function I
    f=func1;
end
if nf==2
    f1=sin(2*pi*z1);
    f2=exp(z2)+sin(2*pi*z2);
    func2=[f1;f2];            %Nonparametric function II
    f=func2;
end

%Variance functions---------------------------------------------
vfnc=(0.15*(1+0.4*(2*vf-7)*(z-0.5)))*0.2;
% Generating model----------------------------------------------
f=rescale(f,-2,2);
e=randn(n,1);                 %Random error terms with N(0,1)
y=x*beta+f+sqrt(vfnc).*e;         %Semiparametric regression model
if ss==1
    y50=y;
 
end
if ss==2
    y150=y;
    
end
if ss==3
    y300=y;
    
end
%%
if ss==1 && vf==1 && nf==1
    f5011=f+sqrt(vfnc).*e;
end
if ss==2 && vf==1 && nf==1
    f15011=f+sqrt(vfnc).*e;
end
if ss==3 && vf==1 && nf==1
    f30011=f+sqrt(vfnc).*e;
end

if ss==1 && vf==1 && nf==2
    f5021=f+sqrt(vfnc).*e;
end
if ss==2 && vf==1 && nf==2
    f15021=f+sqrt(vfnc).*e;
end
if ss==3 && vf==1 && nf==2
    f30021=f+sqrt(vfnc).*e;
end

if ss==1 && vf==2 && nf==1
    f5012=f+sqrt(vfnc).*e;
end
if ss==2 && vf==2 && nf==1
    f15012=f+sqrt(vfnc).*e;
end
if ss==3 && vf==2 && nf==1
    f30012=f+sqrt(vfnc).*e;
end

if ss==1 && vf==2 && nf==2
    f5022=f+sqrt(vfnc).*e;
end
if ss==2 && vf==2 && nf==2
    f15022=f+sqrt(vfnc).*e;
end
if ss==3 && vf==2 && nf==2
    f30022=f+sqrt(vfnc).*e;
end

if ss==1 && vf==3 && nf==1
    f5013=f+sqrt(vfnc).*e;
end
if ss==2 && vf==3 && nf==1
    f15013=f+sqrt(vfnc).*e;
end
if ss==3 && vf==3 && nf==1
    f30013=f+sqrt(vfnc).*e;
end

if ss==1 && vf==3 && nf==2
    f5023=f+sqrt(vfnc).*e;
end
if ss==2 && vf==3 && nf==2
    f15023=f+sqrt(vfnc).*e;
end
if ss==3 && vf==3 && nf==2
    f30023=f+sqrt(vfnc).*e;
end
%%

if ss==1 && nf==1
    func1_50=rescale(func1,-2,2);
    yr=rescale(y,-2.99,3.6);
end
if ss==1 && nf==2
    func2_50=rescale(func2,-2,2);
    yr=rescale(y,-3.3,3.5);
end
if ss==2 && nf==1
    func1_150=rescale(func1,-2,2);
    yr=rescale(y,-3.2,3.85);
end
if ss==2 && nf==2
    func2_150=rescale(func2,-2,2);
    yr=rescale(y,-3.6,3.85);
end
if ss==3 && nf==1
    func1_300=rescale(func1,-2,2);
    yr=rescale(y,-3.25,3.9);
end
if ss==3 && nf==2
    func2_300=rescale(func2,-2,2);
    yr=rescale(y,-3.7,3.9);
end
%%
%%
%SMOOTHING SPLINE ESTIMATION
[ssbetahat,ssfhat,ssyhat,sbias,svarbeta,ssmde,sadjR,sGMSE]=spss(x,z,y);
S(ss).betahat(:,sj,vf)=ssbetahat;
S(ss).fhat(:,sj,vf)=ssfhat;
S(ss).yhat(:,sj,vf)=ssyhat;
S(ss).bias(:,sj,vf)=sbias;
S(ss).varbeta(:,sj,vf)=diag(svarbeta);
S(ss).smde(:,sj,vf)=ssmde;
S(ss).adjR(:,sj,vf)=sadjR;
S(ss).GMSE(:,sj,vf)=sGMSE;
%%
%P-TTLS ESTIMATION
[rbetahat,rfhat,ryhat,rankpttls,HatR,deg,rbias,rvarbeta,rsmde,radjR,rGMSE]=sppttls(x,z,yr,y);
R(ss).betahat(:,sj,vf)=rbetahat;
R(ss).fhat(:,sj,vf)=rfhat;
R(ss).yhat(:,sj,vf)=ryhat;
R(ss).bias(:,sj,vf)=rbias;
R(ss).varbeta(:,sj,vf)=diag(rvarbeta);
R(ss).smde(:,sj,vf)=rsmde;
R(ss).adjR(:,sj,vf)=radjR;
R(ss).GMSE(:,sj,vf)=rGMSE;

aa = sprintf('%d th simulation for sample size %d , var.func %d and reg. function %d is over.',sj,n,vf,nf);
disp(aa)
         end
        end
    end
if nf==1;
W1 = whos;
for ii = 1:length(W1)
    Model1.(W1(ii).name) = eval(W1(ii).name);
end
end
if nf==2;
W2 = whos;
for ii = 1:length(W2)
    Model2.(W2(ii).name) = eval(W2(ii).name);
end
end
end
%clearvars -except Model1 Model2 S R

%%
%RESULTS FOR SIMULATION STUDY---------------------------------------
%PARAMETRIC COMPONENT
%%
%Scenario 1: n=50, vf=1, function1
sc1.betaHatR=mean(Model1.R(1).betahat(:,:,1)')';
sc1.betaHatS=mean(Model1.S(1).betahat(:,:,1)')';
sc1.BiasR=abs(Model1.beta-sc1.betaHatR);
sc1.BiasS=abs(Model1.beta-sc1.betaHatS);
sc1.VarBetaR=mean(Model1.R(1).varbeta(:,:,1)')';
sc1.VarBetaS=mean(Model1.S(1).varbeta(:,:,1)')';
sc1.SmdeR=mean(Model1.R(1).smde(:,:,1));
sc1.SmdeS=mean(Model1.S(1).smde(:,:,1));
sc1.RER=sc1.SmdeR/sc1.SmdeS;
sc1.RES=sc1.SmdeS/sc1.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc1.BiasSQR=mean(abs(Model1.beta-sc1.betaHatR).^2);
sc1.BiasSQS=mean(abs(Model1.beta-sc1.betaHatS).^2);
sc1.var_smdeR=mean(sc1.VarBetaR);
sc1.var_smdeS=mean(sc1.VarBetaS);
sc1.plotSMDER=sc1.BiasSQR+sc1.var_smdeR;
sc1.plotSMDES=sc1.BiasSQS+sc1.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc1.fHatR=mean(Model1.R(1).fhat(:,:,1)')';
sc1.fHatS=mean(Model1.S(1).fhat(:,:,1)')';
sc1.MseR=mean((Model1.func1_50-sc1.fHatR).^2);
sc1.MseS=mean((Model1.func1_50-sc1.fHatS).^2);
sc1.RoMseR=sc1.MseR/sc1.MseS;
sc1.RoMseS=sc1.MseS/sc1.MseR;
%%
%MODEL PERFORMANCE
sc1.radjR2=mean(Model1.R(1).adjR(:,:,1));
sc1.rGMSE=mean(Model1.R(1).GMSE(:,:,1));
sc1.sadjR2=mean(Model1.S(1).adjR(:,:,1));
sc1.sGMSE=mean(Model1.S(1).GMSE(:,:,1));
%%
%Scenario 2: n=50, vf=2, function1
sc2.betaHatR=mean(Model1.R(1).betahat(:,:,2)')';
sc2.betaHatS=mean(Model1.S(1).betahat(:,:,2)')';
sc2.BiasR=abs(Model1.beta-sc2.betaHatR);
sc2.BiasS=abs(Model1.beta-sc2.betaHatS);
sc2.VarBetaR=mean(Model1.R(1).varbeta(:,:,2)')';
sc2.VarBetaS=mean(Model1.S(1).varbeta(:,:,2)')';
sc2.SmdeR=mean(Model1.R(1).smde(:,:,2));
sc2.SmdeS=mean(Model1.S(1).smde(:,:,2));
sc2.RER=sc2.SmdeR/sc2.SmdeS;
sc2.RES=sc2.SmdeS/sc2.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc2.BiasSQR=mean(abs(Model1.beta-sc2.betaHatR).^2);
sc2.BiasSQS=mean(abs(Model1.beta-sc2.betaHatS).^2);
sc2.var_smdeR=mean(sc2.VarBetaR);
sc2.var_smdeS=mean(sc2.VarBetaS);
sc2.plotSMDER=sc2.BiasSQR+sc2.var_smdeR;
sc2.plotSMDES=sc2.BiasSQS+sc2.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc2.fHatR=mean(Model1.R(1).fhat(:,:,2)')';
sc2.fHatS=mean(Model1.S(1).fhat(:,:,2)')';
sc2.MseR=mean((Model1.func1_50-sc2.fHatR).^2);
sc2.MseS=mean((Model1.func1_50-sc2.fHatS).^2);
sc2.RoMseR=sc2.MseR/sc2.MseS;
sc2.RoMseS=sc2.MseS/sc2.MseR;
%%
%MODEL PERFORMANCE
sc2.radjR2=mean(Model1.R(1).adjR(:,:,2));
sc2.rGMSE=mean(Model1.R(1).GMSE(:,:,2));
sc2.sadjR2=mean(Model1.S(1).adjR(:,:,2));
sc2.sGMSE=mean(Model1.S(1).GMSE(:,:,2));
%%
%Scenario 3: n=50, vf=3, function1
sc3.betaHatR=mean(Model1.R(1).betahat(:,:,3)')';
sc3.betaHatS=mean(Model1.S(1).betahat(:,:,3)')';
sc3.BiasR=abs(Model1.beta-sc3.betaHatR);
sc3.BiasS=abs(Model1.beta-sc3.betaHatS);
sc3.VarBetaR=mean(Model1.R(1).varbeta(:,:,3)')';
sc3.VarBetaS=mean(Model1.S(1).varbeta(:,:,3)')';
sc3.SmdeR=mean(Model1.R(1).smde(:,:,3));
sc3.SmdeS=mean(Model1.S(1).smde(:,:,3));
sc3.RER=sc3.SmdeR/sc3.SmdeS;
sc3.RES=sc3.SmdeS/sc3.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc3.BiasSQR=mean(abs(Model1.beta-sc3.betaHatR).^2);
sc3.BiasSQS=mean(abs(Model1.beta-sc3.betaHatS).^2);
sc3.var_smdeR=mean(sc3.VarBetaR);
sc3.var_smdeS=mean(sc3.VarBetaS);
sc3.plotSMDER=sc3.BiasSQR+sc3.var_smdeR;
sc3.plotSMDES=sc3.BiasSQS+sc3.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc3.fHatR=mean(Model1.R(1).fhat(:,:,3)')';
sc3.fHatS=mean(Model1.S(1).fhat(:,:,3)')';
sc3.MseR=mean((Model1.func1_50-sc3.fHatR).^3);
sc3.MseS=mean((Model1.func1_50-sc3.fHatS).^3);
sc3.RoMseR=sc3.MseR/sc3.MseS;
sc3.RoMseS=sc3.MseS/sc3.MseR;
%%
%MODEL PERFORMANCE
sc3.radjR2=mean(Model1.R(1).adjR(:,:,3));
sc3.rGMSE=mean(Model1.R(1).GMSE(:,:,3));
sc3.sadjR2=mean(Model1.S(1).adjR(:,:,3));
sc3.sGMSE=mean(Model1.S(1).GMSE(:,:,3));
%%
%---------------------------------------------------------
%Scenario 4: n=150, vf=1, function1
sc4.betaHatR=mean(Model1.R(2).betahat(:,:,1)')';
sc4.betaHatS=mean(Model1.S(2).betahat(:,:,1)')';
sc4.BiasR=abs(Model1.beta-sc4.betaHatR);
sc4.BiasS=abs(Model1.beta-sc4.betaHatS);
sc4.VarBetaR=mean(Model1.R(2).varbeta(:,:,1)')';
sc4.VarBetaS=mean(Model1.S(2).varbeta(:,:,1)')';
sc4.SmdeR=mean(Model1.R(2).smde(:,:,1));
sc4.SmdeS=mean(Model1.S(2).smde(:,:,1));
sc4.RER=sc4.SmdeR/sc4.SmdeS;
sc4.RES=sc4.SmdeS/sc4.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc4.BiasSQR=mean(abs(Model1.beta-sc4.betaHatR).^2);
sc4.BiasSQS=mean(abs(Model1.beta-sc4.betaHatS).^2);
sc4.var_smdeR=mean(sc4.VarBetaR);
sc4.var_smdeS=mean(sc4.VarBetaS);
sc4.plotSMDER=sc4.BiasSQR+sc4.var_smdeR;
sc4.plotSMDES=sc4.BiasSQS+sc4.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc4.fHatR=mean(Model1.R(2).fhat(:,:,1)')';
sc4.fHatS=mean(Model1.S(2).fhat(:,:,1)')';
sc4.MseR=mean((Model1.func1_150-sc4.fHatR).^2);
sc4.MseS=mean((Model1.func1_150-sc4.fHatS).^2);
sc4.RoMseR=sc4.MseR/sc4.MseS;
sc4.RoMseS=sc4.MseS/sc4.MseR;
%%
%MODEL PERFORMANCE
sc4.radjR2=mean(Model1.R(2).adjR(:,:,1));
sc4.rGMSE=mean(Model1.R(2).GMSE(:,:,1));
sc4.sadjR2=mean(Model1.S(2).adjR(:,:,1));
sc4.sGMSE=mean(Model1.S(2).GMSE(:,:,1));
%%
%Scenario 5: n=150, vf=2, function1
sc5.betaHatR=mean(Model1.R(2).betahat(:,:,2)')';
sc5.betaHatS=mean(Model1.S(2).betahat(:,:,2)')';
sc5.BiasR=abs(Model1.beta-sc5.betaHatR);
sc5.BiasS=abs(Model1.beta-sc5.betaHatS);
sc5.VarBetaR=mean(Model1.R(2).varbeta(:,:,2)')';
sc5.VarBetaS=mean(Model1.S(2).varbeta(:,:,2)')';
sc5.SmdeR=mean(Model1.R(2).smde(:,:,2));
sc5.SmdeS=mean(Model1.S(2).smde(:,:,2));
sc5.RER=sc5.SmdeR/sc5.SmdeS;
sc5.RES=sc5.SmdeS/sc5.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc5.BiasSQR=mean(abs(Model1.beta-sc5.betaHatR).^2);
sc5.BiasSQS=mean(abs(Model1.beta-sc5.betaHatS).^2);
sc5.var_smdeR=mean(sc5.VarBetaR);
sc5.var_smdeS=mean(sc5.VarBetaS);
sc5.plotSMDER=sc5.BiasSQR+sc5.var_smdeR;
sc5.plotSMDES=sc5.BiasSQS+sc5.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc5.fHatR=mean(Model1.R(2).fhat(:,:,2)')';
sc5.fHatS=mean(Model1.S(2).fhat(:,:,2)')';
sc5.MseR=mean((Model1.func1_150-sc5.fHatR).^2);
sc5.MseS=mean((Model1.func1_150-sc5.fHatS).^2);
sc5.RoMseR=sc5.MseR/sc5.MseS;
sc5.RoMseS=sc5.MseS/sc5.MseR;
%%
%MODEL PERFORMANCE
sc5.radjR2=mean(Model1.R(2).adjR(:,:,2));
sc5.rGMSE=mean(Model1.R(2).GMSE(:,:,2));
sc5.sadjR2=mean(Model1.S(2).adjR(:,:,2));
sc5.sGMSE=mean(Model1.S(2).GMSE(:,:,2));
%%
%Scenario 6: n=150, vf=3, function1
sc6.betaHatR=mean(Model1.R(2).betahat(:,:,3)')';
sc6.betaHatS=mean(Model1.S(2).betahat(:,:,3)')';
sc6.BiasR=abs(Model1.beta-sc6.betaHatR);
sc6.BiasS=abs(Model1.beta-sc6.betaHatS);
sc6.VarBetaR=mean(Model1.R(2).varbeta(:,:,3)')';
sc6.VarBetaS=mean(Model1.S(2).varbeta(:,:,3)')';
sc6.SmdeR=mean(Model1.R(2).smde(:,:,3));
sc6.SmdeS=mean(Model1.S(2).smde(:,:,3));
sc6.RER=sc6.SmdeR/sc6.SmdeS;
sc6.RES=sc6.SmdeS/sc6.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc6.BiasSQR=mean(abs(Model1.beta-sc6.betaHatR).^2);
sc6.BiasSQS=mean(abs(Model1.beta-sc6.betaHatS).^2);
sc6.var_smdeR=mean(sc6.VarBetaR);
sc6.var_smdeS=mean(sc6.VarBetaS);
sc6.plotSMDER=sc6.BiasSQR+sc6.var_smdeR;
sc6.plotSMDES=sc6.BiasSQS+sc6.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc6.fHatR=mean(Model1.R(2).fhat(:,:,3)')';
sc6.fHatS=mean(Model1.S(2).fhat(:,:,3)')';
sc6.MseR=mean((Model1.func1_150-sc6.fHatR).^3);
sc6.MseS=mean((Model1.func1_150-sc6.fHatS).^3);
sc6.RoMseR=sc6.MseR/sc6.MseS;
sc6.RoMseS=sc6.MseS/sc6.MseR;
%%
%MODEL PERFORMANCE
sc6.radjR2=mean(Model1.R(2).adjR(:,:,3));
sc6.rGMSE=mean(Model1.R(2).GMSE(:,:,3));
sc6.sadjR2=mean(Model1.S(2).adjR(:,:,3));
sc6.sGMSE=mean(Model1.S(2).GMSE(:,:,3));
%%
%-------------------------------------------------------------
%Scenario 7: n=300, vf=1, function1
sc7.betaHatR=mean(Model1.R(3).betahat(:,:,1)')';
sc7.betaHatS=mean(Model1.S(3).betahat(:,:,1)')';
sc7.BiasR=abs(Model1.beta-sc7.betaHatR);
sc7.BiasS=abs(Model1.beta-sc7.betaHatS);
sc7.VarBetaR=mean(Model1.R(3).varbeta(:,:,1)')';
sc7.VarBetaS=mean(Model1.S(3).varbeta(:,:,1)')';
sc7.SmdeR=mean(Model1.R(3).smde(:,:,1));
sc7.SmdeS=mean(Model1.S(3).smde(:,:,1));
sc7.RER=sc7.SmdeR/sc7.SmdeS;
sc7.RES=sc7.SmdeS/sc7.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc7.BiasSQR=mean(abs(Model1.beta-sc7.betaHatR).^2);
sc7.BiasSQS=mean(abs(Model1.beta-sc7.betaHatS).^2);
sc7.var_smdeR=mean(sc7.VarBetaR);
sc7.var_smdeS=mean(sc7.VarBetaS);
sc7.plotSMDER=sc7.BiasSQR+sc7.var_smdeR;
sc7.plotSMDES=sc7.BiasSQS+sc7.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc7.fHatR=mean(Model1.R(3).fhat(:,:,1)')';
sc7.fHatS=mean(Model1.S(3).fhat(:,:,1)')';
sc7.MseR=mean((Model1.func1_300-sc7.fHatR).^2);
sc7.MseS=mean((Model1.func1_300-sc7.fHatS).^2);
sc7.RoMseR=sc7.MseR/sc7.MseS;
sc7.RoMseS=sc7.MseS/sc7.MseR;
%%
%MODEL PERFORMANCE
sc7.radjR2=mean(Model1.R(3).adjR(:,:,1));
sc7.rGMSE=mean(Model1.R(3).GMSE(:,:,1));
sc7.sadjR2=mean(Model1.S(3).adjR(:,:,1));
sc7.sGMSE=mean(Model1.S(3).GMSE(:,:,1));
%%
%Scenario 8: n=300, vf=2, function1
sc8.betaHatR=mean(Model1.R(3).betahat(:,:,2)')';
sc8.betaHatS=mean(Model1.S(3).betahat(:,:,2)')';
sc8.BiasR=abs(Model1.beta-sc8.betaHatR);
sc8.BiasS=abs(Model1.beta-sc8.betaHatS);
sc8.VarBetaR=mean(Model1.R(3).varbeta(:,:,2)')';
sc8.VarBetaS=mean(Model1.S(3).varbeta(:,:,2)')';
sc8.SmdeR=mean(Model1.R(3).smde(:,:,2));
sc8.SmdeS=mean(Model1.S(3).smde(:,:,2));
sc8.RER=sc8.SmdeR/sc8.SmdeS;
sc8.RES=sc8.SmdeS/sc8.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc8.BiasSQR=mean(abs(Model1.beta-sc8.betaHatR).^2);
sc8.BiasSQS=mean(abs(Model1.beta-sc8.betaHatS).^2);
sc8.var_smdeR=mean(sc8.VarBetaR);
sc8.var_smdeS=mean(sc8.VarBetaS);
sc8.plotSMDER=sc8.BiasSQR+sc8.var_smdeR;
sc8.plotSMDES=sc8.BiasSQS+sc8.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc8.fHatR=mean(Model1.R(3).fhat(:,:,2)')';
sc8.fHatS=mean(Model1.S(3).fhat(:,:,2)')';
sc8.MseR=mean((Model1.func1_300-sc8.fHatR).^2);
sc8.MseS=mean((Model1.func1_300-sc8.fHatS).^2);
sc8.RoMseR=sc8.MseR/sc8.MseS;
sc8.RoMseS=sc8.MseS/sc8.MseR;
%%
%MODEL PERFORMANCE
sc8.radjR2=mean(Model1.R(3).adjR(:,:,2));
sc8.rGMSE=mean(Model1.R(3).GMSE(:,:,2));
sc8.sadjR2=mean(Model1.S(3).adjR(:,:,2));
sc8.sGMSE=mean(Model1.S(3).GMSE(:,:,2));
%%
%Scenario 9: n=300, vf=3, function1
sc9.betaHatR=mean(Model1.R(3).betahat(:,:,3)')';
sc9.betaHatS=mean(Model1.S(3).betahat(:,:,3)')';
sc9.BiasR=abs(Model1.beta-sc9.betaHatR);
sc9.BiasS=abs(Model1.beta-sc9.betaHatS);
sc9.VarBetaR=mean(Model1.R(3).varbeta(:,:,3)')';
sc9.VarBetaS=mean(Model1.S(3).varbeta(:,:,3)')';
sc9.SmdeR=mean(Model1.R(3).smde(:,:,3));
sc9.SmdeS=mean(Model1.S(3).smde(:,:,3));
sc9.RER=sc9.SmdeR/sc9.SmdeS;
sc9.RES=sc9.SmdeS/sc9.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc9.BiasSQR=mean(abs(Model1.beta-sc9.betaHatR).^2);
sc9.BiasSQS=mean(abs(Model1.beta-sc9.betaHatS).^2);
sc9.var_smdeR=mean(sc9.VarBetaR);
sc9.var_smdeS=mean(sc9.VarBetaS);
sc9.plotSMDER=sc9.BiasSQR+sc9.var_smdeR;
sc9.plotSMDES=sc9.BiasSQS+sc9.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc9.fHatR=mean(Model1.R(3).fhat(:,:,3)')';
sc9.fHatS=mean(Model1.S(3).fhat(:,:,3)')';
sc9.MseR=mean((Model1.func1_300-sc9.fHatR).^3);
sc9.MseS=mean((Model1.func1_300-sc9.fHatS).^3);
sc9.RoMseR=sc9.MseR/sc9.MseS;
sc9.RoMseS=sc9.MseS/sc9.MseR;
%%
%MODEL PERFORMANCE
sc9.radjR2=mean(Model1.R(3).adjR(:,:,3));
sc9.rGMSE=mean(Model1.R(3).GMSE(:,:,3));
sc9.sadjR2=mean(Model1.S(3).adjR(:,:,3));
sc9.sGMSE=mean(Model1.S(3).GMSE(:,:,3));
%%
%-------------------------------------------------------------------
%MODEL2--------------------------------------------------------------
%%
%Scenario 10: n=50, vf=1, function2
sc10.betaHatR=mean(Model2.R(1).betahat(:,:,1)')';
sc10.betaHatS=mean(Model2.S(1).betahat(:,:,1)')';
sc10.BiasR=abs(Model2.beta-sc10.betaHatR);
sc10.BiasS=abs(Model2.beta-sc10.betaHatS);
sc10.VarBetaR=mean(Model2.R(1).varbeta(:,:,1)')';
sc10.VarBetaS=mean(Model2.S(1).varbeta(:,:,1)')';
sc10.SmdeR=mean(Model2.R(1).smde(:,:,1));
sc10.SmdeS=mean(Model2.S(1).smde(:,:,1));
sc10.RER=sc10.SmdeR/sc10.SmdeS;
sc10.RES=sc10.SmdeS/sc10.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc10.BiasSQR=mean(abs(Model2.beta-sc10.betaHatR).^2);
sc10.BiasSQS=mean(abs(Model2.beta-sc10.betaHatS).^2);
sc10.var_smdeR=mean(sc10.VarBetaR);
sc10.var_smdeS=mean(sc10.VarBetaS);
sc10.plotSMDER=sc10.BiasSQR+sc10.var_smdeR;
sc10.plotSMDES=sc10.BiasSQS+sc10.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc10.fHatR=mean(Model2.R(1).fhat(:,:,1)')';
sc10.fHatS=mean(Model2.S(1).fhat(:,:,1)')';
sc10.MseR=mean((Model2.func1_50-sc10.fHatR).^2);
sc10.MseS=mean((Model2.func1_50-sc10.fHatS).^2);
sc10.RoMseR=sc10.MseR/sc10.MseS;
sc10.RoMseS=sc10.MseS/sc10.MseR;
%%
%MODEL PERFORMANCE
sc10.radjR2=mean(Model2.R(1).adjR(:,:,1));
sc10.rGMSE=mean(Model2.R(1).GMSE(:,:,1));
sc10.sadjR2=mean(Model2.S(1).adjR(:,:,1));
sc10.sGMSE=mean(Model2.S(1).GMSE(:,:,1));
%%
%Scenario 11: n=50, vf=2, function2
sc11.betaHatR=mean(Model2.R(1).betahat(:,:,2)')';
sc11.betaHatS=mean(Model2.S(1).betahat(:,:,2)')';
sc11.BiasR=abs(Model2.beta-sc11.betaHatR);
sc11.BiasS=abs(Model2.beta-sc11.betaHatS);
sc11.VarBetaR=mean(Model2.R(1).varbeta(:,:,2)')';
sc11.VarBetaS=mean(Model2.S(1).varbeta(:,:,2)')';
sc11.SmdeR=mean(Model2.R(1).smde(:,:,2));
sc11.SmdeS=mean(Model2.S(1).smde(:,:,2));
sc11.RER=sc11.SmdeR/sc11.SmdeS;
sc11.RES=sc11.SmdeS/sc11.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc11.BiasSQR=mean(abs(Model2.beta-sc11.betaHatR).^2);
sc11.BiasSQS=mean(abs(Model2.beta-sc11.betaHatS).^2);
sc11.var_smdeR=mean(sc11.VarBetaR);
sc11.var_smdeS=mean(sc11.VarBetaS);
sc11.plotSMDER=sc11.BiasSQR+sc11.var_smdeR;
sc11.plotSMDES=sc11.BiasSQS+sc11.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc11.fHatR=mean(Model2.R(1).fhat(:,:,2)')';
sc11.fHatS=mean(Model2.S(1).fhat(:,:,2)')';
sc11.MseR=mean((Model2.func1_50-sc11.fHatR).^2);
sc11.MseS=mean((Model2.func1_50-sc11.fHatS).^2);
sc11.RoMseR=sc11.MseR/sc11.MseS;
sc11.RoMseS=sc11.MseS/sc11.MseR;
%%
%MODEL PERFORMANCE
sc11.radjR2=mean(Model2.R(1).adjR(:,:,2));
sc11.rGMSE=mean(Model2.R(1).GMSE(:,:,2));
sc11.sadjR2=mean(Model2.S(1).adjR(:,:,2));
sc11.sGMSE=mean(Model2.S(1).GMSE(:,:,2));
%%
%Scenario 12: n=50, vf=3, function2
sc12.betaHatR=mean(Model2.R(1).betahat(:,:,3)')';
sc12.betaHatS=mean(Model2.S(1).betahat(:,:,3)')';
sc12.BiasR=abs(Model2.beta-sc12.betaHatR);
sc12.BiasS=abs(Model2.beta-sc12.betaHatS);
sc12.VarBetaR=mean(Model2.R(1).varbeta(:,:,3)')';
sc12.VarBetaS=mean(Model2.S(1).varbeta(:,:,3)')';
sc12.SmdeR=mean(Model2.R(1).smde(:,:,3));
sc12.SmdeS=mean(Model2.S(1).smde(:,:,3));
sc12.RER=sc12.SmdeR/sc12.SmdeS;
sc12.RES=sc12.SmdeS/sc12.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc12.BiasSQR=mean(abs(Model2.beta-sc12.betaHatR).^2);
sc12.BiasSQS=mean(abs(Model2.beta-sc12.betaHatS).^2);
sc12.var_smdeR=mean(sc12.VarBetaR);
sc12.var_smdeS=mean(sc12.VarBetaS);
sc12.plotSMDER=sc12.BiasSQR+sc12.var_smdeR;
sc12.plotSMDES=sc12.BiasSQS+sc12.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc12.fHatR=mean(Model2.R(1).fhat(:,:,3)')';
sc12.fHatS=mean(Model2.S(1).fhat(:,:,3)')';
sc12.MseR=mean((Model2.func1_50-sc12.fHatR).^3);
sc12.MseS=mean((Model2.func1_50-sc12.fHatS).^3);
sc12.RoMseR=sc12.MseR/sc12.MseS;
sc12.RoMseS=sc12.MseS/sc12.MseR;
%%
%MODEL PERFORMANCE
sc12.radjR2=mean(Model2.R(1).adjR(:,:,3));
sc12.rGMSE=mean(Model2.R(1).GMSE(:,:,3));
sc12.sadjR2=mean(Model2.S(1).adjR(:,:,3));
sc12.sGMSE=mean(Model2.S(1).GMSE(:,:,3));
%%
%---------------------------------------------------------
%Scenario 13: n=150, vf=1, function2
sc13.betaHatR=mean(Model2.R(2).betahat(:,:,1)')';
sc13.betaHatS=mean(Model2.S(2).betahat(:,:,1)')';
sc13.BiasR=abs(Model2.beta-sc13.betaHatR);
sc13.BiasS=abs(Model2.beta-sc13.betaHatS);
sc13.VarBetaR=mean(Model2.R(2).varbeta(:,:,1)')';
sc13.VarBetaS=mean(Model2.S(2).varbeta(:,:,1)')';
sc13.SmdeR=mean(Model2.R(2).smde(:,:,1));
sc13.SmdeS=mean(Model2.S(2).smde(:,:,1));
sc13.RER=sc13.SmdeR/sc13.SmdeS;
sc13.RES=sc13.SmdeS/sc13.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc13.BiasSQR=mean(abs(Model2.beta-sc13.betaHatR).^2);
sc13.BiasSQS=mean(abs(Model2.beta-sc13.betaHatS).^2);
sc13.var_smdeR=mean(sc13.VarBetaR);
sc13.var_smdeS=mean(sc13.VarBetaS);
sc13.plotSMDER=sc13.BiasSQR+sc13.var_smdeR;
sc13.plotSMDES=sc13.BiasSQS+sc13.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc13.fHatR=mean(Model2.R(2).fhat(:,:,1)')';
sc13.fHatS=mean(Model2.S(2).fhat(:,:,1)')';
sc13.MseR=mean((Model2.func1_150-sc13.fHatR).^2);
sc13.MseS=mean((Model2.func1_150-sc13.fHatS).^2);
sc13.RoMseR=sc13.MseR/sc13.MseS;
sc13.RoMseS=sc13.MseS/sc13.MseR;
%%
%MODEL PERFORMANCE
sc13.radjR2=mean(Model2.R(2).adjR(:,:,1));
sc13.rGMSE=mean(Model2.R(2).GMSE(:,:,1));
sc13.sadjR2=mean(Model2.S(2).adjR(:,:,1));
sc13.sGMSE=mean(Model2.S(2).GMSE(:,:,1));
%%
%Scenario 14: n=150, vf=2, function2
sc14.betaHatR=mean(Model2.R(2).betahat(:,:,2)')';
sc14.betaHatS=mean(Model2.S(2).betahat(:,:,2)')';
sc14.BiasR=abs(Model2.beta-sc14.betaHatR);
sc14.BiasS=abs(Model2.beta-sc14.betaHatS);
sc14.VarBetaR=mean(Model2.R(2).varbeta(:,:,2)')';
sc14.VarBetaS=mean(Model2.S(2).varbeta(:,:,2)')';
sc14.SmdeR=mean(Model2.R(2).smde(:,:,2));
sc14.SmdeS=mean(Model2.S(2).smde(:,:,2));
sc14.RER=sc14.SmdeR/sc14.SmdeS;
sc14.RES=sc14.SmdeS/sc14.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc14.BiasSQR=mean(abs(Model2.beta-sc14.betaHatR).^2);
sc14.BiasSQS=mean(abs(Model2.beta-sc14.betaHatS).^2);
sc14.var_smdeR=mean(sc14.VarBetaR);
sc14.var_smdeS=mean(sc14.VarBetaS);
sc14.plotSMDER=sc14.BiasSQR+sc14.var_smdeR;
sc14.plotSMDES=sc14.BiasSQS+sc14.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc14.fHatR=mean(Model2.R(2).fhat(:,:,2)')';
sc14.fHatS=mean(Model2.S(2).fhat(:,:,2)')';
sc14.MseR=mean((Model2.func1_150-sc14.fHatR).^2);
sc14.MseS=mean((Model2.func1_150-sc14.fHatS).^2);
sc14.RoMseR=sc14.MseR/sc14.MseS;
sc14.RoMseS=sc14.MseS/sc14.MseR;
%%
%MODEL PERFORMANCE
sc14.radjR2=mean(Model2.R(2).adjR(:,:,2));
sc14.rGMSE=mean(Model2.R(2).GMSE(:,:,2));
sc14.sadjR2=mean(Model2.S(2).adjR(:,:,2));
sc14.sGMSE=mean(Model2.S(2).GMSE(:,:,2));
%%
%Scenario 15: n=150, vf=3, function2
sc15.betaHatR=mean(Model2.R(2).betahat(:,:,3)')';
sc15.betaHatS=mean(Model2.S(2).betahat(:,:,3)')';
sc15.BiasR=abs(Model2.beta-sc15.betaHatR);
sc15.BiasS=abs(Model2.beta-sc15.betaHatS);
sc15.VarBetaR=mean(Model2.R(2).varbeta(:,:,3)')';
sc15.VarBetaS=mean(Model2.S(2).varbeta(:,:,3)')';
sc15.SmdeR=mean(Model2.R(2).smde(:,:,3));
sc15.SmdeS=mean(Model2.S(2).smde(:,:,3));
sc15.RER=sc15.SmdeR/sc15.SmdeS;
sc15.RES=sc15.SmdeS/sc15.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc15.BiasSQR=mean(abs(Model2.beta-sc15.betaHatR).^2);
sc15.BiasSQS=mean(abs(Model2.beta-sc15.betaHatS).^2);
sc15.var_smdeR=mean(sc15.VarBetaR);
sc15.var_smdeS=mean(sc15.VarBetaS);
sc15.plotSMDER=sc15.BiasSQR+sc15.var_smdeR;
sc15.plotSMDES=sc15.BiasSQS+sc15.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc15.fHatR=mean(Model2.R(2).fhat(:,:,3)')';
sc15.fHatS=mean(Model2.S(2).fhat(:,:,3)')';
sc15.MseR=mean((Model2.func1_150-sc15.fHatR).^3);
sc15.MseS=mean((Model2.func1_150-sc15.fHatS).^3);
sc15.RoMseR=sc15.MseR/sc15.MseS;
sc15.RoMseS=sc15.MseS/sc15.MseR;
%%
%MODEL PERFORMANCE
sc15.radjR2=mean(Model2.R(2).adjR(:,:,3));
sc15.rGMSE=mean(Model2.R(2).GMSE(:,:,3));
sc15.sadjR2=mean(Model2.S(2).adjR(:,:,3));
sc15.sGMSE=mean(Model2.S(2).GMSE(:,:,3));
%%
%-------------------------------------------------------------
%Scenario 16: n=300, vf=1, function2
sc16.betaHatR=mean(Model2.R(3).betahat(:,:,1)')';
sc16.betaHatS=mean(Model2.S(3).betahat(:,:,1)')';
sc16.BiasR=abs(Model2.beta-sc16.betaHatR);
sc16.BiasS=abs(Model2.beta-sc16.betaHatS);
sc16.VarBetaR=mean(Model2.R(3).varbeta(:,:,1)')';
sc16.VarBetaS=mean(Model2.S(3).varbeta(:,:,1)')';
sc16.SmdeR=mean(Model2.R(3).smde(:,:,1));
sc16.SmdeS=mean(Model2.S(3).smde(:,:,1));
sc16.RER=sc16.SmdeR/sc16.SmdeS;
sc16.RES=sc16.SmdeS/sc16.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc16.BiasSQR=mean(abs(Model2.beta-sc16.betaHatR).^2);
sc16.BiasSQS=mean(abs(Model2.beta-sc16.betaHatS).^2);
sc16.var_smdeR=mean(sc16.VarBetaR);
sc16.var_smdeS=mean(sc16.VarBetaS);
sc16.plotSMDER=sc16.BiasSQR+sc16.var_smdeR;
sc16.plotSMDES=sc16.BiasSQS+sc16.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc16.fHatR=mean(Model2.R(3).fhat(:,:,1)')';
sc16.fHatS=mean(Model2.S(3).fhat(:,:,1)')';
sc16.MseR=mean((Model2.func1_300-sc16.fHatR).^2);
sc16.MseS=mean((Model2.func1_300-sc16.fHatS).^2);
sc16.RoMseR=sc16.MseR/sc16.MseS;
sc16.RoMseS=sc16.MseS/sc16.MseR;
%%
%MODEL PERFORMANCE
sc16.radjR2=mean(Model2.R(3).adjR(:,:,1));
sc16.rGMSE=mean(Model2.R(3).GMSE(:,:,1));
sc16.sadjR2=mean(Model2.S(3).adjR(:,:,1));
sc16.sGMSE=mean(Model2.S(3).GMSE(:,:,1));
%%
%Scenario 17: n=300, vf=2, function2
sc17.betaHatR=mean(Model2.R(3).betahat(:,:,2)')';
sc17.betaHatS=mean(Model2.S(3).betahat(:,:,2)')';
sc17.BiasR=abs(Model2.beta-sc17.betaHatR);
sc17.BiasS=abs(Model2.beta-sc17.betaHatS);
sc17.VarBetaR=mean(Model2.R(3).varbeta(:,:,2)')';
sc17.VarBetaS=mean(Model2.S(3).varbeta(:,:,2)')';
sc17.SmdeR=mean(Model2.R(3).smde(:,:,2));
sc17.SmdeS=mean(Model2.S(3).smde(:,:,2));
sc17.RER=sc17.SmdeR/sc17.SmdeS;
sc17.RES=sc17.SmdeS/sc17.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc17.BiasSQR=mean(abs(Model2.beta-sc17.betaHatR).^2);
sc17.BiasSQS=mean(abs(Model2.beta-sc17.betaHatS).^2);
sc17.var_smdeR=mean(sc17.VarBetaR);
sc17.var_smdeS=mean(sc17.VarBetaS);
sc17.plotSMDER=sc17.BiasSQR+sc17.var_smdeR;
sc17.plotSMDES=sc17.BiasSQS+sc17.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc17.fHatR=mean(Model2.R(3).fhat(:,:,2)')';
sc17.fHatS=mean(Model2.S(3).fhat(:,:,2)')';
sc17.MseR=mean((Model2.func1_300-sc17.fHatR).^2);
sc17.MseS=mean((Model2.func1_300-sc17.fHatS).^2);
sc17.RoMseR=sc17.MseR/sc17.MseS;
sc17.RoMseS=sc17.MseS/sc17.MseR;
%%
%MODEL PERFORMANCE
sc17.radjR2=mean(Model2.R(3).adjR(:,:,2));
sc17.rGMSE=mean(Model2.R(3).GMSE(:,:,2));
sc17.sadjR2=mean(Model2.S(3).adjR(:,:,2));
sc17.sGMSE=mean(Model2.S(3).GMSE(:,:,2));
%%
%Scenario 18: n=300, vf=3, function2
sc18.betaHatR=mean(Model2.R(3).betahat(:,:,3)')';
sc18.betaHatS=mean(Model2.S(3).betahat(:,:,3)')';
sc18.BiasR=abs(Model2.beta-sc18.betaHatR);
sc18.BiasS=abs(Model2.beta-sc18.betaHatS);
sc18.VarBetaR=mean(Model2.R(3).varbeta(:,:,3)')';
sc18.VarBetaS=mean(Model2.S(3).varbeta(:,:,3)')';
sc18.SmdeR=mean(Model2.R(3).smde(:,:,3));
sc18.SmdeS=mean(Model2.S(3).smde(:,:,3));
sc18.RER=sc18.SmdeR/sc18.SmdeS;
sc18.RES=sc18.SmdeS/sc18.SmdeR;
%FOR DECOMPOSITION VAR-BIAS
sc18.BiasSQR=mean(abs(Model2.beta-sc18.betaHatR).^2);
sc18.BiasSQS=mean(abs(Model2.beta-sc18.betaHatS).^2);
sc18.var_smdeR=mean(sc18.VarBetaR);
sc18.var_smdeS=mean(sc18.VarBetaS);
sc18.plotSMDER=sc18.BiasSQR+sc18.var_smdeR;
sc18.plotSMDES=sc18.BiasSQS+sc18.var_smdeS;
%%
%NONPARAMETRIC COMPONENT
sc18.fHatR=mean(Model2.R(3).fhat(:,:,3)')';
sc18.fHatS=mean(Model2.S(3).fhat(:,:,3)')';
sc18.MseR=mean((Model2.func1_300-sc18.fHatR).^3);
sc18.MseS=mean((Model2.func1_300-sc18.fHatS).^3);
sc18.RoMseR=sc18.MseR/sc18.MseS;
sc18.RoMseS=sc18.MseS/sc18.MseR;
%%
%MODEL PERFORMANCE
sc18.radjR2=mean(Model2.R(3).adjR(:,:,3));
sc18.rGMSE=mean(Model2.R(3).GMSE(:,:,3));
sc18.sadjR2=mean(Model2.S(3).adjR(:,:,3));
sc18.sGMSE=mean(Model2.S(3).GMSE(:,:,3));
%%
%---------------------------------------------------------------------

