

% This code is the Bayesian joint decision-estimation algorithm for
% determination, separation, and tracking of an unknown time-varying number of maneuvering sources 

clear all
close all
echo off
clc

%-----------------------Parameters

c=[0 1;1 0];

alpha=[1 1;1 1];
beta=[0.8 0.7;0.8 0.7];


Est_d=8;
%-----------------------

%=======================Trajectory

N=32;
T=0.5;

% Pd = 1;

%________________

Ns=8;

xs(1,1)=0; ys(1,1)=0;
xs(2,1)=100; ys(2,1)=0;
xs(3,1)=100; ys(3,1)=100;
xs(4,1)=0; ys(4,1)=100;
xs(5,1)=0; ys(5,1)=50;
xs(6,1)=50; ys(6,1)=0;
xs(7,1)=100; ys(7,1)=50;
xs(8,1)=50; ys(8,1)=100;

% %----------------------------Extra sensors
xs(9,1)=25; ys(9,1)=0;
xs(10,1)=75; ys(10,1)=0;
xs(11,1)=25; ys(11,1)=100;
xs(12,1)=75; ys(12,1)=100;
xs(13,1)=0; ys(13,1)=25;
xs(14,1)=0; ys(14,1)=75;
xs(15,1)=100; ys(15,1)=25;
xs(16,1)=100; ys(16,1)=75;

% xs(9,1)=80; ys(9,1)=100;
%--------------------------------

%________________TPM_H

TPM_H=[0.99 0.01;
    0.01 0.99];
%________________

P_X_NEW = [10 0 0 0;
           0 1 0 0;
           0 0 10 0;
           0 0 0 1];
       
       sigmaInipNEW = sqrt(10);
       sigmaInivNEW = sqrt(1);

%_______________ Model base

F1=[1 T 0 0;
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];

F2=blkdiag(F1,F1);

Qd =(0.003)^2;

G1=[T^2/2;T;T^2/2;T];

% Q = Qd;

G2=[T^2/2;T;T^2/2;T;T^2/2;T;T^2/2;T];

%_______________

Nmd=3;

Nm = Nmd^4;
Fv = zeros(8,Nm);

i_Nm = 0;

for i4=1:Nmd

    ay2=(i4-2);

    for i3=1:Nmd

        ax2=(i3-2);

        for i2=1:Nmd

            ay1=(i2-2);

            for i1=1:Nmd

                ax1=(i1-2);

                i_Nm = i_Nm + 1;

                %                 ax1=0; ay1=0; ax2=0; ay2=0;

                Fv(:,i_Nm)= [ax1*T^2/2;ax1*T;ay1*T^2/2;ay1*T;...
                    ax2*T^2/2;ax2*T;ay2*T^2/2;ay2*T];

                %-------------------------

                Mm(:,i_Nm)=[ax1;ay1;ax2;ay2];

                %-------------------------

            end
        end
    end
end

%-----------------------

dgl=0.9;
odgl=0.05;

for i=1:Nm

    for j=1:Nm

        p=1;
        for ic=1:4

            if abs(Mm(ic,i)-Mm(ic,j))==2

                p=0;

            end

        end

        if p~=0

            p=1;
            for ic=1:4

                if abs(Mm(ic,i)-Mm(ic,j))==0

                    p=p*dgl;

                elseif abs(Mm(ic,i)-Mm(ic,j))==1

                    p=p*odgl;

                end

            end

        end

        TPM_m_un(i,j)=p;

    end

end

%----------------------TPM normalization

Rtpm = sum(TPM_m_un,2);
MRtpm = diag(Rtpm);

TPM_m = inv(MRtpm)*TPM_m_un;

%----------------------
% %_____________

% ofd=0.0001;
%
% TPM_m =(1-(Nm)*ofd)*eye(Nm);
%
% TPM_m = TPM_m + ofd*ones(Nm,Nm);

% %____________

%===================================

ph2=zeros(1,Nm);
pb=zeros(1,Nm);
L=zeros(1,Nm);

%_____________

phtt=zeros(Nm,N);
phttT=zeros(Nm,N);
%______________

%================================================Iteration

iteration=1;

%_______________

c_JPM=zeros(1,N);

D_JPM_MC = zeros(iteration,N);
%_______________

%-----------------------
%++++++++++++++++++++++++++++

Row_D = zeros(2,N);

SEs1p = zeros(1,N);
SEs2p = zeros(1,N);

SEph1 = zeros(1,N);
SEvh1 = zeros(1,N);

SEph1_nr = zeros(1,N);
SEs1p_nr = zeros(1,N);

SEph2 = zeros(1,N);
SEvh2 = zeros(1,N);


%++++++++++++++++++++++++++++

% SEp1=zeros(iteration,N);
% SEv1=zeros(iteration,N);
% 
% SEp2=zeros(iteration,N);
% SEv2=zeros(iteration,N);
% 
% SEs1p=zeros(iteration,N);
% SEs2p=zeros(iteration,N);
% 
% SEph1=zeros(iteration,N);
% SEvh1=zeros(iteration,N);
% 
% SEph2=zeros(iteration,N);
% SEvh2=zeros(iteration,N);

TraceP1 = zeros(1,N);
TraceP2 = zeros(2,N);

pHkp1_t = zeros(2,N);

Sn = zeros(1,N);

tk1 = 10;
tk2 = 21;

%-----------------------

tic

for iter=1:iteration

    iter

    %---------------------------Trajectories

    %++++++++++++++ 1st
%         X1t0 = [75;-3;95;-2];
%         X2t0 = [30;3.1;20;3.2];
    %++++++++++++++

    %++++++++++++++2nd
    X1t0 = [75;-3;95;-2];
    X2t0 = [-61;15.3;-7;3.2];
    %++++++++++++++


    for i=0:9
        if i==0

%             v=D*randn(4,1);

            X1t(:,1)=F1*X1t0;%+v;
            X2t(:,1)=F1*X2t0;%+v;

        else

%             v=D*randn(4,1);

            X1t(:,i+1)=F1*X1t(:,i);%+v;
            X2t(:,i+1)=F1*X2t(:,i);%+v;

        end

    end

    for i=10:25

%         v=D*randn(4,1);

        X1t(:,i+1)=F1*X1t(:,i) + Fv(1:4,15);%+v;
        X2t(:,i+1)=F1*X2t(:,i) + Fv(5:8,15);%+v;

    end

    for i=26:N-1

%         v=D*randn(4,1);

        X1t(:,i+1)=F1*X1t(:,i);% + Fv(1:4,25);%+v;
        X2t(:,i+1)=F1*X2t(:,i);% + Fv(5:8,25);%+v;

    end

    %======================= End of Trajectory

    %=======================Measurements Noise

    sigmaz=0.2;          % also 1 , 0.1
    R=sigmaz^2*eye(Ns);

    %=======================

    %======================= Initialization Z0

    e_0=[0.2 0.5;0.5 0.2];

    %____________

    sigmaInip1=5;
    sigmaIniv1=0.5;

    x1_0 = X1t0 + diag([sigmaInip1 sigmaIniv1 sigmaInip1 sigmaIniv1])*randn(4,1);

    sigmaInip2=5;
    sigmaIniv2=0.5;

    x2_0 = [10;2;10;2] + diag([sigmaInip2 sigmaIniv2 sigmaInip2 sigmaIniv2])*randn(4,1);

    %______________

    X1_0=x1_0;

    P1x_0 = [sigmaInip1^2 0.01 0 0;
        0.01 sigmaIniv1^2 0 0;
        0 0 sigmaInip2^2 0.01;
        0 0 0.01 sigmaIniv2^2];

    X2_0=[x1_0;x2_0];

    P2x_0 = blkdiag(P1x_0,P1x_0);

    %----------------------------

    for im=1:Nm

        X_mk_Hk_1(:,im)=X1_0;
        P_mk_Hk_1(:,:,im) = P1x_0;

        X_mk_Hk_2(:,im)=X2_0;
        P_mk_Hk_2(:,:,im) = P2x_0;

    end

    pHk=[0.5 0.5];

    pmk_Hk_1 = 0.01*ones(1,Nm);
    pmk_Hk_2 = 0.01*ones(1,Nm);
    pmk_Hk_1(1,41) = pmk_Hk_1(1,41) + 1 - 0.01*Nm;
    pmk_Hk_2(1,41) = pmk_Hk_2(1,41) + 1 - 0.01*Nm;

    e_k = e_0;

    c_k = alpha.*c + beta.*e_k;

    % Cz_k=Cz_0;

    %__________________________

    dd=1;
    DDa=1;
    sa=DDa;

    Xh1=zeros(4,N);
    Xh2=zeros(8,N);

    Xhh1=zeros(4,N);
    Xhh2=zeros(8,N);

    S_hat1=zeros(1,N);
    S_hat2=zeros(2,N);

    S_hat1p=zeros(1,N);
    S_hat2p=zeros(2,N);

    d_JPM=zeros(1,N);

    Dei = zeros(1,N);

    Di = zeros(1,N);
    
    %============================================================ Zk to Zk+1

    for tk=1:N

        %         tk
        %---------------------Measurements

        %---------------------------------Gaussian

        varss1=3500;  % also 1 , 0.1
        varss2=3500;  % also 1 , 0.1

        %---------------for data generation
        Cst = [varss1 0;
            0 varss2];

        us1t=220;    %  400
        us2t=220;    %  400
        %--------------
        %-------------- for processing

        Cs = [varss1/(sa^2) 0;
            0 varss2/(sa^2)];

        us1=us1t/sa;
        us2=us2t/sa;

        us=[us1;us2];

        %-------------

        Ss1t(1,tk) = us1t + sqrt(varss1)*randn;
        Ss2t(1,tk) = us2t + sqrt(varss2)*randn;

        %----------------------------------

        %---------------------Measurements

        if tk<tk1

            Sn(1,tk) = 1;
            for im=1:Ns

                Z(im,tk) =  Ss1t(1,tk)*(1/sqrt((X1t(1,tk)-xs(im,1))^2+(X1t(3,tk)-ys(im,1))^2+dd)) +  ...
                    sigmaz*randn(1,1);

            end

        elseif tk>tk1-1 && tk<tk2

            Sn(1,tk) = 2;

            for im=1:Ns

                Z(im,tk) =  Ss1t(1,tk)*(1/sqrt((X1t(1,tk)-xs(im,1))^2+(X1t(3,tk)-ys(im,1))^2+dd)) +  ...
                    Ss2t(1,tk)*(1/sqrt((X2t(1,tk)-xs(im,1))^2+(X2t(3,tk)-ys(im,1))^2+dd)) + ...
                    sigmaz*randn(1,1);

            end

        elseif tk>tk2-1

            Sn(1,tk) = 1;

            for im=1:Ns

                Z(im,tk) =  Ss1t(1,tk)*(1/sqrt((X1t(1,tk)-xs(im,1))^2+(X1t(3,tk)-ys(im,1))^2+dd)) +  ...
                    sigmaz*randn(1,1);

            end

        end


        Zk=Z(:,tk);

        %_________________________

        %============================= XHkp1 and pHz_kp1 part

        for i_Hkp1=1:2

            %______________

            f_Hkp1(1,i_Hkp1) = 0;
            %______________

            pHkp1_k(1,i_Hkp1) = sum(TPM_H(:,i_Hkp1)'.*pHk(1,:));

            %________________

            c1 = sum(TPM_H(:,i_Hkp1)'.*pHk(1,:));
            pHk_Hkp1(1,:) = 1/c1*TPM_H(:,i_Hkp1)'.*pHk(1,:);

            %________________

            if i_Hkp1==1

                for i_mkp1=1:Nm

                    %_____________

                    pmkp1_k = 0;

                    f_Hkp1mkp1 = 0;
                    %_____________

                    for i_Hk=1:2

                        if i_Hk==1

                            %_________________

                            c1 = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:));
                            pmk_mkp1(1,:) = 1/c1*TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:);

                            M_pmk_mkp1 = zeros(4,Nm);
                            M_pmk_mkp1(1,:) = pmk_mkp1(1,:);
                            M_pmk_mkp1 = cumsum(M_pmk_mkp1,1);
                            %_________________

                            pmkp1_Hk_k = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:));
                            %_________________

                            %++++++++++++++++++++++++++++++++++++ Updating
                            %for Hkp1, Hk, mkp1, mk

                            for i_mk = 1:Nm

                                Xp1 = F1*X_mk_Hk_1(:,i_mk) + Fv(1:4,i_mkp1);
                                Pp1 = F1*P_mk_Hk_1(:,:,i_mk)*F1' + G1*Qd*G1';

                                %-----------------
                                for im=1:Ns

                                    Mb1(im,1)=(DDa/sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd));

                                end

                                %-----------------------

                                for im=1:Ns

                                    Hh1(im,1)=-DDa*(Xp1(1,1)-xs(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                                    Hh1(im,2)=0;
                                    Hh1(im,3)=-DDa*(Xp1(3,1)-ys(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                                    Hh1(im,4)=0;

                                end

                                %-----------------------------------------

                                Covh1 = Hh1*Pp1*Hh1';
                                vars1 = Cs(1,1);
                                mus1 = us(1,1);
                                muh1 = Mb1;
                                Zkp = muh1*us(1,1);

                                CvZk = Covh1*vars1 + Covh1*mus1^2 + muh1*muh1'*vars1 + R;

                                I_CvZk=inv(CvZk);

                                f_Hkp1_mkp1_Hk_mk(1,i_mk) = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));

                            end

                            f_Hkp1mkp1Hk = sum(f_Hkp1_mkp1_Hk_mk.*pmk_mkp1);

                            %=================================================== I

                            XM1 = sum(X_mk_Hk_1.*M_pmk_mkp1,2);
                            %____________

                            PM1 = zeros(4,4);
                            %____________

                            for i=1:Nm

                                PM1 = PM1 + pmk_mkp1(1,i)*(P_mk_Hk_1(:,:,i)+(X_mk_Hk_1(:,i)-XM1)*(X_mk_Hk_1(:,i)-XM1)');

                            end

                            Xp1 = F1*XM1 + Fv(1:4,i_mkp1);
                            Pp1 = F1*PM1*F1'+ G1*Qd*G1';

                            Xp1_Hkp1mkp1Hk(:,i_Hk) = Xp1;
                            Pp1_Hkp1mkp1Hk(:,:,i_Hk) = Pp1;

                            %                             %=================================================== II
                            %
                            %                             %-----------------------Zkp1 prediction
                            %                             %----------------------S Estimation
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Mb1(im,1)=(DDa/sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd));
                            %
                            %                             end
                            %
                            %
                            %                             %-----------------------
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Hh1(im,1)=-DDa*(Xp1(1,1)-xs(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh1(im,2)=0;
                            %                                 Hh1(im,3)=-DDa*(Xp1(3,1)-ys(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh1(im,4)=0;
                            %
                            %                             end
                            %
                            %                             %-----------------------------------------
                            %
                            %                             Covh1 = Hh1*Pp1*Hh1';
                            %                             vars1 = Cs(1,1);
                            %                             mus1 = us(1,1);
                            %                             muh1 = Mb1;
                            %                             Zkp = muh1*us(1,1);
                            %
                            %                             CvZk = Covh1*vars1 + Covh1*mus1^2 + muh1*muh1'*vars1 + R;
                            %
                            %                             I_CvZk=inv(CvZk);
                            %
                            %                             f_Hkp1mkp1Hk = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));
                            %
                            %                             %----------------------------

                            %=================================================

                        else % else of i_Hk

                            %_________________

                            c1 = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:));
                            pmk_mkp1(1,:) = 1/c1*TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:);

                            M_pmk_mkp1 = zeros(4,Nm);
                            M_pmk_mkp1(1,:) = pmk_mkp1(1,:);
                            M_pmk_mkp1 = cumsum(M_pmk_mkp1,1);
                            %_________________

                            pmkp1_Hk_k = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:));
                            %_________________

                            %++++++++++++++++++++++++++++++++++++ Updating
                            %for Hkp1, Hk, mkp1, mk

                            for i_mk = 1:Nm

                                Xpd1 = F1*X_mk_Hk_2(1:4,i_mk) + Fv(1:4,i_mkp1);
                                Ppd1 = F1*P_mk_Hk_2(1:4,1:4,i_mk)*F1' + G1*Qd*G1';

                                Xpd2 = F1*X_mk_Hk_2(5:8,i_mk) + Fv(5:8,i_mkp1);
                                Ppd2 = F1*P_mk_Hk_2(5:8,5:8,i_mk)*F1' + G1*Qd*G1';

                                %+++++++++++++++++++++++++++++

                                ini11min = Xpd1(1,1)>0;
                                ini11max = Xpd1(1,1)<100;
                                ini13min = Xpd1(3,1)>0;
                                ini13max = Xpd1(3,1)<100;
                                ini21min = Xpd2(1,1)>0;
                                ini21max = Xpd2(1,1)<100;
                                ini23min = Xpd2(3,1)>0;
                                ini23max = Xpd2(3,1)<100;

                                ini = ini11min + ini11max + ini13min + ini13max + ini21min + ini21max + ini23min + ini23max;
                                if  ini==8

                                    d1 = sqrt((Xpd1(1,1)-100)^2 + (Xpd1(3,1)-0)^2);
                                    d2 = sqrt((Xpd2(1,1)-100)^2 + (Xpd2(3,1)-0)^2);
                                    W1 = d1/(d1 + d2);
                                    W2 = d2/(d1 + d2);

                                else

                                    d1 = sqrt((Xpd1(1,1)-100)^2 + (Xpd1(3,1)-0)^2);
                                    d2 = sqrt((Xpd2(1,1)-100)^2 + (Xpd2(3,1)-0)^2);
                                    
                                    if d2<10
                                        %.....................
%                                         x2_0 = [10;2;10;2] + diag([sigmaInip2 sigmaIniv2 sigmaInip2 sigmaIniv2])*randn(4,1);
%                                         Xpd2 = x2_0;
                                        
%                                         for im=1:Nm
%                                             
%                                             X_mk_Hk_2(5:8,im)=x2_0;
%                                             P_mk_Hk_2(5:8,5:8,im) = P1x_0;
%                                             
%                                         end
                                        
                                        %.....................
%                                         d2 = sqrt((Xpd2(1,1)-100)^2 + (Xpd2(3,1)-0)^2);
                                        
                                        W1 = d1/(d1 + d2);
                                        W2 = d2/(d1 + d2);
                                        
                                    else
                                    
                                        W1 = 1;
                                        W2 = 0;
                                        
                                    end
                                    
                                end

                                %+++++++++++++++++++++++++++++

                                Xp1 = W1*Xpd1 + W2*Xpd2;
                                Pp1 = W1*Ppd1 + W2*Ppd2;

                                %-----------------
                                for im=1:Ns

                                    Mb1(im,1)=(DDa/sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd));

                                end

                                %-----------------------

                                for im=1:Ns

                                    Hh1(im,1)=-DDa*(Xp1(1,1)-xs(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                                    Hh1(im,2)=0;
                                    Hh1(im,3)=-DDa*(Xp1(3,1)-ys(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                                    Hh1(im,4)=0;

                                end

                                %-----------------------------------------

                                Covh1 = Hh1*Pp1*Hh1';
                                vars1 = Cs(1,1);
                                mus1 = us(1,1);
                                muh1 = Mb1;
                                Zkp = muh1*us(1,1);

                                CvZk = Covh1*vars1 + Covh1*mus1^2 + muh1*muh1'*vars1 + R;

                                I_CvZk=inv(CvZk);

                                f_Hkp1_mkp1_Hk_mk(1,i_mk) = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));

                            end

                            f_Hkp1mkp1Hk = sum(f_Hkp1_mkp1_Hk_mk.*pmk_mkp1);

                            %----------------------------
                            %===================================================I

                            %+++++++++++++ Equaly likely disappearance

                            %---------------- 2nd disappeared

                            XMm1 = sum(X_mk_Hk_2(1:4,:).*M_pmk_mkp1,2);
                            %____________

                            PMm1 = zeros(4,4);
                            %____________

                            for i=1:Nm

                                PMm1 = PMm1 + pmk_mkp1(1,i)*(P_mk_Hk_2(1:4,1:4,i)+(X_mk_Hk_2(1:4,i)-XMm1)*(X_mk_Hk_2(1:4,i)-XMm1)');

                            end
                            Xpm1 = F1*XMm1 + Fv(1:4,i_mkp1);
                            Ppm1 = F1*PMm1*F1'+ G1*Qd*G1';

                            %----------------- 1st disappeared

                            XMm2 = sum(X_mk_Hk_2(5:8,:).*M_pmk_mkp1,2);
                            %____________

                            PMm2 = zeros(4,4);
                            %____________

                            for i=1:Nm

                                PMm2 = PMm2 + pmk_mkp1(1,i)*(P_mk_Hk_2(5:8,5:8,i)+(X_mk_Hk_2(5:8,i)-XMm2)*(X_mk_Hk_2(5:8,i)-XMm2)');

                            end
                            Xpm2 = F1*XMm2 + Fv(5:8,i_mkp1);
                            Ppm2 = F1*PMm2*F1'+ G1*Qd*G1';
                            %+++++++++++++++++++++++++++++

                            ini11min = Xpm1(1,1)>0;
                            ini11max = Xpm1(1,1)<100;
                            ini13min = Xpm1(3,1)>0;
                            ini13max = Xpm1(3,1)<100;
                            ini21min = Xpm2(1,1)>0;
                            ini21max = Xpm2(1,1)<100;
                            ini23min = Xpm2(3,1)>0;
                            ini23max = Xpm2(3,1)<100;

                            ini = ini11min + ini11max + ini13min + ini13max + ini21min + ini21max + ini23min + ini23max;
                            if  ini==8

                                d1 = sqrt((Xpm1(1,1)-100)^2 + (Xpm1(3,1)-0)^2);
                                d2 = sqrt((Xpm2(1,1)-100)^2 + (Xpm2(3,1)-0)^2);
                                W11 = d1/(d1 + d2);
                                W22 = d2/(d1 + d2);

                            else
                                
                                d1 = sqrt((Xpd1(1,1)-100)^2 + (Xpd1(3,1)-0)^2);
                                d2 = sqrt((Xpd2(1,1)-100)^2 + (Xpd2(3,1)-0)^2);
                                
                                if d2<10
                                    %.....................
                                    %      x2_0 = [10;2;10;2] + diag([sigmaInip2 sigmaIniv2 sigmaInip2 sigmaIniv2])*randn(4,1);
                                    %        Xpd2 = x2_0;
                                    
                                    %     for im=1:Nm
                                    %
                                    %          X_mk_Hk_2(5:8,im)=x2_0;
                                    %         P_mk_Hk_2(5:8,5:8,im) = P1x_0;
                                    %
                                    %                 end
                                    
                                    %.....................
                                    %          d2 = sqrt((Xpd2(1,1)-100)^2 + (Xpd2(3,1)-0)^2);
                                    
                                    W11 = d1/(d1 + d2);
                                    W22 = d2/(d1 + d2);
                                    
                                else
                                    
                                    W11 = 1;
                                    W22 = 0;
                                    
                                end
                                
                               
                            end

                            %+++++++++++++++++++++++++++++

                            Xp1 = W11*Xpm1 + W22*Xpm2;
                            Pp1 = W11*Ppm1 + W22*Ppm2;

                            Xp1_Hkp1mkp1Hk(:,i_Hk) = Xp1;
                            Pp1_Hkp1mkp1Hk(:,:,i_Hk) = Pp1;

                            %                             %===================================================II
                            %
                            %                             %-----------------------Zkp1 prediction
                            %                             %----------------------S Estimation
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Mb1(im,1)=(DDa/sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd));
                            %
                            %                             end
                            %
                            %                             %-----------------------
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Hh1(im,1)=-DDa*(Xp1(1,1)-xs(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh1(im,2)=0;
                            %                                 Hh1(im,3)=-DDa*(Xp1(3,1)-ys(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh1(im,4)=0;
                            %
                            %                             end
                            %
                            %                             %--------------------------------------------
                            %
                            %                             Covh1 = Hh1*Pp1*Hh1';
                            %                             vars1 = Cs(1,1);
                            %                             mus1 = us(1,1);
                            %                             muh1 = Mb1;
                            %                             Zkp = muh1*us(1,1);
                            %
                            %                             CvZk = Covh1*vars1 + Covh1*mus1^2 + muh1*muh1'*vars1 + R;
                            %
                            %                             I_CvZk=inv(CvZk);
                            %
                            %                             f_Hkp1mkp1Hk = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));
                            %
                            %                             %----------------------------
                            %                             %=================================================

                        end % End of i_Hk if

                        %-----------------------------II

                        f_Hkp1mkp1 = f_Hkp1mkp1 + f_Hkp1mkp1Hk*pHk_Hkp1(1,i_Hk);

                        %-----------------------------III

                        pmkp1_k = pmkp1_k + pmkp1_Hk_k*pHk_Hkp1(1,i_Hk);

                        %-----------------------------
                        %_______________________

                    end % end of i_Hk for

                    %______________

                    f_Hkp1(1,i_Hkp1) = f_Hkp1(1,i_Hkp1) + f_Hkp1mkp1*pmkp1_k;

                    %______________

                    %-----------------------------I
                    %_______________

                    Xp1_Hkp1mkp1 = zeros(4,1);
                    %_______________

                    for i=1:2

                        Xp1_Hkp1mkp1 = Xp1_Hkp1mkp1 + Xp1_Hkp1mkp1Hk(:,i)*pHk_Hkp1(1,i);

                    end

                    %____________

                    Pp1_Hkp1mkp1 = zeros(4,4);
                    %____________

                    for i=1:2

                        Pp1_Hkp1mkp1 = Pp1_Hkp1mkp1 +...
                            pHk_Hkp1(1,i)*...
                            (Pp1_Hkp1mkp1Hk(:,:,i)+(Xp1_Hkp1mkp1Hk(:,i)-Xp1_Hkp1mkp1)*(Xp1_Hkp1mkp1Hk(:,i)-Xp1_Hkp1mkp1)');

                    end

                    %-----------------------------------------End of I

                    %===========================for JPM
                    %===========================

                    if tk>1

                        if Di(1,tk-1)==1

                            Xp1_Hkp1mkp1_JPM(:,i_mkp1) = Xp1_Hkp1mkp1;
                            Pp1_Hkp1mkp1_JPM(:,:,i_mkp1) = Pp1_Hkp1mkp1;

                            pmkp1_k_JPM(1,i_mkp1) = pmkp1_k;

                        end
                    end

                    %===========================

                    %==============================Updating I

                    Xp1 = Xp1_Hkp1mkp1;
                    Pp1 = Pp1_Hkp1mkp1;

                    for im=1:Ns

                        Mb1(im,1)=(DDa/sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd));

                    end

                    %----------------------

                    %-----------------------

                    for im=1:Ns

                        Hh1(im,1)=-DDa*(Xp1(1,1)-xs(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                        Hh1(im,2)=0;
                        Hh1(im,3)=-DDa*(Xp1(3,1)-ys(im,1))/(sqrt((Xp1(1,1)-xs(im,1))^2+(Xp1(3,1)-ys(im,1))^2+dd))^3;
                        Hh1(im,4)=0;

                    end

                    %-----------------------

                    Covh1 = Hh1*Pp1*Hh1';
                    vars1 = Cs(1,1);
                    mus1 = us(1,1);
                    muh1 = Mb1 ;
                    Zkp1 = muh1*us(1,1);

                    Poo1 = Covh1*vars1 + Covh1*mus1^2 + muh1*muh1'*vars1 + R;
                    Pxo1 = Pp1*Hh1'*us(1,1);
                    Pox1 = Pxo1';

                    X1_Hkp1mkp1(:,i_mkp1) = Xp1 + Pxo1*inv(Poo1)*(Zk - Zkp1);
                    P1_Hkp1mkp1(:,:,i_mkp1) = Pp1 - Pxo1*inv(Poo1)*Pox1;

                    %----------------------------

                    %======================================== Updating II

                    %----------------------------

                    pmkp1_un(1,i_mkp1) = f_Hkp1mkp1*pmkp1_k;

                    %----------------------------

                    %========================================

                end % end of i_mkp1

                %================================II

                c1=sum(pmkp1_un);
                pmkp1 = pmkp1_un/c1;

                %================================
                %_____________

                X1_Hkp1 = zeros(4,1);
                %_____________

                for i=1:Nm

                    X1_Hkp1 = X1_Hkp1 + X1_Hkp1mkp1(:,i)*pmkp1(1,i);

                end

                %____________

                P1_Hkp1 = zeros(4,4);
                %____________

                for i=1:Nm

                    P1_Hkp1 = P1_Hkp1 +...
                        pmkp1(1,i)*...
                        (P1_Hkp1mkp1(:,:,i)+(X1_Hkp1mkp1(:,i)-X1_Hkp1)*(X1_Hkp1mkp1(:,i)-X1_Hkp1)');

                end

                %____________
                %
                %=================================== Start of JPM
                %===================================

                if tk>1

                    if Di(1,tk-1)==1
                        %_______________

                        Xp1_Hkp1_JPM = zeros(4,1);
                        %_______________

                        for i=1:Nm

                            Xp1_Hkp1_JPM = Xp1_Hkp1_JPM + Xp1_Hkp1mkp1_JPM(:,i)*pmkp1_k_JPM(1,i);

                        end

                        %____________

                        Pp1_Hkp1_JPM = zeros(4,4);
                        %____________

                        for i=1:Nm

                            Pp1_Hkp1_JPM = Pp1_Hkp1_JPM +...
                                pmkp1_k_JPM(1,i)*...
                                (Pp1_Hkp1mkp1_JPM(:,:,i)+(Xp1_Hkp1mkp1_JPM(:,i)-Xp1_Hkp1_JPM)*(Xp1_Hkp1mkp1_JPM(:,i)-Xp1_Hkp1_JPM)');

                        end

                        %------------------------------

                        for im=1:Ns

                            Mb1_JPM(im,1)=(DDa/sqrt((Xp1_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp1_Hkp1_JPM(3,1)-ys(im,1))^2+dd));

                        end

                        %___________________

                        for im=1:Ns

                            Hh1_JPM(im,1)=-DDa*(Xp1_Hkp1_JPM(1,1)-xs(im,1))/(sqrt((Xp1_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp1_Hkp1_JPM(3,1)-ys(im,1))^2+dd))^3;
                            Hh1_JPM(im,2)=0;
                            Hh1_JPM(im,3)=-DDa*(Xp1_Hkp1_JPM(3,1)-ys(im,1))/(sqrt((Xp1_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp1_Hkp1_JPM(3,1)-ys(im,1))^2+dd))^3;
                            Hh1_JPM(im,4)=0;

                        end

                        %-----------------------

                        Covh1_JPM = Hh1_JPM*Pp1_Hkp1_JPM*Hh1_JPM';
                        mus1_JPM = us(1,1);
                        vars1_JPM = Cs(1,1);
                        muh1_JPM = Mb1_JPM ;

                        PZkp1_k_JPM = Covh1_JPM*mus1_JPM^2 + Covh1_JPM*vars1_JPM + muh1_JPM*muh1_JPM'*vars1_JPM + R;

                        Zkp1_k_Mean = muh1_JPM*mus1_JPM;

                        %___________________
                        %------------------------------ Zkp1_k

                        Nj=10;

                        if det(PZkp1_k_JPM)<inf

                            Lamb=eig(PZkp1_k_JPM);
                            cmp=Lamb>0;
                            icmp=sum(cmp);

                            Dig = 0;
                            for ii=1:Ns

                                Dig = Dig + PZkp1_k_JPM(ii,ii);
                            end
                            IDig = imag(Dig);

                            if icmp==Ns && IDig==0

                                %----------------------
                                c_JPM(1,tk)=c_JPM(1,tk)+1;
                                %----------------------

                                D_PZkp1_k = chol(PZkp1_k_JPM);

                                Zkp1_k_JPM = diag(Zkp1_k_Mean)*ones(Ns,Nj) + D_PZkp1_k*randn(Ns,Nj);

                                %----------------------------- d calculation

                                Zk_JPM = diag(Zk)*ones(Ns,Nj);

                                d_JPM1 = sqrt(sum((Zk_JPM-Zkp1_k_JPM).^2,1));
                                d_JPM(1,tk) = (mean(d_JPM1));

% %                                 d_JPM1 = (mean((Zk_JPM-Zkp1_k_JPM).^2,2));
% %                                 d_JPM(1,tk) = sqrt(mean(d_JPM1));

                            end
                            %------------------------------

                        end

                    end
                end

                %================================
                %================================ End of JPM

            else % else of i_Hkp1

                for i_mkp1=1:Nm

                    %_____________

                    pmkp1_k = 0;

                    f_Hkp1mkp1 = 0;
                    %_____________

                    for i_Hk=1:2

                        if i_Hk==1

                            %_________________

                            c1 = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:));
                            pmk_mkp1(1,:) = 1/c1*TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:);

                            M_pmk_mkp1 = zeros(4,Nm);
                            M_pmk_mkp1(1,:) = pmk_mkp1(1,:);
                            M_pmk_mkp1 = cumsum(M_pmk_mkp1,1);
                            %_________________

                            pmkp1_Hk_k = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_1(1,:));
                            %_________________

                            %++++++++++++++++++++++++++++++++++++ Updating
                            %for Hkp1, Hk, mkp1, mk

                            for i_mk = 1:Nm

                                Xp21 = F1*X_mk_Hk_1(:,i_mk) + Fv(1:4,i_mkp1);
                                Pp21 = F1*P_mk_Hk_1(:,:,i_mk)*F1' + G1*Qd*G1';

                                %______________

                                if tk==1

%                                     sigmaInipNEW=50;
%                                     sigmaInivNEW=5;
                                    
                                    X2t0g = X2t0 + diag([sigmaInipNEW sigmaInivNEW sigmaInipNEW sigmaInivNEW])*randn(4,1);
                                    
                                    X_NEW = X2t0g;

                                else
                                    
%                                     sigmaInipNEW=50;
%                                     sigmaInivNEW=5;
                                    
                                    X2tNEWg = X2t(:,tk-1) + diag([sigmaInipNEW sigmaInivNEW sigmaInipNEW sigmaInivNEW])*randn(4,1);
                                    
                                    X_NEW = X2tNEWg;

%                                     X_NEW = X2t(:,tk-1);

                                end
                                %_________________

                                Xp22=F1*X_NEW;
                                Pp22=F1*P_X_NEW*F1'+ G1*Qd*G1';

                                %-----------------

                                Xp2 = [Xp21;Xp22];
                                Pp2 = blkdiag(Pp21,Pp22);

                                for im=1:Ns

                                    Mb2(im,1)=(DDa/sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd));

                                    Mb2(im,2)=(DDa/sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd));

                                end

                                %-----------------------

                                for im=1:Ns

                                    Hh21(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                                    Hh21(im,2)=0;
                                    Hh21(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                                    Hh21(im,4)=0;
                                    Hh21(im,5:8)=zeros(1,4);

                                    Hh22(im,1:4)=zeros(1,4);
                                    Hh22(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                                    Hh22(im,6)=0;
                                    Hh22(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                                    Hh22(im,8)=0;

                                end

                                %-----------------------------------------

                                Covh21 = Hh21*Pp2*Hh21' ;
                                Covh22 = Hh22*Pp2*Hh22' ;
                                Covhm12 = Hh21*Pp2*Hh22';
                                Covhm21 = Hh22*Pp2*Hh21';
                                vars1 = Cs(1,1);
                                vars2 = Cs(2,2);
                                mus1 = us(1,1);
                                mus2 = us(2,1);
                                muh21 = Mb2(:,1) ;
                                muh22 = Mb2(:,2) ;

                                Poo2 = Covh21*vars1 + Covh21*mus1^2 + muh21*muh21'*vars1 +  ...
                                    Covhm12*mus1*mus2 + Covhm21*mus1*mus2 + Covh22*vars2 + Covh22*mus2^2 + muh22*muh22'*vars2 + R;
                                CvZk = Poo2;

                                Zkp = (muh21)*us(1,1) + (muh22)*us(2,1);

                                I_CvZk=inv(CvZk);

                                f_Hkp1_mkp1_Hk_mk(1,i_mk) = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));

                            end

                            f_Hkp1mkp1Hk = sum(f_Hkp1_mkp1_Hk_mk.*pmk_mkp1);

                            %----------------------------

                            %=================================================== I

                            %...................................... x1 part

                            XM21 = sum(X_mk_Hk_1.*M_pmk_mkp1,2);
                            %____________

                            PM21 = zeros(4,4);
                            %____________

                            for i=1:Nm

                                PM21 = PM21 + pmk_mkp1(1,i)*(P_mk_Hk_1(:,:,i)+(X_mk_Hk_1(:,i)-XM21)*(X_mk_Hk_1(:,i)-XM21)');

                            end

                            Xp21 = F1*XM21 + Fv(1:4,i_mkp1);
                            Pp21 = F1*PM21*F1'+ G1*Qd*G1';

                            %.................................... x2 part

                            %_________________

                            if tk==1
                                
%                                 sigmaInipNEW=50;
%                                 sigmaInivNEW=5;
                                
                                X2t0g = X2t0 + diag([sigmaInipNEW sigmaInivNEW sigmaInipNEW sigmaInivNEW])*randn(4,1);
                                
                                X_NEW = X2t0g;
                                
                            else
                                
%                                 sigmaInipNEW=50;
%                                 sigmaInivNEW=5;
                                
                                X2tNEWg = X2t(:,tk-1) + diag([sigmaInipNEW sigmaInivNEW sigmaInipNEW sigmaInivNEW])*randn(4,1);
                                
                                X_NEW = X2tNEWg;
                                
                                %                                     X_NEW = X2t(:,tk-1);
                                
                            end
                            %_________________

                            XM22 = X_NEW;  % at this state at k
                            PM22 = P_X_NEW;

                            Xp22=F1*XM22;
                            Pp22=F1*PM22*F1'+ G1*Qd*G1';

                            %..........................

                            Xp2 = [Xp21;Xp22];
                            Pp2 = [Pp21 zeros(4,4);
                                zeros(4,4) Pp22];

                            Xp2_Hkp1mkp1Hk(:,i_Hk) = Xp2;
                            Pp2_Hkp1mkp1Hk(:,:,i_Hk) = Pp2;

                            %....................................

                            %                             %=================================================== II
                            %
                            %                             %----------------------S Estimation
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Mb2(im,1)=(DDa/sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd));
                            %
                            %                                 Mb2(im,2)=(DDa/sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd));
                            %
                            %                             end
                            %
                            %                             %-----------------------
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Hh21(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh21(im,2)=0;
                            %                                 Hh21(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh21(im,4)=0;
                            %                                 Hh21(im,5:8)=zeros(1,4);
                            %
                            %                                 Hh22(im,1:4)=zeros(1,4);
                            %                                 Hh22(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh22(im,6)=0;
                            %                                 Hh22(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh22(im,8)=0;
                            %
                            %                             end
                            %
                            %                             %----------------------------
                            %
                            %                             Covh21 = Hh21*Pp2*Hh21';
                            %                             Covh22 = Hh22*Pp2*Hh22';
                            %                             Covhm12 = Hh21*Pp2*Hh22';
                            %                             Covhm21 = Hh22*Pp2*Hh21';
                            %                             vars1 = Cs(1,1);
                            %                             vars2 = Cs(2,2);
                            %                             mus1 = us(1,1);
                            %                             mus2 = us(2,1);
                            %                             muh21 = Mb2(:,1);
                            %                             muh22 = Mb2(:,2);
                            %
                            %                             Poo2 = Covh21*vars1 + Covh21*mus1^2 + muh21*muh21'*vars1 +  ...
                            %                                 Covhm12*mus1*mus2 + Covhm21*mus1*mus2 + Covh22*vars2 + Covh22*mus2^2 + muh22*muh22'*vars2 + R;
                            %                             CvZk = Poo2;
                            %                             %
                            %                             %                             %                             Pxo = Pp2*Hh2(:,:,1)*us(1,1) + Pp2*Hh2(:,:,2)*us(2,1);
                            %                             %                             %                             Pox = Pxo';
                            %                             Zkp = (muh21)*us(1,1) + (muh22)*us(2,1);
                            %
                            %                             I_CvZk=inv(CvZk);
                            %
                            %                             f_Hkp1mkp1Hk = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));
                            %
                            %                             %----------------------------
                            %
                            %                             %=================================================

                        else % else of i_Hk

                            %_________________

                            c1 = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:));
                            pmk_mkp1(1,:) = 1/c1*TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:);

                            M_pmk_mkp1 = zeros(8,Nm);
                            M_pmk_mkp1(1,:) = pmk_mkp1(1,:);
                            M_pmk_mkp1 = cumsum(M_pmk_mkp1,1);
                            %_________________

                            pmkp1_Hk_k = sum(TPM_m(:,i_mkp1)'.*pmk_Hk_2(1,:));
                            %_________________

                            %++++++++++++++++++++++++++++++++++++ Updating
                            %for Hkp1, Hk, mkp1, mk

                            for i_mk =1:Nm

                                Xp2 = F2*X_mk_Hk_2(:,i_mk) + Fv(:,i_mkp1);
                                Pp2 = F2*P_mk_Hk_2(:,:,i_mk)*F2' + G2*Qd*G2';


                                for im=1:Ns

                                    Mb2(im,1)=(DDa/sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd));

                                    Mb2(im,2)=(DDa/sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd));

                                end

                                %-----------------------

                                for im=1:Ns

                                    Hh21(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                                    Hh21(im,2)=0;
                                    Hh21(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                                    Hh21(im,4)=0;
                                    Hh21(im,5:8)=zeros(1,4);

                                    Hh22(im,1:4)=zeros(1,4);
                                    Hh22(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                                    Hh22(im,6)=0;
                                    Hh22(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                                    Hh22(im,8)=0;

                                end

                                %-----------------------------------------

                                Covh21 = Hh21*Pp2*Hh21' ;
                                Covh22 = Hh22*Pp2*Hh22' ;
                                Covhm12 = Hh21*Pp2*Hh22';
                                Covhm21 = Hh22*Pp2*Hh21';
                                vars1 = Cs(1,1);
                                vars2 = Cs(2,2);
                                mus1 = us(1,1);
                                mus2 = us(2,1);
                                muh21 = Mb2(:,1) ;
                                muh22 = Mb2(:,2) ;

                                Poo2 = Covh21*vars1 + Covh21*mus1^2 + muh21*muh21'*vars1 +  ...
                                    Covhm12*mus1*mus2 + Covhm21*mus1*mus2 + Covh22*vars2 + Covh22*mus2^2 + muh22*muh22'*vars2 + R;
                                CvZk = Poo2;

                                Zkp = (muh21)*us(1,1) + (muh22)*us(2,1);

                                I_CvZk=inv(CvZk);

                                f_Hkp1_mkp1_Hk_mk(1,i_mk) = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));

                            end

                            f_Hkp1mkp1Hk = sum(f_Hkp1_mkp1_Hk_mk.*pmk_mkp1);

                            %----------------------------

                            %===================================================I

                            %..................................... x1 part

                            XM2 = sum(X_mk_Hk_2(:,:).*M_pmk_mkp1,2);
                            %____________

                            PM2 = zeros(8,8);
                            %____________

                            for i=1:Nm

                                PM2 = PM2 + pmk_mkp1(1,i)*(P_mk_Hk_2(:,:,i)+(X_mk_Hk_2(:,i)-XM2)*(X_mk_Hk_2(:,i)-XM2)');

                            end

                            Xp2 = F2*XM2 + Fv(:,i_mkp1);
                            Pp2 = F2*PM2*F2' + G2*Qd*G2';

                            %                             %..................................... x1 part
                            %
                            %                             XM21 = sum(X_mk_Hk_2(1:4,:).*M_pmk_mkp1,2);
                            %                             %____________
                            %
                            %                             PM21 = zeros(4,4);
                            %                             %____________
                            %
                            %                             for i=1:Nm
                            %
                            %                                 PM21 = PM21 + pmk_mkp1(1,i)*(P_mk_Hk_2(1:4,1:4,i)+(X_mk_Hk_2(1:4,i)-XM21)*(X_mk_Hk_2(1:4,i)-XM21)');
                            %
                            %                             end
                            %
                            %                             Xp21 = F1*XM21 + Fv(1:4,i_mkp1);
                            %                             Pp21 = F1*PM21*F1'+Q1;
                            %
                            %                             %................................... x2 part
                            %
                            %                             XM22 = sum(X_mk_Hk_2(5:8,:).*M_pmk_mkp1,2);
                            %                             %____________
                            %
                            %                             PM22 = zeros(4,4);
                            %                             %____________
                            %
                            %                             for i=1:Nm
                            %
                            %                                 PM22 = PM22 + pmk_mkp1(1,i)*(P_mk_Hk_2(5:8,5:8,i)+(X_mk_Hk_2(5:8,i)-XM22)*(X_mk_Hk_2(5:8,i)-XM22)');
                            %
                            %                             end
                            %
                            %                             Xp22 = F1*XM22 + Fv(5:8,i_mkp1);
                            %                             Pp22 = F1*PM22*F1'+Q1;
                            %
                            %                             %..........................
                            %
                            %                             Xp2 = [Xp21;Xp22];
                            %                             Pp2 = [Pp21 zeros(4,4);
                            %                                 zeros(4,4) Pp22];

                            Xp2_Hkp1mkp1Hk(:,i_Hk) = Xp2;
                            Pp2_Hkp1mkp1Hk(:,:,i_Hk) = Pp2;

                            %....................................

                            %                             %===================================================II
                            %
                            %                             %----------------------S Estimation
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Mb2(im,1)=(DDa/sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd));
                            %
                            %                                 Mb2(im,2)=(DDa/sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd));
                            %
                            %                             end
                            %
                            %                             %----------------------
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Hh21(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh21(im,2)=0;
                            %                                 Hh21(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh21(im,4)=0;
                            %                                 Hh21(im,5:8)=zeros(1,4);
                            %
                            %                                 Hh22(im,1:4)=zeros(1,4);
                            %                                 Hh22(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh22(im,6)=0;
                            %                                 Hh22(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh22(im,8)=0;
                            %
                            %                             end
                            %
                            %                             for im=1:Ns
                            %
                            %                                 Hh2(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh2(im,2)=0;
                            %                                 Hh2(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh2(im,4)=0;
                            %                                 Hh2(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh2(im,6)=0;
                            %                                 Hh2(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                            %                                 Hh2(im,8)=0;
                            %
                            %                             end
                            %
                            %                             %----------------------------
                            %
                            %                             Covh21 = Hh21*Pp2*Hh21';
                            %                             Covh22 = Hh22*Pp2*Hh22';
                            %                             Covhm21 = Hh21*Pp2*Hh22';
                            %                             Covhm12 = Hh22*Pp2*Hh21';
                            %                             vars1 = Cs(1,1);
                            %                             vars2 = Cs(2,2);
                            %                             mus1 = us(1,1);
                            %                             mus2 = us(2,1);
                            %                             muh21 = Mb2(:,1);
                            %                             muh22 = Mb2(:,2);
                            %
                            %                             Poo2 = Covh21*vars1 + Covh21*mus1^2 + muh21*muh21'*vars1 + Covhm12*mus1*mus2 + Covhm21*mus1*mus2 + ...
                            %                                 Covh22*vars2 + Covh22*mus2^2 + muh22*muh22'*vars2 + R;
                            %                             CvZk = Poo2;
                            %                             %
                            %                             %                             %                             Pxo = Pp2*Hh2(:,:,1)*us(1,1) + Pp2*Hh2(:,:,2)*us(2,1);
                            %                             %                             %                             Pox = Pxo';
                            %                             Zkp = (muh21)*us(1,1) + (muh22)*us(2,1);
                            %
                            %                             I_CvZk=inv(CvZk);
                            %
                            %                             %----------------------------
                            %
                            %                             f_Hkp1mkp1Hk = 1/(sqrt(2*pi*det(CvZk)))*exp(-1/2*(Zk-Zkp)'*I_CvZk*(Zk-Zkp));
                            %
                            %                             %----------------------------
                            %
                            %                             %=================================================

                        end % End of i_Hk if

                        %_______________________

                        %-----------------------------II

                        f_Hkp1mkp1 = f_Hkp1mkp1 + f_Hkp1mkp1Hk*pHk_Hkp1(1,i_Hk);

                        %-----------------------------III

                        pmkp1_k = pmkp1_k + pmkp1_Hk_k*pHk_Hkp1(1,i_Hk);

                        %-----------------------------
                        %_______________________

                    end % end of i_Hk for

                    %______________

                    f_Hkp1(1,i_Hkp1) = f_Hkp1(1,i_Hkp1) + f_Hkp1mkp1*pmkp1_k;

                    %______________

                    %__________________

                    Xp2_Hkp1mkp1 = zeros(8,1);
                    %__________________

                    for i=1:2

                        Xp2_Hkp1mkp1 = Xp2_Hkp1mkp1 + Xp2_Hkp1mkp1Hk(:,i)*pHk_Hkp1(1,i);

                    end

                    %____________

                    Pp2_Hkp1mkp1 = zeros(8,8);
                    %____________

                    for i=1:2

                        Pp2_Hkp1mkp1 = Pp2_Hkp1mkp1 +...
                            pHk_Hkp1(1,i)*...
                            (Pp2_Hkp1mkp1Hk(:,:,i)+(Xp2_Hkp1mkp1Hk(:,i)-Xp2_Hkp1mkp1)*(Xp2_Hkp1mkp1Hk(:,i)-Xp2_Hkp1mkp1)');

                    end

                    %--------------------------- End of I

                    %===========================for JPM
                    %===========================

                    if tk>1

                        if Di(1,tk-1)==2

                            Xp2_Hkp1mkp1_JPM(:,i_mkp1) = Xp2_Hkp1mkp1;
                            Pp2_Hkp1mkp1_JPM(:,:,i_mkp1) = Pp2_Hkp1mkp1;

                            pmkp1_k_JPM(1,i_mkp1) = pmkp1_k;

                        end
                    end

                    %===========================
                    %===========================

                    %==============================Updating I

                    Xp2 = Xp2_Hkp1mkp1;
                    Pp2 = Pp2_Hkp1mkp1;

                    %----------------------S Estimation

                    for im=1:Ns

                        Mb2(im,1)=(DDa/sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd));

                        Mb2(im,2)=(DDa/sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd));

                    end

                    %----------------------

                    Hh21 = zeros(Ns,8);
                    Hh22 = zeros(Ns,8);

                    for im=1:Ns

                        Hh21(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                        Hh21(im,2)=0;
                        Hh21(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                        Hh21(im,4)=0;
                        Hh21(im,5:8)=zeros(1,4);

                        Hh22(im,1:4)=zeros(1,4);
                        Hh22(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                        Hh22(im,6)=0;
                        Hh22(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                        Hh22(im,8)=0;

                    end

                    Hhk2 = zeros(Ns,8);
                    for im=1:Ns

                        Hhk2(im,1)=-DDa*(Xp2(1,1)-xs(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                        Hhk2(im,2)=0;
                        Hhk2(im,3)=-DDa*(Xp2(3,1)-ys(im,1))/(sqrt((Xp2(1,1)-xs(im,1))^2+(Xp2(3,1)-ys(im,1))^2+dd))^3;
                        Hhk2(im,4)=0;
                        Hhk2(im,5)=-DDa*(Xp2(5,1)-xs(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                        Hhk2(im,6)=0;
                        Hhk2(im,7)=-DDa*(Xp2(7,1)-ys(im,1))/(sqrt((Xp2(5,1)-xs(im,1))^2+(Xp2(7,1)-ys(im,1))^2+dd))^3;
                        Hhk2(im,8)=0;

                    end

                    %----------------------------

                    Covh21 = [Hhk2(:,1:4) zeros(Ns,4)]*Pp2*[Hhk2(:,1:4) zeros(Ns,4)]';
                    Covh22 = [zeros(Ns,4) Hhk2(:,5:8)]*Pp2*[zeros(Ns,4) Hhk2(:,5:8)]';
                    Covhm12 = [Hhk2(:,1:4) zeros(Ns,4)]*Pp2*[zeros(Ns,4) Hhk2(:,5:8)]';
                    Covhm21 = [zeros(Ns,4) Hhk2(:,5:8)]*Pp2*[Hhk2(:,1:4) zeros(Ns,4)]';
                    vars1 = Cs(1,1);
                    vars2 = Cs(2,2);
                    mus1 = us(1,1);
                    mus2 = us(2,1);
                    muh21 = Mb2(:,1);
                    muh22 = Mb2(:,2);

                    Poo2 = Covh21*mus1^2' + Covh22*mus2^2 + Cs(1,1)*Covh21 + Cs(2,2)*Covh22 + Covhm21*mus1*mus2 + Covhm12*mus1*mus2 + ...
                        muh21*muh21'*vars1 + muh22*muh22'*vars2 + R;

                    Pxo2 = Pp2*([Hhk2(:,1:4) zeros(Ns,4)]*us(1,1))' + Pp2*([zeros(Ns,4) Hhk2(:,5:8)]*us(2,1))';
                    Pox2 = Pxo2';
                    Zkp2 = (muh21)*us(1,1) + (muh22)*us(2,1);

                    X2_Hkp1mkp1(:,i_mkp1) = Xp2 + Pxo2*inv(Poo2)*(Zk - Zkp2);
                    P2_Hkp1mkp1(:,:,i_mkp1) = Pp2 - Pxo2*inv(Poo2)*Pox2;

                    %========================================

                    %----------------------------

                    pmkp1_un(1,i_mkp1) = f_Hkp1mkp1*pmkp1_k;

                    %----------------------------

                    %========================================

                end % end of i_mkp1

                %================================II

                c1=sum(pmkp1_un);
                pmkp1 = pmkp1_un/c1;

                %================================

                %_______________

                X2_Hkp1 = zeros(8,1);
                %_______________

                for i=1:Nm

                    X2_Hkp1 = X2_Hkp1 + X2_Hkp1mkp1(:,i)*pmkp1(1,i);

                end

                %____________

                P2_Hkp1 = zeros(8,8);
                %____________

                for i=1:Nm

                    P2_Hkp1 = P2_Hkp1 +...
                        pmkp1(1,i)*...
                        (P2_Hkp1mkp1(:,:,i)+(X2_Hkp1mkp1(:,i)-X2_Hkp1)*(X2_Hkp1mkp1(:,i)-X2_Hkp1)');

                end

                %____________

                %================================

                %=================================== Start of JPM
                %===================================

                if tk>1

                    if Di(1,tk-1)==2
                        %_______________

                        Xp2_Hkp1_JPM = zeros(8,1);
                        %_______________

                        for i=1:Nm

                            Xp2_Hkp1_JPM = Xp2_Hkp1_JPM + Xp2_Hkp1mkp1_JPM(:,i)*pmkp1_k_JPM(1,i);

                        end

                        %____________

                        Pp2_Hkp1_JPM = zeros(8,8);
                        %____________

                        for i=1:Nm

                            Pp2_Hkp1_JPM = Pp2_Hkp1_JPM +...
                                pmkp1_k_JPM(1,i)*...
                                (Pp2_Hkp1mkp1_JPM(:,:,i)+(Xp2_Hkp1mkp1_JPM(:,i)-Xp2_Hkp1_JPM)*(Xp2_Hkp1mkp1_JPM(:,i)-Xp2_Hkp1_JPM)');

                        end

                        %------------------------------

                        for im=1:Ns

                            Mb2_JPM(im,1)=(DDa/sqrt((Xp2_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(3,1)-ys(im,1))^2+dd));

                            Mb2_JPM(im,2)=(DDa/sqrt((Xp2_Hkp1_JPM(5,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(7,1)-ys(im,1))^2+dd));

                        end

                        %___________________

                        for im=1:Ns

                            Hh2_JPM(im,1)=-DDa*(Xp2_Hkp1_JPM(1,1)-xs(im,1))/(sqrt((Xp2_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(3,1)-ys(im,1))^2+dd))^3;
                            Hh2_JPM(im,2)=0;
                            Hh2_JPM(im,3)=-DDa*(Xp2_Hkp1_JPM(3,1)-ys(im,1))/(sqrt((Xp2_Hkp1_JPM(1,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(3,1)-ys(im,1))^2+dd))^3;
                            Hh2_JPM(im,4)=0;
                            Hh2_JPM(im,5)=-DDa*(Xp2_Hkp1_JPM(5,1)-xs(im,1))/(sqrt((Xp2_Hkp1_JPM(5,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(7,1)-ys(im,1))^2+dd))^3;
                            Hh2_JPM(im,6)=0;
                            Hh2_JPM(im,7)=-DDa*(Xp2_Hkp1_JPM(7,1)-ys(im,1))/(sqrt((Xp2_Hkp1_JPM(5,1)-xs(im,1))^2+(Xp2_Hkp1_JPM(7,1)-ys(im,1))^2+dd))^3;
                            Hh2_JPM(im,8)=0;
                        end

                        %==========================================================================

                        Covh21_JPM = [Hh2_JPM(:,1:4) zeros(Ns,4)]*Pp2_Hkp1_JPM*[Hh2_JPM(:,1:4) zeros(Ns,4)]';
                        Covh22_JPM = [zeros(Ns,4) Hh2_JPM(:,5:8)]*Pp2_Hkp1_JPM*[zeros(Ns,4) Hh2_JPM(:,5:8)]';
                        Covhm12_JPM = [Hh2_JPM(:,1:4) zeros(Ns,4)]*Pp2_Hkp1_JPM*[zeros(Ns,4) Hh2_JPM(:,5:8)]';
                        Covhm21_JPM = [zeros(Ns,4) Hh2_JPM(:,5:8)]*Pp2_Hkp1_JPM*[Hh2_JPM(:,1:4) zeros(Ns,4)]';
                        vars1_JPM = Cs(1,1);
                        vars2_JPM = Cs(2,2);
                        mus1_JPM = us(1,1);
                        mus2_JPM = us(2,1);
                        muh21_JPM = Mb2_JPM(:,1);
                        muh22_JPM = Mb2_JPM(:,2);

                        Zkp1_k_Mean = muh21_JPM*mus1_JPM + muh22_JPM*mus2_JPM;

                        PZkp1_k_JPM = Covh21_JPM*vars1_JPM + Covh21_JPM*mus1_JPM^2 + muh21_JPM*muh21_JPM'*vars1_JPM +  ...
                            Covhm12_JPM*mus1_JPM*mus2_JPM + Covhm21_JPM*mus1_JPM*mus2_JPM + Covh22_JPM*vars2_JPM + Covh22_JPM*mus2_JPM^2 + muh22_JPM*muh22_JPM'*vars2_JPM + R;

                        %___________________
                        %------------------------------ Zkp1_k

                        Nj=10;

                        if det(PZkp1_k_JPM)<inf

                            Lamb=eig(PZkp1_k_JPM);
                            cmp=Lamb>0;
                            icmp=sum(cmp);

                            Dig = 0;
                            for ii=1:Ns

                                Dig = Dig + PZkp1_k_JPM(ii,ii);

                            end
                            IDig = imag(Dig);

                            if icmp==Ns && IDig==0

                                %----------------------
                                c_JPM(1,tk)=c_JPM(1,tk)+1;
                                %----------------------

                                D_PZkp1_k = chol(PZkp1_k_JPM);

                                Zkp1_k_JPM = diag(Zkp1_k_Mean)*ones(Ns,Nj) + D_PZkp1_k*randn(Ns,Nj);

                                %----------------------------- d calculation

                                Zk_JPM = diag(Zk)*ones(Ns,Nj);


                                d_JPM1 = sqrt(sum((Zk_JPM-Zkp1_k_JPM).^2,1));
                                d_JPM(1,tk) = (mean(d_JPM1));

% %                                 d_JPM1 = mean((Zk_JPM-Zkp1_k_JPM).^2,2);
% %                                 d_JPM(1,tk) = sqrt(mean(d_JPM1));

                                %------------------------------

                            end

                        end

                    end

                end

                %======================================================

                %================================
                %================================ End of JPM

            end % End of i_Hkp1 if

        end % End of i_Hkp1 for

        %============================= End of XHkp1 part

        %------------------------------ pHkp1  (j = 1 , 2)
        c1 = sum(f_Hkp1.*pHkp1_k);
        pHkp1 =  f_Hkp1.*pHkp1_k/c1;
        %------------------------------
        
        %+++++++++++++++++++++++++++++++++++++++
        
        pHkp1_t(:,tk) = pHkp1_t(:,tk) + pHkp1'; 
        
        TraceP1(1,tk) = TraceP1(1,tk) + trace(P1_Hkp1);
        
        TraceP2(1,tk) = TraceP2(1,tk) + trace(P2_Hkp1(1:4,1:4));
        TraceP2(2,tk) = TraceP2(2,tk) + trace(P2_Hkp1(5:8,5:8));
        
        %+++++++++++++++++++++++++++++++++++++++

        %============================

        for j=1:2

            for i=1:2

                %----------------mse part

                if j==1

                    if i==1

                        ep1(i,j) = trace(P1_Hkp1);
%                         ep2(i,j) = 0;

                    else

                        ep1(i,j) = Est_d;
%                         ep2(i,j) = Est_d;
                        
%                         if ep1(i,j)>19
%                             
%                             ep1(i,j) = 18;
%                             
%                         end

                    end

                else

                    if i==1

%                         ep1(i,j) = 0;%trace(P2xj_kp1)/j;
                        ep1(i,j) = Est_d;

                    else

                        ep1(i,j) = trace(P2_Hkp1)/j;
%                         ep2(i,j) = 0;
                        
                        
%                         if ep1(i,j)>19
%                             
%                             ep1(i,j) = 18;
%                             
%                         end

                    end

                end

                %------------ End of mse part

            end % End of i

        end % End of j

        e_kp1 = ep1;% + ep2;

        %================================== Decision

        %+++++++++++++++++++++++++

        Cz_kp1(1,1) = (c(1,1) + beta(1,1)*e_kp1(1,1))*pHkp1(1,1) + ...
                      (c(1,2) + beta(1,2)*e_kp1(1,2))*pHkp1(1,2);
                  
                  Cz_kp1(2,1) = (c(2,1) + beta(2,1)*e_kp1(2,1))*pHkp1(1,1) + ...
                                (c(2,2) + beta(2,2)*e_kp1(2,2))*pHkp1(1,2);
        %+++++++++++++++++++++++++
        
%         c_kp1 = alpha.*c + beta.*e_kp1;
% 
%         for i=1:2
% 
%             Cz_kp1(i,1) = sum(c_kp1(i,:).*pHkp1(1,:),2);
% 
%         end

        %-------------------

        [Vin,Inin]=min(Cz_kp1);
        ci1=Inin-1;

        %===================================

        if ci1==0  % D1

            %_____________Decision output

            Di(1,tk)=1;

            Xh1(:,tk) = X1_Hkp1(:,1);
            Ph1 = P1_Hkp1;

            %-----------------------------------S Estimation

            Mb1=zeros(Ns,1);

            for im=1:Ns

                Mb1(im,1)=(DDa/sqrt((Xh1(1,tk)-xs(im,1))^2+(Xh1(3,tk)-ys(im,1))^2+dd));

            end

            %==================================================

            Csh1=Cs(1,1)-Cs(1,1)*Mb1'*inv(R+Mb1*Cs(1,1)*Mb1')*Mb1*Cs(1,1);

            csh1 = us(1,1) + Csh1*Mb1'*inv(R)*(Zk-Mb1*us(1,1));

            S_hat1p(1,tk)=sa*csh1;

            %================================================== 2nd Step of
            %Filtering

            Hhh1 = zeros(Ns,4);

            for im=1:Ns

                Hhh1(im,1)=-DDa*(Xh1(1,tk)-xs(im,1))/(sqrt((Xh1(1,tk)-xs(im,1))^2+(Xh1(3,tk)-ys(im,1))^2+dd))^3;
                Hhh1(im,2)=0;
                Hhh1(im,3)=-DDa*(Xh1(3,tk)-ys(im,1))/(sqrt((Xh1(1,tk)-xs(im,1))^2+(Xh1(3,tk)-ys(im,1))^2+dd))^3;
                Hhh1(im,4)=0;

            end

            %===============

            Covh1 = Hhh1*Ph1*Hhh1';
            Poo1 = Covh1*(S_hat1p(1,tk)/sa)^2 + R;

            Pxo1 = Ph1*Hhh1'*S_hat1p(1,tk)/sa;
            Pox1 = Pxo1';
            muh1 = Mb1 ;
            Zkp1 = muh1*S_hat1p(1,tk)/sa;

            Xhh1(:,tk) = Xh1(:,tk) + Pxo1*inv(Poo1)*(Zk - Zkp1);
            Phh1 = Ph1 - Pxo1*inv(Poo1)*Pox1;

            %==================================================

        else  % D2

            %_____________Decision output

            Di(1,tk)=2;

            %_____________

            Xh2(:,tk) = X2_Hkp1;
            Ph2 = P2_Hkp1;

            %-----------------------
            %-----------------------------------S Estimation

            for im=1:Ns

                Mb2(im,1)=(DDa/sqrt((Xh2(1,tk)-xs(im,1))^2+(Xh2(3,tk)-ys(im,1))^2+dd));

                Mb2(im,2)=(DDa/sqrt((Xh2(5,tk)-xs(im,1))^2+(Xh2(7,tk)-ys(im,1))^2+dd));

            end

            %==================================================

            Csh2=Cs-Cs*Mb2'*inv(R+Mb2*Cs*Mb2')*Mb2*Cs;

            csh2=us+Csh2*Mb2'*inv(R)*(Zk-Mb2*us);

            S_hat2p(:,tk)=sa*csh2;

            %================================================== 2nd Step of
            %Filtering

            Hh2 = zeros(Ns,8,1);
            for im=1:Ns

                Hh2(im,1,1)=-DDa*(Xh2(1,tk)-xs(im,1))/(sqrt((Xh2(1,tk)-xs(im,1))^2+(Xh2(3,tk)-ys(im,1))^2+dd))^3;
                Hh2(im,2,1)=0;
                Hh2(im,3,1)=-DDa*(Xh2(3,tk)-ys(im,1))/(sqrt((Xh2(1,tk)-xs(im,1))^2+(Xh2(3,tk)-ys(im,1))^2+dd))^3;
                Hh2(im,4,1)=0;
                Hh2(im,5:8,1)=zeros(1,4);

                Hh2(im,1:4,2)=zeros(1,4);
                Hh2(im,5,2)=-DDa*(Xh2(5,tk)-xs(im,1))/(sqrt((Xh2(5,tk)-xs(im,1))^2+(Xh2(7,tk)-ys(im,1))^2+dd))^3;
                Hh2(im,6,2)=0;
                Hh2(im,7,2)=-DDa*(Xh2(7,tk)-ys(im,1))/(sqrt((Xh2(5,tk)-xs(im,1))^2+(Xh2(7,tk)-ys(im,1))^2+dd))^3;
                Hh2(im,8,2)=0;

            end

            %-----------------------------

            Covh21 = Hh2(:,:,1)*Ph2*Hh2(:,:,1)';
            Covh22 = Hh2(:,:,2)*Ph2*Hh2(:,:,2)';
            Covhm12 = Hh2(:,:,1)*Ph2*Hh2(:,:,2)';
            Covhm21 = Hh2(:,:,2)*Ph2*Hh2(:,:,1)';

            Poo2 = Covh21*(S_hat2p(1,tk)/sa)^2 + Covh22*(S_hat2p(2,tk)/sa)^2 + R + Covhm21*(S_hat2p(1,tk)/sa)*(S_hat2p(2,tk)/sa) + ...
                Covhm12*(S_hat2p(1,tk)/sa)*(S_hat2p(2,tk)/sa);

            Pxo2 = Ph2*Hh2(:,:,1)'*S_hat2p(1,tk)/sa + Ph2*Hh2(:,:,2)'*S_hat2p(2,tk)/sa;
            Pox2 = Pxo2';
            muh21 = Mb2(:,1);
            muh22 = Mb2(:,2);
            Zkp2 = muh21*S_hat2p(1,tk)/sa + muh22*S_hat2p(2,tk)/sa;

            Xhh2(:,tk) = Xh2(:,tk) + Pxo2*inv(Poo2)*(Zk - Zkp2);
            Phh2 = Ph2 - Pxo2*inv(Poo2)*Pox2;

            %==================================================

        end

        %==================

        Dei(1,tk) = abs(Sn(1,tk)-Di(1,tk));

        %==================

        %------------------------preparing for the next time sample

        X_mk_Hk_1 = X1_Hkp1mkp1;
        P_mk_Hk_1 = P1_Hkp1mkp1;

        X_mk_Hk_2 = X2_Hkp1mkp1;
        P_mk_Hk_2 = P2_Hkp1mkp1;

        pHk = pHkp1;

%         e_k=e_kp1;
%         c_k = c_kp1;
% 
%         Cz_k =Cz_kp1;

        %------------------------

        %==================================

    end % End of for tk

    %-------------------------

    D_JPM_MC(iter,:) = d_JPM;

    d_JPM=zeros(1,N);

    %-------------------------

    DDi(iter,:) = Di;

    DDei(iter,:) = Dei;

    %--------------------------------

    for i=1:N
        
        if i<10
            
            if Di(1,i)==1
                
                Row_D(1,i) = Row_D(1,i) + 1;
                
                SEs1p(1,i)= SEs1p(1,i) + abs(S_hat1p(1,i)-Ss1t(1,i));
                
                SEph1(1,i) = SEph1(1,i) + sqrt((X1t(1,i)-Xhh1(1,i))^2 + (X1t(3,i)-Xhh1(3,i))^2);
                SEvh1(1,i) = SEvh1(1,i) + sqrt((X1t(2,i)-Xhh1(2,i))^2 + (X1t(4,i)-Xhh1(4,i))^2);
                
%                 %.......................
%                 Sh1d1_1(1,i) = Sh1d1_1(1,i) + SEph1(1,i);
%                 %.......................
                
            else
                
                Row_D(1,i) = Row_D(1,i) + 1;
                
                %-----------------------
                
                SEph1_c(1,1) = sqrt((X1t(1,i)-Xhh2(1,i))^2 + (X1t(3,i)-Xhh2(3,i))^2);
                SEph1_c(2,1) = sqrt((X1t(1,i)-Xhh2(5,i))^2 + (X1t(3,i)-Xhh2(7,i))^2);
                
                [SEph1_r,ind_r] = min(SEph1_c);
                
               
                SEvh1_c(1,1) = sqrt((X1t(2,i)-Xhh2(2,i))^2 + (X1t(4,i)-Xhh2(4,i))^2);
                SEvh1_c(2,1) = sqrt((X1t(2,i)-Xhh2(6,i))^2 + (X1t(4,i)-Xhh2(8,i))^2);
                
                SEvh1_r = SEvh1_c(ind_r,1);
                
                SEs1p_c(1,1) = abs(S_hat2p(1,i)-Ss1t(1,i));
                SEs1p_c(2,1) = abs(S_hat2p(2,i)-Ss1t(1,i));
               
%                 %......................
%                 Sh1d2_1(1,i) = Sh1d2_1(1,i) + SEph1_c(1,1);
%                 Sh1d2_2(1,i) = Sh1d2_2(1,i) + SEph1_c(2,1);
%                 %......................
                
                SEs1_r = SEs1p_c(ind_r,1);
                
                %-----------------------
                
                SEs1p(1,i) = SEs1p(1,i) + SEs1_r;
                
                SEph1(1,i) = SEph1(1,i) + SEph1_r;
                SEvh1(1,i) = SEvh1(1,i) + SEvh1_r;
                
            end
            
        elseif i>9 && i<21
            
            if Di(1,i)==1
                
                %-------------------------
                
                SEph2_c(1,1) = sqrt((X1t(1,i)-Xhh1(1,i))^2 + (X1t(3,i)-Xhh1(3,i))^2);
                SEph2_c(2,1) = sqrt((X2t(1,i)-Xhh1(1,i))^2 + (X2t(3,i)-Xhh1(3,i))^2);
                
%                 %.........................
%                 Sh2d1_1(1,i) = Sh2d1_1(1,i) + SEph2_c(1,1);
%                 Sh2d1_2(1,i) = Sh2d1_2(1,i) + SEph2_c(2,1);
%                 %.........................
                
                [SEph_r,ind_r] = min(SEph2_c);
                
                Row_D(ind_r,i) = Row_D(ind_r,i) + 1;
                
                SEvh2_c(1,1) = sqrt((X1t(2,i)-Xhh1(2,i))^2 + (X1t(4,i)-Xhh1(4,i))^2);
                SEvh2_c(2,1) = sqrt((X2t(2,i)-Xhh1(2,i))^2 + (X2t(4,i)-Xhh1(4,i))^2);
                
                SEvh_r = SEvh2_c(ind_r,1);
                
                SEs2p_c(1,1) = abs(S_hat1p(1,i)-Ss1t(1,i));
                SEs2p_c(2,1) = abs(S_hat1p(1,i)-Ss2t(1,i));
                
                SEs_r = SEs2p_c(ind_r,1);
                
                %-----------------------
                
                if ind_r==1
                    
                    SEs1p(1,i) = SEs1p(1,i) + SEs_r;
                    
                    SEph1(1,i) = SEph1(1,i) + SEph_r;
                    SEvh1(1,i) = SEvh1(1,i) + SEvh_r;
                    
                else
                    
                    SEs2p(1,i) = SEs2p(1,i) + SEs_r;
                    
                    SEph2(1,i) = SEph2(1,i) + SEph_r;
                    SEvh2(1,i) = SEvh2(1,i) + SEvh_r;
                
                end
                
            else
                
                Row_D(1,i) = Row_D(1,i) + 1;
                Row_D(2,i) = Row_D(2,i) + 1;
                
                SEs1p(1,i) = SEs1p(1,i) + abs(S_hat2p(1,i)-Ss1t(1,i));
                SEs2p(1,i) = SEs2p(1,i) + abs(S_hat2p(2,i)-Ss2t(1,i));
                
                SEph1(1,i) = SEph1(1,i) + sqrt((X1t(1,i)-Xhh2(1,i))^2 + (X1t(3,i)-Xhh2(3,i))^2);
                SEvh1(1,i) = SEvh1(1,i) + sqrt((X1t(2,i)-Xhh2(2,i))^2 + (X1t(4,i)-Xhh2(4,i))^2);
                
                SEph2(1,i) = SEph2(1,i) + sqrt((X2t(1,i)-Xhh2(5,i))^2 + (X2t(3,i)-Xhh2(7,i))^2);
                SEvh2(1,i) = SEvh2(1,i) + sqrt((X2t(2,i)-Xhh2(6,i))^2 + (X2t(4,i)-Xhh2(8,i))^2);
                
%                  %.........................
%                 Sh2d2_1(1,i) = Sh2d2_1(1,i) + SEph1(1,i);
%                 Sh2d2_2(1,i) = Sh2d2_2(1,i) + SEph2(2,i);
%                 %.........................
                
            end
            
        elseif i>20
            
            if Di(1,i)==1
                
                Row_D(1,i) = Row_D(1,i) + 1;
                
                SEs1p(1,i)= SEs1p(1,i) + abs(S_hat1p(1,i)-Ss1t(1,i));
                
                SEph1(1,i) = SEph1(1,i) + sqrt((X1t(1,i)-Xhh1(1,i))^2 + (X1t(3,i)-Xhh1(3,i))^2);
                SEvh1(1,i) = SEvh1(1,i) + sqrt((X1t(2,i)-Xhh1(2,i))^2 + (X1t(4,i)-Xhh1(4,i))^2);
                
%                  %.........................
%                 Sh1d1_1(1,i) = Sh1d1_1(1,i) + SEph1(1,i);
%                
%                 %.........................
                
            else
                
                Row_D(1,i) = Row_D(1,i) + 1;
                
                %-----------------------
                
                SEph1_c(1,1) = sqrt((X1t(1,i)-Xhh2(1,i))^2 + (X1t(3,i)-Xhh2(3,i))^2);
                SEph1_c(2,1) = sqrt((X1t(1,i)-Xhh2(5,i))^2 + (X1t(3,i)-Xhh2(7,i))^2);
                
                [SEph1_r,ind_r] = min(SEph1_c);
                
                SEvh1_c(1,1) = sqrt((X1t(2,i)-Xhh2(2,i))^2 + (X1t(4,i)-Xhh2(4,i))^2);
                SEvh1_c(2,1) = sqrt((X1t(2,i)-Xhh2(6,i))^2 + (X1t(4,i)-Xhh2(8,i))^2);
                
                SEvh1_r = SEvh1_c(ind_r,1);
                
                SEs1p_c(1,1) = abs(S_hat2p(1,i)-Ss1t(1,i));
                SEs1p_c(2,1) = abs(S_hat2p(2,i)-Ss1t(1,i));
                
                SEs1_r = SEs1p_c(ind_r,1);
                
%                  %.........................
%                 Sh1d2_1(1,i) = Sh1d2_1(1,i) + SEph1(1,i);
%                 Sh1d2_2(1,i) = Sh1d2_2(1,i) + SEph2(2,i);
%                 %.........................
                %-----------------------
                
                SEs1p(1,i) = SEs1p(1,i) + SEs1_r;
                
                %.....................
                
                if ind_r==1
                    
                    SEs1p_nr(1,i) = SEs1p_nr(1,i) + SEs1p_c(1,1);
                
                else
                    
                    SEs1p_nr(1,i) = SEs1p_nr(1,i) + SEs1p_c(2,1);
                end
                
                %...
                
                if ind_r==1
                    
                    SEph1_nr(1,i) = SEph1_nr(1,i) + SEph1_c(1,1);
                
                else
                    
                    SEph1_nr(1,i) = SEph1_nr(1,i) + SEph1_c(2,1);
                end
                %.....................
                
                SEph1(1,i) = SEph1(1,i) + SEph1_r;
                SEvh1(1,i) = SEvh1(1,i) + SEvh1_r;
                
            end
            
        end
        
    end

    %--------------------------------
    %==========================================================

end % End of iteration

TperLoop = toc/iteration
%--------------------------

for i=1:N

    D_JPM(1,i) = sum(D_JPM_MC(:,i),1)/c_JPM(1,i);

end

%--------------------------

DFi=mean(DDi,1);

DFei=mean(DDei,1);

%+++++++++++++++++++++++++++++++++++

AEEs1p =zeros(1,N);
AEEph1 =zeros(1,N);
AEEvh1 =zeros(1,N);

AEEs1p_nr =zeros(1,N);
AEEph1_nr =zeros(1,N);


for i=1:N
    
    if Row_D(1,i)>0
        
        AEEs1p(1,i) = SEs1p(1,i)/Row_D(1,i);
        
        AEEph1(1,i) = SEph1(1,i)/Row_D(1,i);
        
        AEEvh1(1,i) = SEvh1(1,i)/Row_D(1,i);
    
        
        AEEs1p_nr = SEs1p_nr(1,i)/Row_D(1,i);
        
        AEEph1_nr = SEph1_nr(1,i)/Row_D(1,i);
        
    else
   
        AEEs1p(1,i) = SEs1p(1,i)/Row_D(1,i);
        
        AEEph1(1,i) = SEph1(1,i)/Row_D(1,i);
        
        AEEvh1(1,i) = SEvh1(1,i)/Row_D(1,i);
   
        
        AEEs1p_nr = SEs1p_nr(1,i)/Row_D(1,i);
        
        AEEph1_nr = SEph1_nr(1,i)/Row_D(1,i);
        
    end
    
end

AEEs2p =zeros(1,N);
AEEph2 =zeros(1,N);
AEEvh2 =zeros(1,N);

for i=10:20
    
    if Row_D(2,i)>0
        
        AEEs2p(1,i) = SEs2p(1,i)/Row_D(2,i);
        
        AEEph2(1,i) = SEph2(1,i)/Row_D(2,i);
        
        AEEvh2(1,i) = SEvh2(1,i)/Row_D(2,i);
        
    else
   
        AEEs2p(1,i) = 0;
        
        AEEph2(1,i) = 0;
        
        AEEvh2(1,i) = 0;
   
    end
        
    
end

%+++++++++++++++++++++++++++++++++++

%-------------Signals
% 
% RMSEs1p= (mean(SEs1p,1));
% 
% RMSEs2p= (mean(SEs2p,1));
% 
% % %-------------
% 
% RMSEp1= (mean(SEp1,1));
% RMSEv1= (mean(SEv1,1));
% 
% RMSEph1= (mean(SEph1,1));
% RMSEvh1= (mean(SEvh1,1));
% 
% RMSEp2= (mean(SEp2,1));
% RMSEv2= (mean(SEv2,1));
% 
% RMSEph2= (mean(SEph2,1));
% RMSEvh2= (mean(SEvh2,1));
% 
% 
% %+++++++++++++++++++++++++++++++++
% 
% RMSEs2p(1,1:9) = zeros(1,9);
% RMSEph2(1,1:9) = zeros(1,9);
% RMSEvh2(1,1:9) = zeros(1,9);
% 
% RMSEs2p(1,21:32) = zeros(1,12);
% RMSEph2(1,21:32) = zeros(1,12);
% RMSEvh2(1,21:32) = zeros(1,12);

%+++++++++++++++++++++++++++++++++
%------------------------

figure(1)
hold on
plot(X1t(1,:),X1t(3,:))
hold on
plot(X2t(1,:),X2t(3,:),'r')
hold on
plot(Xhh1(1,:),Xhh1(3,:),'o')
hold on
plot(Xhh2(1,:),Xhh2(3,:),'.r')
hold on
plot(Xhh2(5,:),Xhh2(7,:),'D')
axis([0 100 0 100])

figure(2)
hold on
plot(AEEph1)
hold on
plot(AEEph2,'r--')
% hold on
% plot(AEEph1_nr,'k-.')

figure(3)
hold on
plot(AEEvh1)
hold on
plot(AEEvh2,'r--')

figure(4)
hold on
plot(AEEs1p)
hold on
plot(AEEs2p,'r--')
% hold on
% plot(AEEs1p_nr,'k-.')

figure(5)
hold on
plot(D_JPM)

% figure(6)
% hold on
% plot(c_JPM)

figure(7)
plot(DFi)

figure(8)
plot(DFei)

% figure(9)
% hold on 
% plot(Row_D(1,:))
% hold on
% plot(Row_D(2,:),'r')


%--------------------------

% ATraceP2 = mean(TraceP2,1);

%--------------------------

% figure(10)
% hold on
% plot(TraceP1/iteration)
% hold on
% plot(TraceP2(1,:)/iteration,'r')
% hold on
% plot(TraceP2(2,:)/iteration,'k')
% hold on
% plot(ATraceP2/iteration,'g')


figure(11)
hold on
plot(pHkp1_t(1,:)/iteration)
hold on
plot(pHkp1_t(2,:)/iteration,'r')


