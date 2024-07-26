clc;clear
syms X;
%INPUTS
%NELEM_array = [1 10 20 40 80 100]; %no. of elements
NELEM_array = [1 2 5 10 100]; %no. of elements
XLEN  = 1; %length of bar
IPVAL = 3; %Approximation order
AE    = [1 0] ; %Linear Variation of AE 1+2*x
C     = [0 0 0]; %quadratic variation of distributed springs 1+2*x+x^2
T     = [0 1 0 0]; %quadratic variation of distributes tractions 0+1*x+0*x^2+0*sin(pi*x/L)

%Type of BC(Dirichlet-1 or Neumann-2 or Robin-3)
BCTYPE = [1 2];% A B
BC = [0 0];

% BCTYPE = [2 3];
% BC = [1/pi 10 0];% A B C

%Type of Approximating Function(Lagrange-1 or Hierarchic-2)
TYPE = 1;

X_p=0:1/(100):XLEN;
folder_path='/Users/nikhil/Desktop/FEM/A1/Q2/';

U_fem_p=zeros(length(NELEM_array),length(X_p));
DU_fem_p=zeros(length(NELEM_array),length(X_p));
error_fem_p=zeros(length(NELEM_array),length(X_p));
Energy1=zeros(length(NELEM_array),1);
Energy2=zeros(length(NELEM_array),1);


 %Exact Solution
AE_ = AE(1)+AE(2)*X;
C_ = C(1)+C(2)*X+C(3)*X^2;      % stiffness coefficient
T_ = T(1)+T(2)*X+T(3)*X^2+T(4)*sin(pi*X/XLEN);      % distributed load
a= 0;      % left end of the bar
b= XLEN;      % right end of the bar
XLEN=b-a;

syms Uexact(X)
diff_eq = diff(AE_*Uexact,X,2) - C_*Uexact + T_ == 0;
if(BCTYPE(1)==1)
    COND1 = subs(Uexact,X,a) == BC(1);
elseif(BCTYPE(1)==2)
    COND1 = subs(diff(Uexact,X,1),X,a) == BC(1);

end
if(BCTYPE(2)==1)
    COND2 = subs(Uexact,X,b) == BC(2);
elseif(BCTYPE(2)==2)
    COND2 = subs(diff(Uexact,X,1),X,b) == BC(2);
elseif(BCTYPE(2)==3)
    COND2 = subs(diff(Uexact,X,1),X,b) == BC(2)*(BC(3)-subs(Uexact,X,b));

end
Uexact = dsolve(diff_eq, COND1, COND2);
Uexact
Uexact_eval = double(subs(Uexact,X,X_p));
DUexact_eval = double(subs(diff(Uexact),X,X_p));


for P=1:length(NELEM_array)
    NELEM=NELEM_array(P);
    %PRE-PROCESSING
    
    NDOF=IPVAL+1;
    
    %Mesh Generation
    h_I=XLEN/NELEM; %Uniform Mesh
    X_node=zeros(NELEM+1,1);
    for i=1:NELEM+1
        X_node(i)=(i-1)*h_I;
    end
   
    %Fix the order of Integration
    NINT=int32(((2+4+4)+1)/2);
    
    %Load Integration Points
    X_i=[-0.93246951 -0.66120939 -0.23861919 0.23861919 0.66120939 0.93246951];
    W_i=[0.17132449 0.36076157 0.46791393 0.46791393 0.36076157 0.17132449];
    
    %load shape functions
    [N,DN]=Shape_Function(TYPE,NINT,X_i,IPVAL);
    [EF,EK]= Element_Cal(1,X_node,X_i,W_i,NINT,NDOF,AE,C,T,N,DN,XLEN);
    [GLOBF,GLOBK]=Local_Global(IPVAL,NELEM,NDOF,X_node,X_i,W_i,NINT,AE,C,T,N,DN,XLEN);
    [GLOBF, GLOBK]=Boundary_Conditions(GLOBF,GLOBK,BCTYPE,BC,NELEM,IPVAL);
    U=GLOBK\GLOBF;
    
    %FEM Solution
    
    [U_FEM,DU_FEM]=solution(TYPE,NELEM,IPVAL,X_node,h_I,U);
    [U_FEM_p,DU_FEM_p]=plotting(NELEM,IPVAL,U_FEM,DU_FEM,X_p,X_node);
    
    %Error
    for i=1:length(X_p)
        error(i)=(Uexact_eval(i)-U_FEM_p(i));
    end
    
    f1=@(x)(AE(1)+AE(2)*x)*subs(diff(Uexact),X,x)^2+(C(1)+C(2)*x+C(3)*x^2)*subs(Uexact,X,x)^2;
    Energy_exact=double(Integrate(f1,0,XLEN,10));
    Energy_FEM=0;
    for i=1:NELEM
        f2=@(x)(AE(1)+AE(2)*x)*subs(DU_FEM(i),X,x)^2+(C(1)+C(2)*x+C(3)*x^2)*subs(Uexact,X,x)^2;
        Energy_FEM=Energy_FEM+double(Integrate(f2,X_node(i),X_node(i+1),10));
    end
    Energy1(P)=Energy_exact;
    Energy2(P)=Energy_FEM;

    for m=1:length(X_p)
        U_fem_p(P,m)=U_FEM_p(m);
        DU_fem_p(P,m)=DU_FEM_p(m);
        error_fem_p(P,m)=error(m);
    end
    LP(P)=0;
    SEE(P)=0;
    for p=1:NELEM
        f3=@(x) (subs(Uexact,X,x)-subs(U_FEM(p),X,x))^2;
        f4=@(x) subs(Uexact,X,x)-subs(U_FEM(p),X,x);
        LP(P)=LP(P)+double(Integrate(f3,X_node(p),X_node(p+1),10));
        SEE(P)=SEE(P)+double(Integrate(f4,X_node(p),X_node(p+1),10));
    end
    LP(P)=LP(P)^0.5;

end

figure;
plot(X_p,Uexact_eval,'LineWidth',2);
hold on;
plot(X_p,U_fem_p(1,:),'LineWidth',2);
plot(X_p,U_fem_p(2,:),'LineWidth',2);
plot(X_p,U_fem_p(3,:),'LineWidth',2);
plot(X_p,U_fem_p(4,:),'LineWidth',2);
plot(X_p,U_fem_p(5,:),'LineWidth',2);
%plot(X_p,U_fem_p(6,:),'LineWidth',2);

xlabel('x');
ylabel('u');
title('Finite Element Solution');
text1=['u_{Exact}'];
text2=['u_{FEM} for ',num2str(NELEM_array(1)),' elements'];
text3=['u_{FEM} for ',num2str(NELEM_array(2)),' elements'];
text4=['u_{FEM} for ',num2str(NELEM_array(3)),' elements'];
text5=['u_{FEM} for ',num2str(NELEM_array(4)),' elements'];
text6=['u_{FEM} for ',num2str(NELEM_array(5)),' elements'];
%text7=['u_{FEM} for ',num2str(NELEM_array(6)),' elements'];
legend(text1,text2,text3,text4,text5,text6)%,text7)
grid on;
hold off;
file_name=['U_plot',num2str(IPVAL),'.png'];
saveas(gcf,fullfile(folder_path,file_name));

figure;
plot(X_p,DUexact_eval,'LineWidth',2);
hold on;
plot(X_p,DU_fem_p(1,:),'LineWidth',2);
plot(X_p,DU_fem_p(2,:),'LineWidth',2);
plot(X_p,DU_fem_p(3,:),'LineWidth',2);
plot(X_p,DU_fem_p(4,:),'LineWidth',2);
plot(X_p,DU_fem_p(5,:),'LineWidth',2);
%plot(X_p,DU_fem_p(6,:),'LineWidth',2);
xlabel('x');
ylabel('du/dx');
title('Finite Element Solution');
text1=['du/dx_{Exact}'];
text2=[num2str(NELEM_array(1)),' elements'];
text3=[num2str(NELEM_array(2)),' elements'];
text4=[num2str(NELEM_array(3)),' elements'];
text5=[num2str(NELEM_array(4)),' elements'];
text6=[num2str(NELEM_array(5)),' elements'];
%text7=[num2str(NELEM_array(6)),' elements'];
legend(text1,text2,text3,text4,text5,text6)%,text7)
grid on;
hold off;
file_name=['DU_plot',num2str(IPVAL),'.png'];
saveas(gcf,fullfile(folder_path,file_name));

figure;
plot(X_p,error_fem_p(1,:),'LineWidth',2);
hold on;
plot(X_p,error_fem_p(2,:),'LineWidth',2);
plot(X_p,error_fem_p(3,:),'LineWidth',2);
plot(X_p,error_fem_p(4,:),'LineWidth',2);
plot(X_p,error_fem_p(5,:),'LineWidth',2);
%plot(X_p,error_fem_p(6,:),'LineWidth',2);

xlabel('x');
ylabel('Error');
title('Finite Element Error');
text2=[num2str(NELEM_array(1)),' elements'];
text3=[num2str(NELEM_array(2)),' elements'];
text4=[num2str(NELEM_array(3)),' elements'];
text5=[num2str(NELEM_array(4)),' elements'];
text6=[num2str(NELEM_array(5)),' elements'];
%text7=[num2str(NELEM_array(6)),' elements'];
legend(text2,text3,text4,text5,text6)%,text7);
grid on;
hold off;
file_name=['Error_plot',num2str(IPVAL),'.png'];
saveas(gcf,fullfile(folder_path,file_name));

figure;
plot(NELEM_array,Energy1,'LineWidth',2);
hold on;
plot(NELEM_array,Energy2,'LineWidth',2);
xlabel('NELEM');
ylabel('Strain Energy');
title('Strain Energy vs NELEM');
legend('Exact','FEM')
hold off;
file_name=['Energy_Plot','.png'];
saveas(gcf,fullfile(folder_path,file_name));

figure;
plot(NELEM_array,SEE,'LineWidth',2);
xlabel('NELEM');
ylabel('Strain Energy of Error');
title('Strain Energy of Error vs NELEM');
file_name=['SOE_Plot',num2str(IPVAL),'.png'];
saveas(gcf,fullfile(folder_path,file_name));

figure;
plot(log(NELEM_array),log(LP),'LineWidth',2);
xlabel('log(NELEM)');
ylabel('log (relative error)');
title('log (relative error) vs log (NELEM)');
file_name=['LP_Plot',num2str(IPVAL),'.png'];
saveas(gcf,fullfile(folder_path,file_name));



function [F,DF]=solution(TYPE,NELEM,IPVAL,X_node,h_I,U1)
    for i=1:NELEM
        syms X;
        xk=X_node(i);
        xk1=X_node(i+1);
        a=(xk1+xk)/2;
        b=(xk1-xk)/2;
        X_=(X-a)/b;
        if(TYPE==1)
            %Linear Shape Functions
            if(IPVAL== 1)
                N(1) = 1/2 - X_/2;
                N(2) = X_/2 + 1/2;
                DN(1)= -1/2;
                DN(2)= 1/2;
                
                F(i)=(N(1)*U1(i)+N(2)*U1(i+1));
                DF(i)=(DN(1)*U1(i)+DN(2)*U1(i+1))/b;
                % F(i)=(N(1)*U1(IPVAL*(i-1)+1)+N(2)*U1(IPVAL*(i-1)+2))*XJAC;
                % DF(i)=(DN(1)*U1(IPVAL*(i-1)+1)+N(2)*U1(IPVAL*(i-1)+2))/XJAC;
        
            elseif(IPVAL==2)
                %Quadratic Shape Functions
                N(1)= (X_*(X_-1))/2;
                N(2)=-(X_-1)*(X_+1);
                N(3)= (X_*(X_+1))/2;
                DN(1)= X_ - 1/2;
                DN(2)=-2*X_;
                DN(3)= X_ + 1/2;

                F(i)=N(1)*U1(2*i-1)+N(2)*U1(2*i)+N(3)*U1(2*i+1);
                DF(i)=(DN(1)*U1(2*i-1)+DN(2)*U1(2*i)+DN(3)*U1(2*i+1))/b;


            elseif(IPVAL==3)
                N(1)=-(9*(X_ - 1)*(X_ - 1/3)*(X_ + 1/3))/16;
                N(2)=(27*(X_ - 1)*(X_ + 1)*(X_ - 1/3))/16;
                N(3)=-(27*(X_ - 1)*(X_ + 1)*(X_ + 1/3))/16;
                N(4)=(9*(X_ + 1)*(X_ - 1/3)*(X_ + 1/3))/16;
                DN(1)=(9*X_)/8 - (27*X_^2)/16 + 1/16;
                DN(2)=(81*X_^2)/16 - (9*X_)/8 - 27/16;
                DN(3)=27/16 - (81*X_^2)/16 - (9*X_)/8;
                DN(4)=(27*X_^2)/16 + (9*X_)/8 - 1/16;

                F(i)=N(1)*U1(3*i-2)+N(2)*U1(3*i-1)+N(3)*U1(3*i)+N(4)*U1(3*i+1);
                DF(i)=(DN(1)*U1(3*i-2)+DN(2)*U1(3*i-1)+DN(3)*U1(3*i)+DN(4)*U1(3*i+1))/b;

            elseif(IPVAL==4)
                N(1)= (2*X_*(X_ - 1)*(X_ - 1/2)*(X_ + 1/2))/3;
                N(2)=-(8*X_*(X_- 1)*(X_ + 1)*(X_ - 1/2))/3;
                N(3)= 4*X_^4 - 5*X_^2 + 1;
                N(4)=-(8*X_*(X_ - 1)*(X_ + 1)*(X_ + 1/2))/3;
                N(5)= (2*X_*(X_ + 1)*(X_ - 1/2)*(X_ + 1/2))/3;
                DN(1)= (8*X_^3)/3 - 2*X_^2 - X_/3 + 1/6;
                DN(2)=-(32*X_^3)/3 + 4*X_^2 + (16*X_)/3 - 4/3;
                DN(3)= 2*X_*(8*X_^2 - 5);
                DN(4)=-(32*X_^3)/3 - 4*X_^2 + (16*X_)/3 + 4/3;
                DN(5)= (8*X_^3)/3 + 2*X_^2 - X_/3 - 1/6;

                F(i)=N(1)*U1(4*i-3)+N(2)*U1(4*i-2)+N(3)*U1(4*i-1)+N(4)*U1(4*i)+N(5)*U1(4*i+1);
                DF(i)=(DN(1)*U1(4*i-3)+DN(2)*U1(4*i-2)+DN(3)*U1(4*i-1)+DN(4)*U1(4*i)+DN(5)*U1(4*i+1))/b;

            end

        elseif(TYPE==2)%Hierarchic Shape Functions
            if(IPVAL== 1)
                N(1) = (1-X_)/2;
                N(2) = (1+X_)/2;
                DN(1)= -1/2;
                DN(2)= 1/2;
                F(i)=(N(1)*U1(i)+N(2)*U1(i+1));
                DF(i)=(DN(1)*U1(i)+DN(2)*U1(i+1))/b;
        
            elseif(IPVAL==2)
                %Quadratic Shape Functions
                N(1) = (1-X_)/2;
                N(3) = (1+X_)/2;
                N(2) = 3/(2*sqrt(6))*(X_^2-1);
                DN(1) = -1/2;
                DN(3) = 1/2;
                DN(2) = 3/sqrt(6);
                F(i)=N(1)*U1(2*i-1)+N(2)*U1(2*i)+N(3)*U1(2*i+1);
                DF(i)=(DN(1)*U1(2*i-1)+DN(2)*U1(2*i)+DN(3)*U1(2*i+1))/b;

            elseif(IPVAL==3)
                %Quadratic Shape Functions
                N(1) = (1-X_)/2;
                N(4) = (1+X_)/2;
                N(2) = 3/(2*sqrt(6))*(X_^2-1);
                N(3) = 5/(2*sqrt(10))*(X_*(X_^2-1));
                DN(1) = -1/2;
                DN(4) = 1/2;
                DN(2) = 3/sqrt(6);
                DN(3) = 15/(2*sqrt(10))*X_^2-5/(2*sqrt(10));
                F(i)=N(1)*U1(3*i-2)+N(2)*U1(3*i-1)+N(3)*U1(3*i)+N(4)*U1(3*i+1);
                DF(i)=(DN(1)*U1(3*i-2)+DN(2)*U1(3*i-1)+DN(3)*U1(3*i)+DN(4)*U1(3*i+1))/b;

            elseif(IPVAL==4)
                %Quadratic Shape Functions
                N(1) = (1-X_)/2;
                N(5) = (1+X_)/2;
                N(2) = 3/(2*sqrt(6))*(X_^2-1);
                N(3) = 5/(2*sqrt(10))*(X_*(X_^2-1));
                N(4) = 7/(8*sqrt(14))*(5*X_^4-6*X_^2+1);
                DN(1) = -1/2;
                DN(5) = 1/2;
                DN(2) = 3/sqrt(6);
                DN(3) = 15/(2*sqrt(10))*X_^2-5/(2*sqrt(10));
                DN(4) = 35/(2*sqrt(14))*X_^3-21/(2*sqrt(14))*X_;
                
                F(i)=N(1)*U1(4*i-3)+N(2)*U1(4*i-2)+N(3)*U1(4*i-1)+N(4)*U1(4*i)+N(5)*U1(4*i+1);
                DF(i)=(DN(1)*U1(4*i-3)+DN(2)*U1(4*i-2)+DN(3)*U1(4*i-1)+DN(4)*U1(4*i)+DN(5)*U1(4*i+1))/b;
            end
        end
    end
end

%Integration using Gaussian Quadrature
function I=Integrate(f,A,B,n)

    %Gaussian variables
    x1=-0.932469514;
    x2=-0.661209386;
    x3=-0.2386191860;
    w1=0.171324492;
    w2=0.360761573;
    w3=0.467913935;

    Integral=0;
    diff=(B-A)/n;
    %changing inegration
    for i=1:n
        t=@(x) (diff)/2*x+(A+i*diff-diff/2);
        Integral = Integral + (f(t(x1))*w1 + f(t(x2))*w2 + f(t(x3))*w3 + f(t(-1*x3))*w3 + f(t(-1*x2))*w2 + f(t(-1*x1))*w1 ) * (diff)/2;
    end 
    I=Integral;
end


%Plotting
function [U_p,DU_p]=plotting(NELEM,IPVAL,U_FEM,DU_FEM,X_p,X_node)
    U_p=zeros(length(X_p),1);
    DU_p=zeros(length(X_p),1);
    syms X;
    for i=1:length(X_p)
        for j=1:NELEM
            if(X_p(i)==X_node(NELEM+1))
                U_p(i)=subs(U_FEM(j),X,X_p(i));
                DU_p(i)=subs(DU_FEM(j),X,X_p(i));

            elseif(X_node(j)<=X_p(i) && X_p(i)<X_node(j+1))
                U_p(i)=subs(U_FEM(j),X,X_p(i));
                DU_p(i)=subs(DU_FEM(j),X,X_p(i));
            end
        end
    end
end


%shape function
function [N,DN] = Shape_Function(TYPE,NINT,X_i,IPVAL)
    syms X;
    for i = 1:NINT
        if(TYPE==1)
            %Linear Shape Functions
            if(IPVAL== 1)
                N(1,i) = 1/2 - X_i(i)/2;
                N(2,i) = X_i(i)/2 + 1/2;
                DN(1,i)= -1/2;
                DN(2,i)= 1/2;

            elseif(IPVAL==2)
                %Quadratic Shape Functions
                N(1,i)= (X_i(i)*(X_i(i)-1))/2;
                N(2,i)=-(X_i(i)-1)*(X_i(i)+1);
                N(3,i)= (X_i(i)*(X_i(i)+1))/2;
                DN(1,i)= X_i(i) - 1/2;
                DN(2,i)=-2*X_i(i);
                DN(3,i)= X_i(i) + 1/2;

            elseif(IPVAL==3)
                N(1,i)=-(9*(X_i(i) - 1)*(X_i(i) - 1/3)*(X_i(i) + 1/3))/16;
                N(2,i)=(27*(X_i(i) - 1)*(X_i(i) + 1)*(X_i(i) - 1/3))/16;
                N(3,i)=-(27*(X_i(i) - 1)*(X_i(i) + 1)*(X_i(i) + 1/3))/16;
                N(4,i)=(9*(X_i(i) + 1)*(X_i(i) - 1/3)*(X_i(i) + 1/3))/16;
                DN(1,i)=(9*X_i(i))/8 - (27*X_i(i)^2)/16 + 1/16;
                DN(2,i)=(81*X_i(i)^2)/16 - (9*X_i(i))/8 - 27/16;
                DN(3,i)=27/16 - (81*X_i(i)^2)/16 - (9*X_i(i))/8;
                DN(4,i)=(27*X_i(i)^2)/16 + (9*X_i(i))/8 - 1/16;

            elseif(IPVAL==4)
                N(1,i)= (2*X_i(i)*(X_i(i) - 1)*(X_i(i) - 1/2)*(X_i(i) + 1/2))/3;
                N(2,i)=-(8*X_i(i)*(X_i(i) - 1)*(X_i(i) + 1)*(X_i(i) - 1/2))/3;
                N(3,i)= 4*X_i(i)^4 - 5*X_i(i)^2 + 1;
                N(4,i)=-(8*X_i(i)*(X_i(i) - 1)*(X_i(i) + 1)*(X_i(i) + 1/2))/3;
                N(5,i)= (2*X_i(i)*(X_i(i) + 1)*(X_i(i) - 1/2)*(X_i(i) + 1/2))/3;
                DN(1,i)= (8*X_i(i)^3)/3 - 2*X_i(i)^2 - X_i(i)/3 + 1/6;
                DN(2,i)=-(32*X_i(i)^3)/3 + 4*X_i(i)^2 + (16*X_i(i))/3 - 4/3;
                DN(3,i)= 2*X_i(i)*(8*X_i(i)^2 - 5);
                DN(4,i)=-(32*X_i(i)^3)/3 - 4*X_i(i)^2 + (16*X_i(i))/3 + 4/3;
                DN(5,i)= (8*X_i(i)^3)/3 + 2*X_i(i)^2 - X_i(i)/3 - 1/6;
            end

        elseif(TYPE==2)%Hierarchic Shape Functions
            if(IPVAL== 1)
                N(1,i) = (1-X_i(i))/2;
                N(2,i) = (1+X_i(i))/2;
                DN(1,i)= -1/2;
                DN(2,i)= 1/2;
        
            elseif(IPVAL==2)
                %Quadratic Shape Functions
                N(1,i) = (1-X_i(i))/2;
                N(3,i) = (1+X_i(i))/2;
                N(2,i) = 3/(2*sqrt(6))*(X_i(i)^2-1);
                DN(1,i) = -1/2;
                DN(3,i) = 1/2;
                DN(2,i) = 3/sqrt(6);

            elseif(IPVAL==3)
                %Quadratic Shape Functions
                N(1,i) = (1-X_i(i))/2;
                N(4,i) = (1+X_i(i))/2;
                N(2,i) = 3/(2*sqrt(6))*(X_i(i)^2-1);
                N(3,i) = 5/(2*sqrt(10))*(X_i(i)*(X_i(i)^2-1));
                DN(1,i) = -1/2;
                DN(4,i) = 1/2;
                DN(2,i) = 3/sqrt(6);
                DN(3,i) = 15/(2*sqrt(10))*X_i(i)^2-5/(2*sqrt(10));

            elseif(IPVAL==4)
                %Quadratic Shape Functions
                N(1,i) = (1-X_i(i))/2;
                N(5,i) = (1+X_i(i))/2;
                N(2,i) = 3/(2*sqrt(6))*(X_i(i)^2-1);
                N(3,i) = 5/(2*sqrt(10))*(X_i(i)*(X_i(i)^2-1));
                N(4,i) = 7/(8*sqrt(14))*(5*X_i(i)^4-6*X_i(i)^2+1);
                DN(1,i) = -1/2;
                DN(5,i) = 1/2;
                DN(2,i) = 3/sqrt(6);
                DN(3,i) = 15/(2*sqrt(10))*X_i(i)^2-5/(2*sqrt(10));
                DN(4,i) = 35/(2*sqrt(14))*X_i(i)^3-21/(2*sqrt(14))*X_i(i);
            end
        end
    end
end

function [EF,EK]= Element_Cal(IEL,X,X_i,W_i,NINT,NDOF,AE,C,T,N,DN,XLEN)
    %Initialize
    EF=zeros(NDOF,1);
    EK=zeros(NDOF,NDOF);

    %Jacobian
    X1=X(IEL);
    X2=X(IEL+1);
    XJAC=(X2-X1)/2;

    for l=1:NINT
        XXi = X1 * (1-X_i(l))/2+ X2 * (1+X_i(l))/2;
        AE_XXi = AE(1)+AE(2)*XXi;
        T_XXi = T(1)+T(2)*XXi+T(3)*XXi^2+T(4)*sin(pi*XXi/XLEN);
        C_XXi = C(1)+C(2)*XXi+C(3)*XXi^2;

        for i=1:NDOF
            EF(i)=EF(i)+T_XXi*N(i,l)*W_i(l)*XJAC;

            for j=1:NDOF
                EK(i,j)=EK(i,j)+AE_XXi*DN(i,l)*DN(j,l)*W_i(l)/XJAC+C_XXi*N(i,l)*N(j,l)*W_i(l)*XJAC;
            end
        end
    end
end



function [GLOBF,GLOBK]=Local_Global(IPVAL,NELEM,NDOF,X,X_i,W_i,NINT,AE,C,T,N,DN,XLEN)
    
    for k=1:NELEM
        ITEMP=(k-1)*IPVAL;
        for l=1:NDOF
            JC=ITEMP+l;
            IELDOFS(l,k)=JC;
        end
    end
    %Initialize
    NTOTDOF=NELEM*IPVAL+1;
    GLOBF=zeros(NTOTDOF,1);
    GLOBK=zeros(NTOTDOF,NTOTDOF);

    %Assembly
    for IEL=1:NELEM
        [EF,EK]=Element_Cal(IEL,X,X_i,W_i,NINT,NDOF,AE,C,T,N,DN,XLEN);
        for IROW=1:NDOF
            IROWG=IELDOFS(IROW,IEL);
            GLOBF(IROWG)=GLOBF(IROWG)+EF(IROW);
            for ICOL=1:NDOF
                ICOLG=IELDOFS(ICOL,IEL);
                GLOBK(IROWG,ICOLG)=GLOBK(IROWG,ICOLG)+EK(IROW,ICOL);
            end
        end
    end
end

function [GLOBF, GLOBK]=Boundary_Conditions(GLOBF,GLOBK,BCTYPE,BC,NELEM,IPVAL)
    NTOTDOF=NELEM*IPVAL+1;

%BC at point A
    if(BCTYPE(1)==1)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,1)*BC(1);
        end
        GLOBF(1)=BC(1);
        for j=1:NTOTDOF
            GLOBK(1,j)=0;
            GLOBK(j,1)=0;
        end
        GLOBK(1,1)=1;
    elseif(BCTYPE(1)==2)
        GLOBF(1)=GLOBF(1)-BC(1);
    elseif(BCTYPE(1)==3)
        GLOBK(1,1)=GLOBK(1,1)+BC(1);
        GLOBF(1)=GLOBF(1)+BC(1)*0;
    end

%BC at point B
    if(BCTYPE(2)==1)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,NTOTDOF)*BC(2);
        end
        GLOBF(NTOTDOF)=BC(2);
        for j=1:NTOTDOF
            GLOBK(NTOTDOF,j)=0;
            GLOBK(j,NTOTDOF)=0;
        end
        GLOBK(NTOTDOF,NTOTDOF)=1;
    elseif(BCTYPE(2)==2)
        GLOBF(NTOTDOF)=GLOBF(NTOTDOF)+BC(2);
    elseif(BCTYPE(2)==3)
        GLOBK(NTOTDOF,NTOTDOF)=GLOBK(NTOTDOF,NTOTDOF)+BC(2);
        GLOBF(NTOTDOF)=GLOBF(NTOTDOF)+BC(2)*BC(3);
    end
end

%Point Load at center
function [GLOBF, GLOBK]=Point_Load(GLOBF,GLOBK,NELEM)
    if(rem(NELEM,2)==0)
        N_load=NELEM/2;
    end 
end