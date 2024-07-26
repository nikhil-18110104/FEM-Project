clc;clear
%initializing variables
EIzz=2000/12;
L=0.1;
q=240;
NELEM_array=[1,4,10,50,100];


%Boundary Conditions
% TYPE 1.disp 2.slope 3.shear_force 4.Bending_Moment 5.Moment_at_center
BCTYPE=[1,2,3,4]; %
BC=[0,0,600,10];
folder_path='/Users/nikhil/Desktop/FEM/A3/';
for P=1:length(NELEM_array)
    NELEM=NELEM_array(P);
    h_I=L/NELEM;
    X_node=0:h_I:L;

    [GLOBF,GLOBK]=Local_Global(NELEM,X_node,q,EIzz,h_I);
    [GLOBK,GLOBF]=Boundary_Condition(NELEM,GLOBK,GLOBF,BCTYPE,BC);
    
    U1=GLOBK\GLOBF;
    %FEM solution
    [F,DF,M,V]=solution(NELEM,X_node,h_I,U1,EIzz);
    
    %Exact solution
    syms Uexact(x)
    diff_eq = diff(EIzz*(diff(Uexact,x,2)),x,2) - q == 0;
    
    k=[0 0 L L L/2];
    for i=1:4
        if(BCTYPE(i)==1)
            COND(i) = subs(Uexact,x,k(i)) == BC(i);
        elseif(BCTYPE(i)==2)
            COND(i) = subs(diff(Uexact,x,1),x,k(i)) == BC(i);
        elseif(BCTYPE(i)==3)
            COND(i) = subs(-diff(EIzz*diff(Uexact,x,2),x,1),x,k(i)) == BC(i);
        elseif(BCTYPE(i)==4)
            COND(i) = subs(EIzz*diff(Uexact,x,2),x,k(i)) == -BC(i);
        elseif(BCTYPE(i)==5)
            COND(i) = subs(EIzz*diff(Uexact,x,2),x,k(i+1)) == +BC(i);
        end
    end
    
    Uexact = dsolve(diff_eq, COND(1), COND(2),COND(3),COND(4));
    Uexact=simplify(Uexact);
    
    %plotting
    X_p=0:L/100:L;
    [U_p,DU_p]=plotting(NELEM,F,DF,X_p,X_node);
    [M_p,V_p]=plotting(NELEM,M,V,X_p,X_node);
    Uexact_eval = subs(Uexact,x,X_p);
    DUexact_eval = subs(diff(Uexact),x,X_p);
    Mexact_eval = subs(EIzz*diff(Uexact,x,2),x,X_p);
    Vexact_eval = subs(-diff(EIzz*diff(Uexact,x,2),x),x,X_p);
    
    figure;
    plot(X_p*100,U_p*1000,'LineWidth',2);
    hold on;
    plot(X_p*100,Uexact_eval*1000,'LineWidth',2,'Linestyle','--');
    xlabel('x (cm)');
    ylabel('u (mm)');
    title('Transverse displacement Variation');
    text1=['Exact'];
    text2=['FEM for ',num2str(NELEM_array(P)),' elements'];
    legend(text2,text1);
    grid on;
    file_name=['U_Plot',num2str(P),'.png'];
    saveas(gcf,fullfile(folder_path,file_name));
    
    figure;
    plot(X_p*100,DU_p,'LineWidth',2);
    hold on;
    plot(X_p*100,DUexact_eval,'LineWidth',2,'Linestyle','--');
    xlabel('x (cm)');
    ylabel('du/dx');
    title('Transverse displacement slope Variation');
    text1=['Exact'];
    text2=['FEM for ',num2str(NELEM_array(P)),' elements'];
    legend(text2,text1);
    grid on;
    file_name=['DU_Plot',num2str(P),'.png'];
    saveas(gcf,fullfile(folder_path,file_name));
    
    figure;
    plot(X_p*100,V_p,'LineWidth',2);
    hold on;
    plot(X_p*100,Vexact_eval,'LineWidth',2,'Linestyle','--');
    xlabel('x (cm)');
    ylabel('V (N)');
    title('Shear Force Variation');
    text1=['Exact'];
    text2=['FEM for ',num2str(NELEM_array(P)),' elements'];
    legend(text2,text1);
    grid on;
    file_name=['V_Plot',num2str(P),'.png'];
    saveas(gcf,fullfile(folder_path,file_name));
    
    figure;
    plot(X_p*100,M_p,'LineWidth',2);
    hold on;
    plot(X_p*100,Mexact_eval,'LineWidth',2,'Linestyle','--');
    xlabel('x (cm)');
    ylabel('M (N-m)');
    title('Moment Variation');
    text1=['Exact'];
    text2=['FEM for ',num2str(NELEM_array(P)),' elements'];
    legend(text2,text1);
    grid on;
    file_name=['M_Plot',num2str(P),'.png'];
    saveas(gcf,fullfile(folder_path,file_name));

    figure;
    plot(X_p*100,M_p*(-6),'LineWidth',2);
    hold on;
    plot(X_p*100,Mexact_eval*(-6),'LineWidth',2,'Linestyle','--');
    xlabel('x (cm)');
    ylabel('Bending Stress (MPa)');
    title('Beding stress Variation');
    text1=['Exact'];
    text2=['FEM for ',num2str(NELEM_array(P)),' elements'];
    legend(text2,text1);
    grid on;
    file_name=['BS_Plot',num2str(P),'.png'];
    saveas(gcf,fullfile(folder_path,file_name));
end


function [U_p,DU_p]=plotting(NELEM,F,DF,X_p,X_node)
    U_p=zeros(length(X_p),1);
    DU_p=zeros(length(X_p),1);
    syms X;
    for i=1:length(X_p)
        for j=1:NELEM
            if(X_p(i)==X_node(NELEM+1))
                U_p(i)=subs(F(j),X,X_p(i));
                DU_p(i)=subs(DF(j),X,X_p(i));

            elseif(X_node(j)<=X_p(i) && X_p(i)<X_node(j+1))
                U_p(i)=subs(F(j),X,X_p(i));
                DU_p(i)=subs(DF(j),X,X_p(i));
            end
        end
    end
end

function [F,DF,M,V]=solution(NELEM,X_node,h_I,U1,EIzz)
    for i=1:NELEM
        syms X
        xk=X_node(i);
        xk1=X_node(i+1);
        XJAC=(xk1-xk)/2;
        N1 = 1 - (3*((X-xk)/h_I)^2) + (2*((X-xk)/h_I)^3);
        N2 = -(X-xk)*(1-((X-xk)/h_I))^2;
        N3 = (3*((X-xk)/h_I)^2) - (2*((X-xk)/h_I)^3);
        N4 = -(X-xk)*(((X-xk)/h_I)^2 - (X-xk)/h_I);
        DN1 = -(6*(X - xk)*(h_I - X + xk))/h_I^3;
        DN2 = - ((X - xk)/h_I - 1)^2 - (2*((X - xk)/h_I - 1)*(X - xk))/h_I;
        DN3 = (6*(X - xk)*(h_I - X + xk))/h_I^3;
        DN4 = ((X - xk)*(2*h_I - 3*X + 3*xk))/h_I^2;
        D2N1 = -(6*(h_I - 2*X + 2*xk))/h_I^3;
        D2N2 = (2*(2*h_I - 3*X + 3*xk))/h_I^2;
        D2N3 = (6*(h_I - 2*X + 2*xk))/h_I^3;
        D2N4 = (2*(h_I - 3*X + 3*xk))/h_I^2;
        D3N1 = 12/h_I^3;
        D3N2 = -6/h_I^2;
        D3N3 = -12/h_I^3;
        D3N4 = -6/h_I^2;

        F(i)=(N1*U1(2*i-1)+N2*U1(2*i)+N3*U1(2*i+1)+N4*U1(2*i+2));
        DF(i)=(DN1*U1(2*i-1)+DN2*U1(2*i)+DN3*U1(2*i+1)+DN4*U1(2*i+2));
        M(i)=EIzz*(D2N1*U1(2*i-1)+D2N2*U1(2*i)+D2N3*U1(2*i+1)+D2N4*U1(2*i+2));
        V(i)=-EIzz*(D3N1*U1(2*i-1)+D3N2*U1(2*i)+D3N3*U1(2*i+1)+D3N4*U1(2*i+2));
    end
end


function [GLOBF,GLOBK]=Local_Global(NELEM,X_node,q,EIzz,h_I)
    GLOBK=zeros(2*(NELEM+1));
    GLOBF=zeros(2*NELEM+2,1);
    
    for k=1:NELEM
        %reating shape functions 
        syms X
        xk=X_node(k);
        xk1=X_node(k+1);
        N1 = 1 - (3*((X-xk)/h_I)^2) + (2*((X-xk)/h_I)^3);
        N2 = -(X-xk)*(1-((X-xk)/h_I))^2;
        N3 = (3*((X-xk)/h_I)^2) - (2*((X-xk)/h_I)^3);
        N4 = -(X-xk)*(((X-xk)/h_I)^2 - (X-xk)/h_I);
        DN1 = -(6*(h_I - 2*X + 2*xk))/h_I^3;
        DN2 = (2*(2*h_I - 3*X + 3*xk))/h_I^2;
        DN3 = (6*(h_I - 2*X + 2*xk))/h_I^3;
        DN4 = (2*(h_I - 3*X + 3*xk))/h_I^2;
        N=[N1 N2 N3 N4];
        DN=[DN1 DN2 DN3 DN4];
        EK=zeros(4);
        EF=zeros(4,1);
    
        %Elemental calculations with numerical integration
        for i=1:4
            f=@(x)subs(q*N(i),X,x);
            Fel= Integrate(f,xk,xk1,10);
            EF(i,1)=EF(i,1)+Fel;
            for j=1:4
                f=@(x)subs(EIzz*DN(i)*DN(j),X,x);
                Kel= double(Integrate(f,xk,xk1,10));
                EK(i,j)=EK(i,j)+Kel;
            end
        end
    
        %Assembly
        for i=1:2
            i1=2*i-1;
            i2=2*i;
            I1=2*(k-1)+i1;
            I2=2*(k-1)+i2;
            GLOBF(I1,1)=GLOBF(I1,1)+EF(i1,1);
            GLOBF(I2,1)=GLOBF(I2,1)+EF(i2,1);
            for j=1:2
                j1=2*j-1;
                j2=2*j;
                J1=2*(k-1)+j1;
                J2=2*(k-1)+j2;
                GLOBK(I1,J1)=GLOBK(I1,J1)+EK(i1,j1);
                GLOBK(I1,J2)=GLOBK(I1,J2)+EK(i1,j2);
                GLOBK(I2,J1)=GLOBK(I2,J1)+EK(i2,j1);
                GLOBK(I2,J2)=GLOBK(I2,J2)+EK(i2,j2);
            end
        end
    end
end

function [GLOBK,GLOBF]=Boundary_Condition(NELEM,GLOBK,GLOBF,BCTYPE,BC)
    NTOTDOF=NELEM*2+2;
    %BC 1 at point A
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
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,2)*BC(1);
        end
        GLOBF(2)=BC(1);
        for j=1:NTOTDOF
            GLOBK(2,j)=0;
            GLOBK(j,2)=0;
        end
        GLOBK(2,2)=1;

    elseif(BCTYPE(1)==3)
        GLOBF(1)=GLOBF(1)+BC(1);

    elseif(BCTYPE(1)==4)
        GLOBF(2)=GLOBF(2)+BC(1);
    end

    %BC 2 at point A
    if(BCTYPE(2)==1)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,1)*BC(2);
        end
        GLOBF(2)=BC(2);
        for j=1:NTOTDOF
            GLOBK(1,j)=0;
            GLOBK(j,1)=0;
        end
        GLOBK(1,1)=1;

    elseif(BCTYPE(2)==2)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,2)*BC(2);
        end
        GLOBF(2)=BC(2);
        for j=1:NTOTDOF
            GLOBK(2,j)=0;
            GLOBK(j,2)=0;
        end
        GLOBK(2,2)=1;

    elseif(BCTYPE(2)==3)
        GLOBF(1)=GLOBF(1)+BC(2);

    elseif(BCTYPE(2)==4)
        GLOBF(2)=GLOBF(2)+BC(2);
    end

    %BC 1 at point B
    if(BCTYPE(3)==1)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,NTOTDOF-1)*BC(3);
        end
        GLOBF(NTOTDOF-1)=BC(3);
        for j=1:NTOTDOF
            GLOBK(NTOTDOF-1,j)=0;
            GLOBK(j,NTOTDOF-1)=0;
        end
        GLOBK(NTOTDOF-1,NTOTDOF-1)=1;

    elseif(BCTYPE(3)==2)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,NTOTDOF)*BC(3);
        end
        GLOBF(NTOTDOF)=BC(3);
        for j=1:NTOTDOF
            GLOBK(NTOTDOF,j)=0;
            GLOBK(j,NTOTDOF)=0;
        end
        GLOBK(NTOTDOF,NTOTDOF)=1;
    elseif(BCTYPE(3)==3)
        GLOBF(NTOTDOF-1)=GLOBF(NTOTDOF-1)+BC(3);
    elseif(BCTYPE(3)==4)
        GLOBF(NTOTDOF)=GLOBF(NTOTDOF)+BC(3);
    end

    %BC 2 at point B
    if(BCTYPE(4)==1)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,NTOTDOF-1)*BC(4);
        end
        GLOBF(NTOTDOF-1)=BC(2);
        for j=1:NTOTDOF
            GLOBK(NTOTDOF-1,j)=0;
            GLOBK(j,NTOTDOF-1)=0;
        end
        GLOBK(NTOTDOF-1,NTOTDOF-1)=1;

    elseif(BCTYPE(4)==2)
        for j=1:NTOTDOF
            GLOBF(j)=GLOBF(j)-GLOBK(j,NTOTDOF)*BC(4);
        end
        GLOBF(NTOTDOF)=BC(4);
        for j=1:NTOTDOF
            GLOBK(NTOTDOF,j)=0;
            GLOBK(j,NTOTDOF)=0;
        end
        GLOBK(NTOTDOF,NTOTDOF)=1;
    elseif(BCTYPE(4)==3)
        GLOBF(NTOTDOF-1)=GLOBF(NTOTDOF-1)+BC(4);
    elseif(BCTYPE(4)==4)
        GLOBF(NTOTDOF)=GLOBF(NTOTDOF)+BC(4);
    elseif(BCTYPE(4)==5)
        GLOBF((NTOTDOF/2))=GLOBF((NTOTDOF/2))-BC(4);
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
