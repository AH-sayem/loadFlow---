
clc;
clear all;
close all ;

y = [8.985190-44.835953i, -3.815629+19.078144i, -5.169561+25.847809i, 0;
     -3.815629+19.078144i, 8.985190-44.835953i, 0, -5.169561+25.847809;
     -5.169561+25.847809i, 0, 8.193267-40.863838i, -3.023705+15.118528i;
     0, -5.159561+25.847809, -3.023705+15.118528i, 8.193267-40.863838i]
 
  B=[1.0000,1.0000,NaN,NaN,50.0000,30.9900;2.0000,1.0000,0,0,170.0000,105.3500;
  3.0000,1.0000,0,0,200.0000,123.9400;4.0000,1.0200,318.0000,NaN,80.0000,49.5800 ]; 
 
 z=1./y 

%y=1./z %Generate admittance matrix
for i=1:length(y) %generate Y-BUS Matrix
    for j=1:1:length(z)
        if i==j
            Y(i,j)=sum(y(i,:));
        else
            Y(i,j)=-y(i,j);
        end
    end
end

a=size(B); %compute size of given matrix
n=a(1);
for i=1:n 
    V(B(i,1))=B(i,2);
    Pg(B(i,1))=B(i,3);
    Qg(B(i,1))=B(i,4);
    Pl(B(i,1))=B(i,5);
    Ql(B(i,1))=B(i,6);
end
%e=input('Enter tolerance\n'); 
e=0.000001; %tolerance
r=10;   %maximum number of iteration
Pg(isnan(Pg))=0;
Qg(isnan(Qg))=0;
Pl(isnan(Pl))=0;
Ql(isnan(Ql))=0;
for i=1:n
    if Pg(i)>0 %generation bus voltage is constant
        Vconstant(i)=V(i);
    else 
        Vconstant(i)=0;
    end
end
lamda=angle(V);
for j=1:r
    Vprev=V;
    lamdaprev=lamda;
    fprintf('Iteration %u\n',j);
    for i=2:n 
        current=0;
        for k=1:n %calculating I with generalized equation
            current=current+(Y(i,k)*V(k));
        end
        I(i)=current;
        S(i)=V(i)*conj(I(i)); %calculating S(i)
        f(i)=real(S(i)); %taking real part of S as f
        g(i)=imag(S(i)); %taking imaginary part of S as g
        delP(i)=Pg(i)-Pl(i)-f(i); %calculating delP
        delQ(i)=Qg(i)-Ql(i)-g(i); %calculating delQ
    end
    delP(:,1)=[];  %eliminating 1st row of delP
    for i=1:n    
        if Vconstant(i)~=0 %generation bus Q 
            delQ(:,i)=[]; %eliminates corresponding rows of generation bus
        end
    end
    delQ(:,1)=[];
    PQ=[delP';delQ'];
    for i=2:n
        for k=1:n
            if i~=k %for non diagonal elements
                H(i,k)=(-1)*abs(V(i))*abs(V(k))*abs(Y(i,k))*sin(angle(Y(i,k))-lamda(i)+lamda(k));
                L(i,k)=abs(V(i))*abs(Y(i,k))*cos(angle(Y(i,k))-lamda(i)+lamda(k));
                M(i,k)=(-1)*abs(V(i))*abs(V(k))*abs(Y(i,k))*cos(angle(Y(i,k))-lamda(i)+lamda(k));
                N(i,k)=(-1)*abs(V(i))*abs(Y(i,k))*sin(angle(Y(i,k))-lamda(i)+lamda(k));
            else  %Setting diagonal elements to zero for now
                    H(i,i)=0;
                    L(i,i)=0;
                    M(i,i)=0;
                    N(i,i)=0;
            end
        end
    end
    for i=2:n      
        for k=1:n %determinig the sum of non-diagonal terms of each row
            if i~=k
                H(i,i)=H(i,i)+H(i,k);
                L(i,i)=L(i,i)+L(i,k);
                M(i,i)=M(i,i)+M(i,k);
                N(i,i)=N(i,i)+N(i,k);
            end
        end
        H(i,i)=H(i,i)*(-1); %determining final values of diagonal elements
        L(i,i)=2*abs(V(i))*abs(Y(i,i))*cos(angle(Y(i,i)))+L(i,i);
        M(i,i)=M(i,i)*(-1);
        N(i,i)=-2*abs(V(i))*abs(Y(i,i))*sin(angle(Y(i,i)))+N(i,i);
    end
    for i=1:n    
        if Vconstant(i)~=0 %LMN are edited for generation bus 
            L(:,i)=[];  %eliminating rows of L refering generation bus
            M(i,:)=[];  %eliminating columns of M refering generation bus
            N(i,:)=[];  %eliminating columns of N refering generation bus
            N(:,i)=[];  %eliminating rows of N refering generation bus
        end
    end
    H(1,:)=[];  %eliminating 1st column of H
    H(:,1)=[];  %eliminating 1st rows of H
    L(1,:)=[];  %eliminating 1st column  of L
    L(:,1)=[];  %eliminating 1st rows of L
    M(1,:)=[];  %eliminating 1st column  of M
    M(:,1)=[];  %eliminating 1st row of M
    N(1,:)=[];  %eliminating 1st column  of N
    N(:,1)=[];  %eliminating 1st row of N
    J=[H L;M N];  %forming Jaccobian matrix
    dlamdadv=inv(J)*PQ;  %calucating dellamda-delv
    lamda(1)=0;  %setting slack bus lamda to 1
    g=0; %for tracking number of generation bus
    for i=2:n
        lamda(i)=lamda(i)+dlamdadv(i-1);
        if Vconstant(i)~=0
            V(i)=Vconstant(i);
            g=g+1;   
        else
            V(i)=abs(V(i))+dlamdadv(n+i-g-2);
        end
    end
    for i=1:n
        fprintf('V(%u)=%.6f<%.6f\n',i,V(i),(lamda(i)*180/pi));
        V(i)=complex(abs(V(i))*cos(lamda(i)),abs(V(i))*sin(lamda(i)));
    end
    terminate1=1;
    terminate2=1;
    for i=2:a(1) %tolerance termination condition
        if (Vprev(i)-V(i))<e
            terminate1=1*terminate1;
        else
            terminate1=0;
        end
    end
    for i=2:a(1) %tolerance termination condition
        if (lamdaprev(i)-lamda(i))<e
            terminate2=1*terminate2;
        else
            terminate2=0;
        end
    end
    if terminate1==1 && terminate2==1
        break;
    end
end
fprintf('\nFinal result\n');
for i=1:n
    fprintf('V(%u)=%.6f<%.6f\n',i,abs(V(i)),(lamda(i)*180/pi));
end