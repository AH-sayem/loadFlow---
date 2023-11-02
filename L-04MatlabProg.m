clc ;
close all ;
clear all ;


%zdata=[0 , 1 ,0 ,1.0 ; 0 ,2 ,0 , 0.8 ; 1 , 2 , 0 , 0.4  ; 1 , 3 , 0  , 0.2 ; 2,3, 0 , 0.2 ; 3 , 4, 0 , 0.08]

%z=[0 , 1 ,0 ,1.0 ; 0 ,2 ,0 , 0.8 ; 1 , 2 , 0 , 0.4  ; 1 , 3 , 0  , 0.2 ; 2,3, 0 , 0.2 ; 3 , 4, 0 , 0.08]
%Y= ybus(z)

% zdata for 30 bus =[1,2,0.0192,0.0575; 1,3,0.0452,0.01852;
% 2,4,.0570,0.1737;3,4,.0132,.0379;2,5,.0472,0.1983;2,6,.0581,0.1763;
%4,6,0.0119,.0414;5,7,.0460,.1160;6,7,0.0267,.0820;6,8,.0120,.0420;
%6,9,0.0,.2080;6,10,0.0,5560;9,11,0.0,.2080;9,10,0.0,.1100;4,12,0.0,.2560;
%12,13,0.0,.1400;12,14,.1231,.02559;12,15,.0662,.1304;12,16,.0945,.1987;14,15,.2210,.1997;
%16,17,.0824,.1923;15,18,.1073,.2185;18,19,.0639,.1292;15,18,.1073,.2185;18,19,.0639,.1292;
%19,20,.0340,.0680;10,20,.0936,.2090;10,17,.0324,.0845;10,21,.0348,.0749;10,22,.0727,.1499;
%21,22,.0116,.0236;15,23,.100,.2020;22,24,.1150,.1790;23,24,.1321,.2700;24,25,.1885,.3292;
%25,26,.2544,.3800;25,27,.1093,.2087;28,27,.0,.3960;27,29,.2198,.4253;27,30,.3202,0.6027;
%29,30,.2399,.4533;8,28,.0636,.200;6,28,.0169,.0599]


zdata = input('Entert the Zdata matrix=');

%function [Y]= ybus(zdata)
nl = zdata(:,1) ; nr = zdata(:,2) ; R = zdata (:,3) ; X = zdata(:,4);
nbr= length(zdata(:,1)) ;
nbus = max(max(nl), max(nr));
Z= R+j*X ;

y = ones(nbr,1)./ Z ;
Y = zeros(nbus, nbus);
for k = 1: nbr ;
    if nl(k) >  0 & nr(k) > 0 
        Y(nl(k) , nr(k)) = Y(nl(k) , nr(k))- y(k);
        Y(nr(k) , nl(k))= Y(nl(k) , nr(k) ) ;
    end
end

for n=1 : nbus      % formation of diagonal element 
    for k =1: nbr 
        if nl(k) ==  n | nr(k) == n 
            Y(n,n) = Y(n,n) + y(k);
        else , end 
    end
end
fprintf('YBUS MATRIX:\n')
disp(Y)
%Ibus=[-j*1.1; -j*1.25 ; 0 ; 0];
ZBUS=inv(Y)
%VBUS=ZBUS*Ibus 



