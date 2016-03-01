clear;


function GP = Newtonp2(nbre_point,nbre_iter)
//nbre_point=513;
//nbre_iter=50;
stacksize(20000000);
rac2pi=1/(sqrt(%pi * 2));   

a=-0.5;   b=0.5; 
grille = linspace(a,b,nbre_point)';
Moy=zeros(nbre_point,1);   Var=ones(nbre_point,1);
Moy1=Moy(1:nbre_point-1);  Var1=Var(1:nbre_point-1);
for i=1:nbre_iter;
   C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
   P=cdfnor("PQ",grille,Moy,Var);
   Q=cdfnor("PQ",C,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*C.^2); expdiag=rac2pi.*exp(-0.5.*grille.^2);
   grillem=grille(1:nbre_point-1) ; grillep=grille(2:nbre_point); 
   gradient=grille.*([Q;1]-[0;Q]) + [expo;0]-[0;expo];
   
   D=-0.25.*(grille-[0;grillem]).*([0;expo])+ 0.25.*(grille-[grillep;0]).*([expo;0])+ [Q;1]-[0;Q]
   D=diag(D);
   V1=-0.25.*([0;expo]).*(grille-[0;grillem]);
   U1=-0.25.*([expo;0]).*([grillep;0]-grille);
   U=diag(U1)*diag(ones(nbre_point-1,1),1);
   V=diag(V1)*diag(ones(nbre_point-1,1),-1);
   hessien=U+D+V;
   A=hessien\gradient;
   grille=grille - A; 
end
 C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
 Q=cdfnor("PQ",C,Moy1,Var1);
 expo=rac2pi.*exp(-0.5.*C.^2); 
 distor_local = ([0;expo]).*([0;C]-2*grille)+([expo;0]).*(-[C;0]+2*grille)+(1+grille.^2).*([Q;1]-[0;Q]);
   distor= sqrt(sum(distor_local));
   
   
//C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
proba_cum=cdfnor("PQ",C,Moy1,Var1);
proba=[proba_cum(1);(proba_cum(2:nbre_point-1)-proba_cum(1:nbre_point-2));(1-proba_cum(nbre_point-1))];

//grille(nbre_point+1)=0;
//proba(nbre_point+1)=distor;

GP=[grille,proba];
endfunction



function GP = Lloyd(nbre_point,nbre_iter)
//nbre_point=513;
//nbre_iter=50;
stacksize(20000000);
rac2pi=1/(sqrt(%pi * 2));   

a=-0.5;   b=0.5; 
grille = linspace(a,b,nbre_point)';
Moy=zeros(nbre_point,1);   Var=ones(nbre_point,1);
Moy1=Moy(1:nbre_point-1);  Var1=Var(1:nbre_point-1);
for i=1:nbre_iter;
   C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
   //P=cdfnor("PQ",grille,Moy,Var);
   Q=cdfnor("PQ",C,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*C.^2); 
   grille = ([0;expo]-[expo;0])./([Q;1]-[0;Q]);
   
end
//C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
proba_cum=cdfnor("PQ",C,Moy1,Var1);
proba=[proba_cum(1);(proba_cum(2:nbre_point-1)-proba_cum(1:nbre_point-2));(1-proba_cum(nbre_point-1))];

GP=[grille,proba];
endfunction


pause;





function  dr=driftBS(x,mu,Delta)
dr = x + mu*x*Delta;
endfunction

function  R=volBS(x,sigma,Delta)
R = sigma*x*sqrt(Delta);
endfunction


function [GRILLES,PROBAS]=NewtonBS(nbre_point,nbre_pas,nbre_iter,nbre_iter1,mu,sigma,x0,T)

//t=timer();
 

Delta = T/nbre_pas;
GRILLES = zeros(nbre_point,nbre_pas+1);
PROBAS = zeros(nbre_point,nbre_pas+1);
GRILLES(:,1) = x0;
PROBAS(:,1) = 1;

rac2pi=1/(sqrt(%pi * 2));
moyx0 = driftBS(x0,mu,Delta);
moyin = moyx0;
varx0 = volBS(x0,sigma,Delta);


GP = Newtonp2(nbre_point,nbre_iter1);
 
grille0 = GP(1:nbre_point,1);

grille= varx0*grille0 + moyx0;


Moy1=zeros(nbre_point-1,1);   Var1=ones(nbre_point-1,1);
Moy1=Moy1(1:nbre_point-1);  Var1=Var1(1:nbre_point-1);

for i=1:nbre_iter;
   C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
   xpm = (C-moyx0)/varx0; 
   if (C==[]) then xpm = [];
   end
   
   Q=cdfnor("PQ",xpm,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*xpm.^2); 
   
   grillem=grille(1:nbre_point-1) ; grillep=grille(2:nbre_point); 
   
   gradient=(grille-moyx0).*([Q;1]-[0;Q]) + varx0*([expo;0]-[0;expo]);
   
   U1=-0.25.*expo.*(grillep-grillem)/varx0;
   V1=-0.25.*expo.*(grillep-grillem)/varx0;
   U=diag([U1;0])*diag(ones(nbre_point-1,1),1);
   V=diag([0;V1])*diag(ones(nbre_point-1,1),-1);
   
   D=  [Q;1]-[0;Q] + [U1;0] + [0;V1];
   D=diag(D);
   hessien=U+D+V;
   A=hessien\gradient;
   grille=grille - A; 
   
end
C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
xpm = (C-moyx0)/varx0; 
proba_cum=cdfnor("PQ",xpm,Moy1,Var1);
proba=[proba_cum(1);(proba_cum(2:nbre_point-1)-proba_cum(1:nbre_point-2));(1-proba_cum(nbre_point-1))];




GRILLES(:,2) = grille;
PROBAS(:,2) = proba;

for t=3:nbre_pas+1;
    
   moyx0 = driftBS(GRILLES(:,t-1),mu,Delta);
   varx0 = volBS(GRILLES(:,t-1),sigma,Delta);
   
   MOYX0 = ones(nbre_point-1,1)*moyx0';
   VARX0 = ones(nbre_point-1,1)*varx0';
   PROB = diag(PROBAS(:,t-1));
    
   
  
   Moy1=zeros(nbre_point-1,nbre_point);   Var1=ones(nbre_point-1,nbre_point);



for i=1:nbre_iter;
   C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
   CM = C * ones(1,nbre_point);
   xpm = (CM-MOYX0)./VARX0; 
   Q=cdfnor("PQ",xpm,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*xpm.^2); 
   
   grillem=grille(1:nbre_point-1)*ones(1,nbre_point); 
   grillep=grille(2:nbre_point)* ones(1,nbre_point); 
   grillei = grille * ones(1,nbre_point);
   UN = ones(1,nbre_point);
   ZEROS = zeros(1,nbre_point);
   
   gradient=(grillei-[MOYX0;moyx0']).*([Q;UN]-[ZEROS;Q]) + ([VARX0;varx0']).*([expo;ZEROS]-[ZEROS;expo]);
   
  
   gradient= sum(gradient*PROB,'c');
   
   U1=-0.25.*expo.*(grillep-grillem)./VARX0;
   V1= U1;
   D=  [Q;UN]-[ZEROS;Q] + [U1;ZEROS] + [ZEROS;V1];
   
   U1=sum(U1*PROB,'c');
   V1=sum(V1*PROB,'c');
   D=sum(D*PROB,'c');
   
   U=diag([U1;0])*diag(ones(nbre_point-1,1),1);
   V=diag([0;V1])*diag(ones(nbre_point-1,1),-1);
   D=diag(D);
   
   hessien=U+D+V;
   A=hessien\gradient;
   grille=grille - A; 
 
end
C = 0.5*(grille(1:nbre_point-1)+grille(2:nbre_point));
CM = C * ones(1,nbre_point);


xpm = (CM-MOYX0)./VARX0; 
Q=cdfnor("PQ",xpm,Moy1,Var1);

proba= sum(([Q;UN]-[ZEROS;Q])*PROB,'c');


GRILLES(:,t)=grille;
PROBAS(:,t) = proba;
end


endfunction


tc=timer();
x0=100;
sigma = 0.05;
T=1;
mu=0.03;
nbre_iter=10;
nbre_iter1=10;
nbre_pas = 20;
//nbre_point=50;
K=100;


//for nbre_point = 5:100

[G,P] = NewtonBS(nbre_point,nbre_pas,nbre_iter,nbre_iter1,mu,sigma,x0,T);

str=msprintf("grille%d.txt",nbre_point);

write(str,G,"(X,E20.12)")

end

pause

Prix  = exp(-mu*T)*sum( max(K-G(:,nbre_pas+1),0).*P(:,nbre_pas+1) );

tc=timer();

printf("temps de calcul  : %f s \n ",tc)

//pause



function [GRILLES,PROBAS,Dist]=NewtonBS_DISPAT(NPT,nbre_pas,nbre_iter,nbre_iter1,mu,sig,x0,T)


Delta = T/nbre_pas;
DIST=zeros(NPT(nbre_pas+1),nbre_pas);    
GRILLES = zeros(NPT(nbre_pas+1),nbre_pas+1);
PROBAS = zeros(NPT(nbre_pas+1),nbre_pas+1);
GRILLES(:,1) = x0;
PROBAS(:,1) = 1;

rac2pi=1/(sqrt(%pi * 2));
moyx0 = driftBS(x0,mu,Delta);
moyin = moyx0;
varx0 = volBS(x0,sig,Delta);

GP = Newtonp2(NPT(2),nbre_iter1);
grille0 = GP(1:NPT(2),1);

grille = varx0*grille0 + moyx0;

Moy1=zeros(NPT(2),1);   Var1=ones(NPT(2),1);
Moy1=Moy1(1:NPT(2)-1);  Var1=Var1(1:NPT(2)-1);

for i=1:nbre_iter;
   C = 0.5*(grille(1:NPT(2)-1)+grille(2:NPT(2)));
   xpm = (C-moyx0)/varx0; 
   if (C==[]) then xpm = [];
   end

   Q=cdfnor("PQ",xpm,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*xpm^2); 
 
   grillem=grille(1:NPT(2)-1) ; grillep=grille(2:NPT(2)); 
   
   gradient=(grille-moyx0).*([Q;1]-[0;Q]) + varx0*([expo;0]-[0;expo]);
   
   U1=-0.25.*expo.*(grillep-grillem)/varx0;
   V1=-0.25.*expo.*(grillep-grillem)/varx0;
   U=diag([U1;0])*diag(ones(NPT(2)-1,1),1);
   V=diag([0;V1])*diag(ones(NPT(2)-1,1),-1);
   
   D=  [Q;1]-[0;Q] + [U1;0] + [0;V1];
   D=diag(D);
   hessien=U+D+V;
   A=hessien\gradient;
   grille=grille - A; 
   
end
 
 
C = 0.5*(grille(1:NPT(2)-1)+grille(2:NPT(2)));
xpm = (C-moyx0)/varx0; 
if (C==[]) then xpm = [];
end

proba_cum=cdfnor("PQ",xpm,Moy1,Var1);
proba=[proba_cum(1);(proba_cum(2:NPT(2)-1)-proba_cum(1:NPT(2)-2));(1-proba_cum(NPT(2)-1))];
Q=cdfnor("PQ",xpm,Moy1,Var1);
expo=rac2pi.*exp(-0.5.*xpm^2);


  D0=([Q;1]-[0;Q]).*(((moyx0-grille)^2)+ varx0^2)+(([xpm;0].*[expo;0]-[0;xpm].*[0;expo]).*(-(varx0.^2))) - 2.*varx0.*(moyx0-grille).*([expo;0]-[0;expo]);
  Dist0=sum(D0);
 DIST(1:NPT(2),2)=D0;

 //pause;

GRILLES(1:NPT(2),2) = grille;
PROBAS(1:NPT(2),2) = proba;




for t=3:nbre_pas+1;
    
  
  // if  length(find(proba <0))>0 then  
  //     printf("t: %f s \n ",t);
   //    pause;
  // end
   
   moyx0 = driftBS(GRILLES(1:NPT(t-1),t-1),mu,Delta);
   varx0 = volBS(GRILLES(1:NPT(t-1),t-1),sig,Delta);
   
   MOYX0 = ones(NPT(t)-1,1)*moyx0';
   VARX0 = ones(NPT(t)-1,1)*varx0';
   

   PROB = diag(PROBAS(1:NPT(t-1),t-1));
    
    
    
   //moyx0 = driftMB(GRILLES(1:NPT(t-1),t-1));
  // varx0 = ones(GRILLES(1:NPT(t-1),t-1))*volMB(Delta);
   
   //MOYX0 = ones(NPT(t)-1,1)*moyx0';
   //VARX0 = ones(NPT(t)-1,1)*varx0';
   //PROB = diag(PROBAS(1:NPT(t-1),t-1));
 
   
   //grille = sqrt(Delta)*grille0 + moyin;
   Moy1=zeros(NPT(t)-1,NPT(t-1));   Var1=ones(NPT(t)-1,NPT(t-1));
//Moy1=Moy(1:nbre_point-1);  Var1=Var(1:nbre_point-1);
   
   //moyavant = sum(GRILLES(1:NPT(t-1)).*PROBAS(1:NPT(t-1)));
   //varavant = sum((GRILLES(1:NPT(t-1)).^2).*PROBAS(1:NPT(t-1))) - moyavant^2;
   
   //grille = sqrt((t-1)/nbre_pas) * GP(1:NPT(t),1);
   
   if (NPT(t-1) <> NPT(t)) then 
      //GP = Newtonp2(NPT(t),nbre_iter);
      Tm = (t-1)/nbre_pas;
      nbre_pas1 = Tm/Delta;
      [G,P]=NewtonBS(NPT(t),nbre_pas1,nbre_iter,nbre_iter1,mu,sig,x0,Tm);   
      //NewtonBS(nbre_point,nbre_pas,nbre_iter,nbre_iter1,mu,sigma,x0,T)
      grille =  G(1:NPT(t),nbre_pas1 + 1);
   end
   
  //pause;

for i=1:nbre_iter;
   
   C = 0.5*(grille(1:NPT(t)-1)+grille(2:NPT(t)));
   CM = C * ones(1,NPT(t-1));
   xpm = (CM-MOYX0)./VARX0; 
   //P=cdfnor("PQ",grille,Moy,Var);
   Q=cdfnor("PQ",xpm,Moy1,Var1);
   expo=rac2pi.*exp(-0.5.*xpm.^2); 
   //expdiag=rac2pi.*exp(-0.5.*grille^2);
   
   grillem=grille(1:NPT(t)-1)*ones(1,NPT(t-1)); 
   grillep=grille(2:NPT(t))* ones(1,NPT(t-1)); 
   grillei = grille * ones(1,NPT(t-1));
   UN = ones(1,NPT(t-1));
   ZEROS = zeros(1,NPT(t-1));
   

   
   gradient=(grillei-[MOYX0;moyx0']).*([Q;UN]-[ZEROS;Q]) + ([VARX0;varx0']).*([expo;ZEROS]-[ZEROS;expo]);
   
  
   gradient = sum(gradient*PROB,'c');
   
   U1=-0.25.*expo.*(grillep-grillem)./VARX0;
   V1= U1;
   D=  [Q;UN]-[ZEROS;Q] + [U1;ZEROS] + [ZEROS;V1];
   
   U1=sum(U1*PROB,'c');
   V1=sum(V1*PROB,'c');
   D=sum(D*PROB,'c');
   
   
   
   U=diag([U1;0])*diag(ones(NPT(t)-1,1),1);
   V=diag([0;V1])*diag(ones(NPT(t)-1,1),-1);
   D=diag(D);
   
   hessien=U+D+V;
   
   A=hessien\gradient;
   grille=grille - A; 
 //  pause;
end

 //pause

grillei = grille * ones(1,NPT(t-1));
C = 0.5*(grille(1:NPT(t)-1)+grille(2:NPT(t)));
CM = C * ones(1,NPT(t-1));

xpm = (CM-MOYX0)./VARX0; 
Q=cdfnor("PQ",xpm,Moy1,Var1);
expo=rac2pi.*exp(-0.5.*xpm.^2); 
proba= sum(([Q;UN]-[ZEROS;Q])*PROB,'c');

//PROB1 = diag(PROB(:,t-1));

GRILLES(1:NPT(t),t) = grille;
PROBAS(1:NPT(t),t) = proba;

  //pause;
 //D0=([Q;1]-[0;Q]).*(((moyx0-grille)^2)+ varx0^2)+(([xpm;0].*[expo;0]-[0;xpm].*[0;expo]).*(-(varx0^2))) - 2.*varx0.*(moyx0-grille).*([expo;0]-[0;expo]);
 D1=([Q;UN]-[ZEROS;Q]).*(((grillei-[MOYX0;moyx0']).^2)+[VARX0;varx0'].^2) + ([xpm;ZEROS].*[expo;ZEROS]-[ZEROS;xpm].*[ZEROS;expo]).*(-([VARX0;varx0'].^2)) - 2.*[VARX0;varx0'].*([MOYX0;moyx0']-grillei).*([expo;ZEROS]-[ZEROS;expo]);
 D11=sum(D1*PROB,'c');
 Dist=sqrt(sum(D11));
 DIST(1:NPT(t),t)=D11;
 //Err=sqrt(sum(DIST,'r'));
end
//PG=[grille,proba];
endfunction






//printf("temps de calcul  : %f s \n ",t)

nbre_pas =50;
nbre_iter=10;
nbre_iter1=5;
T=1;
Delta = T/nbre_pas;
x0=100;
mu = 0.15;
sig=0.4;
C=mu + 0.5*sig^2;
L=max(mu,sig);
kappa = 2 + 6*L;

Kp = 16*(L^3)*(4+sqrt(Delta))/sqrt(%pi*2);

t=0:Delta:1;
t(1)=[];

    
a=exp(C*(T-t)/3).*(exp((kappa + Kp)*t)*x0^3 + (exp(kappa*Delta)*L + Kp)*(exp((kappa + Kp)*t) -1)/(kappa +Kp))^(1/3);
a = sqrt(a);

N=zeros(1,nbre_pas+1);
//NPT=zeros(1,nbre_pas+1);



//for nd=1:100
   
   npt = 10001;
   nptnon=400;
   
   for i=1:nbre_pas
    N(i+1) = (a(i)*npt)/sum(a);
   end

   NPT=int(N);

   for i = 1:nbre_pas+1
   if (N(i) - NPT(i) >= 0.5) then NPT(i) = NPT(i) +1; end ;  
   end

NPT(1)=1;
   //NPT(1) = 1;
   nbre_pt = NPT(nbre_pas+1);

   
//for nbre_point=100:nbre_dist
  //  [G,P,Err]=NewtonBS_DISPAT(NPT,nbre_pas,nbre_iter,nbre_iter1,mu,sig,x0,T);
    [G1,P1]=NewtonBS(nptnon,nbre_pas,nbre_iter,nbre_iter1,mu,sig,x0,T);
   
    
//end

//Prix du Put
K=130;

payoff = max(K-G1(:,nbre_pas+1),0);

//payoff_opt = max(K-G(:,nbre_pas+1),0);

put = exp(-mu*T)*sum(payoff.*P1(:,nbre_pas+1));

//put_opt = exp(-mu*T)*sum(payoff_opt.*P(:,nbre_pas+1));









//xset("line style",1); plot2d([5:100],Y,style=2);
//xset("line style",4); plot2d([5:100],0.60987*I,style=5);
//legends(["$N \mapsto \vert Y_0 - \widehat Y_0^N\vert$";"$N \mapsto 0.61/N$"],[[2;1],[5;4]],opt=2)
//xtitle("$N \mapsto \vert Y_0 - \widehat Y_0^N\vert {\textrm{ and } N \mapsto 0.61/N }$","$N=5, \ldots,100$","")

