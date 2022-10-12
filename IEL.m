clear
close all

kbi=3*(10^6);%/M/s , based on the hybridization rate constant fitted by Zhang and Winfree
DG_assoc=1.9;%kcal/mol
DG_bp=-1.7;%kcal/mol
R=8.314;%J/mol*K
T=310.15;%Kelvin
u0=1;%M
DGsp=7.3;%kcal/mol
DGs=5.3;%kcal/mol

kuni=kbi*exp((((-abs(DG_bp)/2)+DG_assoc)*4184)/(R*T));
l=exp(-(abs(DG_bp)-DG_assoc)/(R*T));
kr(1)=kbi*l;
kback(1,1)=kr(1);
pzip=kuni/(kuni+kr(1));

h=input('Insert the toehold length:');
b=input('Insert the branch migration length:');
DGinitial=input('Insert the initial free energy Gibbs:');
DGfinal=input('Insert the final free energy Gibbs:');

u(1)=0.1;

for i=1:10
    kforward(i,1)=kbi*u(i);
    u(i+1)=u(i)+0.5;
end


for j=1:10
    for i=2:h
        kfray=kuni*exp(-(i-1)*abs(DG_bp)/(R*T));
        kr(i)=kfray*(i-1)*(1-pzip);
        kback(1,i)=kr(i);
        kforward(j,i)=kuni;
    end
end


for i=1:10
    DGbox(i,1)=DGinitial;
    %DGinit
    DGbox(i,2)=DGbox(i,1)+DG_assoc+((R*T*log(u0/u(i)))/4184);
    %DGtoehold
    for j=1:h-2
        DGbox(i,j+2)=DGbox(i,j+1)-(R*T*log((kuni)/kr(j))/4184);
    end
    kfirst=0.5*kuni*exp(-DGsp/(R*T));
    kforward(i,1+h)=kfirst;
    %DGfirst, toehold is fully bounded- state C to C+1
    DGbox(i,h+1)=DGbox(i,h)+DGsp;
    %state C+1 to D
    n=1/exp(DG_assoc/(R*T));
    keff=(kuni*exp(-(DGsp-h*abs(DG_bp))/(R*T))*n)/(2*b);
    for k=b-1:-1:1
         r=rem(k,2);
         if r==0
             DGbox(i,h+1+b-k)=DGbox(i,h+b-k)+DGs;
         else
             DGbox(i,h+1+b-k)=DGbox(i,h+b-k)-DGs;
         end
    end
end

%step D to E
for i=1:10
    DGbox(i,h+b)=DGfinal;
    kforward(i,h+1)=kbi*u(i);
end


%Branch migration time step
kbm=kuni*exp(-DGs/(R*T));
tbm=1/kbm;

%calculation pi

for i=1:10
        pforward(i,1)=kforward(i,1)/(kforward(i,1)+kback(1,1));
        pback(i,1)=kback(1,1)/(kforward(i,1)+kback(1,1));
end


for i=2:h
   pforward(1,i)=kforward(1,i)/(kforward(1,i)+kback(1,i));
   pback(1,i)=kback(1,i)/(kforward(1,i)+kback(1,i));
end 

%plot of DGbox
figure(1);
for i=1:10
    plot(DGbox(i,:));
    hold on
end
xlabel('number of steps');
ylabel('DGbox(kcal/mol)');

figure(2);
for i=1:10
    plot(kforward(i,:));
    hold on
end
xlabel('number of steps');
ylabel('kforward(/M/s)');

figure(3);
plot(kback(1,:));
hold on
xlabel('number of steps');
ylabel('kback(/M/s)');
