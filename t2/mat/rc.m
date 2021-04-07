
pkg load symbolic
format long



fid=fopen('data.txt');
for k=1:9
fgetl(fid);
end

fid2=fopen('data_circuit1.m', 'wt');

for k=10:16
line=fgetl(fid);
fprintf(fid2,'%s*1000;\n',line);
end

line=fgetl(fid);
fprintf(fid2,'%s;\n',line);

line=fgetl(fid);
fprintf(fid2,'%s*0.000001;\n',line);

line=fgetl(fid);
fprintf(fid2,'%s /1000;\n',line);

line=fgetl(fid);
fprintf(fid2,'%s *1000;\n',line);


fclose(fid);

fclose(fid2);

data_circuit1

G1=1/R1
G2=1/R2
G3=1/R3
G4=1/R4
G5=1/R5
G6=1/R6
G7=1/R7



 D=[1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 Kd*G6 -1;0 -G3 0 G3+G4+G5 -G5 G6 0;];
 E=[Vs;0;0;0;0;0;0];
 F=D\E
 printf("lelo");
 printf("\n\n\nNodal Method\n\n");
 
V1=F(1,1)
V2=F(2,1)
V3=F(3,1)
V5=F(4,1)
V6=F(5,1)
V7=F(6,1)
V8=F(7,1)

IR1 = (V1-V2)/R1
IR2 = (V2-V3)/R2
IR3 = (V5-V2)/R3
IR4 = V5/R4
IR5 = (V6-V5)/R5
IR6 = (-V7)/R6
IR7 = (V7-V8)/R7
IVs = IR1
Ib = -IR2
Ic = 0
Ikd = IR6


tab_file=fopen('theorcir1_TAB.tex', 'wt');
fprintf(tab_file, "$V1$ & %f $V$\\\\ \\hline\n$V2$ & %f $V$\\\\ \\hline\n$V3$ & %f $V$\\\\ \\hline\n$V5$ & %f $V$\\\\ \\hline\n$V6$ & %f $V$\\\\ \\hline\n$V7$ & %f $V$\\\\ \\hline\nV8 & %f $V$\\\\ \\hline\n$I(R_1)$ & %f $mA$\\\\ \\hline\n$I(R_2)$ & %f $mA$\\\\ \\hline\n$I(R_3)$ & %f $mA$\\\\ \\hline\n$I(R_4)$ & %f $mA$\\\\ \\hline\n$I(R_5)$ & %f $mA$\\\\ \\hline\n$I(R_6)$ & %f $mA$\\\\ \\hline\n$I(R_7)$ & %f $mA$\\\\ \\hline\n$I(V_s)$ & %f $mA$\\\\ \\hline\n$I_b$ & %f $mA$\\\\ \\hline\n$I_c$ & %f $mA$\\\\ \\hline\n$I(K_d)$ & %f $mA$\\\\ \\hline\n",V1,V2,V3,V5,V6,V7,V8,IR1*1000,IR2*1000,IR3*1000,IR4*1000,IR5*1000,IR6*1000,IR7*1000,IVs*1000,Ib*1000,Ic*1000,Ikd*1000);
fclose(tab_file);

tab_file=fopen('data_TAB.tex', 'wt');
fprintf(tab_file, "$R_1$ & %f $k\\Omega$\\\\ \\hline\n$R_2$ & %f $k\\Omega$\\\\ \\hline\n$R_3$ & %f $k\\Omega$\\\\ \\hline\n$R_4$ & %f $k\\Omega$\\\\ \\hline\n$R_5$ & %f $k\\Omega$\\\\ \\hline\n$R_6$ & %f $k\\Omega$\\\\ \\hline\n$R_7$ & %f $k\\Omega$\\\\ \\hline\n$V_S$ & %f $V$\\\\ \\hline\n$C$n& %f $uF$\\\\ \\hline\n$K_b$ & %f $mS$\\\\ \\hline\n$K_d$ & %f $mS$\\\\ \\hline\n",R1/1000,R2/1000,R3/1000,R4/1000,R5/1000,R6/1000,R7/1000,Vs,C/0.000001,Kb*1000,Kd/1000);
fclose(tab_file);


f_net=fopen("circuit1.txt","w");
fprintf(f_net, "*Circuito 1\nR1 1 2 %f \nR2 2 3 %f \nR3 5 2 %f \nR4 5 0 %f \nR5 6 5 %f \nR6 0 7 %f \nR7 4 8 %f \nC 6 8 %fu \nV 1 0 DC %f \nVf 7 4 0 \nH 5 8 Vf %f \nG 6 3 (2,5) %fm", R1,R2,R3,R4,R5,R6,R7,C,Vs,Kd,Kb);
fclose(f_net);

Vx=V6-V8

G=[1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 Kd*G6 -1;0 0 0 0 1 0 -1];

H=[0;0;0;0;0;0;Vx];


I=G\H



printf("\n\n\nNodal Method\n\n");

V1 = I(1,1) %V
V2 = I(2,1) %V
V3 = I(3,1) %V
V5 = I(4,1) %V
V6 = I(5,1) %V
V7 = I(6,1) %V
V8 = I(7,1) %V

Ix = ((V6-V5)/R5) + ((V3-V2)/R2)
REq = abs(Vx/Ix)
TAU = REq*C

tab_file=fopen('theorcir2_TAB.tex', 'wt');
fprintf(tab_file, "$V_x$ & %f $V$\\\\ \\hline\n$I_x$ & %f $mA$\\\\ \\hline\n$R_{Eq}$ & %f $k\\Omega$\\\\ \\hline\n$\\tau$ & %f $s$\\\\ \\hline",Vx, Ix*1000, REq/1000,TAU);

f_net=fopen("circuit2.txt","w");
fprintf(f_net, "*Circuito 2\nR1 1 2 %f \nR2 2 3 %f \nR3 5 2 %f \nR4 5 0 %f \nR5 6 5 %f \nR6 0 7 %f \nR7 4 8 %f \nVs 1 0 0 \nVf 7 4 0 \nH 5 8 Vf %f \nG 6 3 (2,5) %fm", R1,R2,R3,R4,R5,R6,R7,Kd,Kb);
fclose(f_net);


J = [1, 0, 0, 0, 0, 0, 0 ; -G1, G1+G2+G3, -G2, -G3, 0, 0, 0; 0, Kb+G2, -G2, -Kb, 0, 0, 0 ; -G1, G1, 0, G4, 0, G6, 0 ; 0, 0, 0, 0, 0, -G6-G7, G7 ; 0, 0, 0, 1, 0, G6*Kd, -1 ; 0, -G3, 0, G3+G4+G5, -G5, G6, 0]
K = [0; 0; 0; 0; 0; 0; 0]
L = J\K

V_11a = L(1,1)
V_22a = L(2,1)
V_33a = L(3,1)
V_44a = 0
V_55a = L(4,1)
V_66a = L(5,1)
V_77a = L(6,1)
V_88a = L(7,1)

fid = fopen("data_nule_tab.tex","w")
fprintf(fid, "$V_{1}$ & %f \\\\ \\hline \n", V_11a)
fprintf(fid, "$V_{2}$ & %f \\\\ \\hline \n", V_22a)
fprintf(fid, "$V_{3}$ & %f \\\\ \\hline \n", V_33a)
fprintf(fid, "$V_{4}$ & %f \\\\ \\hline \n", V_44a)
fprintf(fid, "$V_{5}$ & %f \\\\ \\hline \n", V_55a)
fprintf(fid, "$V_{6}$ & %f \\\\ \\hline \n", V_66a)
fprintf(fid, "$V_{7}$ & %f \\\\ \\hline \n", V_77a)
fprintf(fid, "$V_{8}$ & %f \\\\ \\hline \n", V_88a)
fclose(fid)

f_net2=fopen("circuit3.txt","w");
fprintf(f_net, "*Circuito 3\nR1 1 2 %f \nR2 2 3 %f \nR3 5 2 %f \nR4 5 0 %f \nR5 6 5 %f \nR6 0 7 %f \nR7 4 8 %f \nC 6 8 %fu \nVs 1 0 0 \nVf 7 4 0 \nH 5 8 Vf %f \nG 6 3 (2,5) %fm", R1,R2,R3,R4,R5,R6,R7,C,Kd,Kb);
fclose(f_net2);

syms t
syms A
syms wn
syms V6n(t)

V6n(t) = A*exp(wn*t)

A = Vx
wn = - (1/TAU)

t=0:1e-5:20e-3;
V6n = A*exp(wn*t)

hf = figure (1);
plot (t*1000, V6n, "g");

xlabel ("t[ms]");
ylabel ("v6_n [V]");
legend ('V6_n(t)','Location','Northeast')
print (hf, "natural.eps", "-depsc");



w = 2*pi*1000


M= [1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 G6*Kd -1;0 -G3 0 G3+G4+G5 -G5-(j*w*C) G6 j*w*C];

N= [-j; 0; 0; 0; 0; 0; 0];

O= M\N;

V1r = abs(O(1,1))
V1teta = angle(O(1,1))
V2r = abs(O(2,1))
V2teta = angle(O(2,1))
V3r = abs(O(3,1))
V3teta = angle(O(3,1))
V5r = abs(O(4,1))
V5teta = angle(O(4,1))
V6r = abs(O(5,1))
V6teta = angle(O(5,1))
V7r = abs(O(6,1))
V7teta = angle(O(6,1))
V8r = abs(O(7,1))
V8teta = angle(O(7,1))

tab_file=fopen("forced_tab.tex","w");

fprintf(tab_file, "$V_{1}$ & %f & %f\\\\ \\hline\n$V_{2}$ & %f & %f\\\\ \\hline\n$V_{3}$ & %f & %f\\\\ \\hline\n$V_{5}$ & %f & %f\\\\ \\hline\n$V_{6}$ & %f & %f\\\\ \\hline\n$V_{7}$ & %f & %f\\\\ \\hline\n$V_{8}$ & %f & %f\\\\ \\hline", V1r,V1teta*180/pi,V2r,V2teta*180/pi,V3r,V3teta*180/pi,
V5r,V5teta*180/pi,V6r,V6teta*180/pi,V7r,V7teta*180/pi,V8r,V8teta*180/pi);

fclose(tab_file);

%% Total solution
f_net=fopen("circuit4.txt","w");
fprintf(f_net, "*Circuito 4\nR1 1 2 %f \nR2 2 3 %f \nR3 5 2 %f \nR4 5 0 %f \nR5 6 5 %f \nR6 0 7 %f \nR7 4 8 %f \nC 6 8 %fu \nVs 1 0 0.0 AC 1.0 -90 SIN(0 1 1000) \nVf 7 4 0 \nH 5 8 Vf %f \nG 6 3 (2,5) %fm", R1,R2,R3,R4,R5,R6,R7,C,Kd,Kb);
fclose(f_net);

V6f = V6r*cos((w*t)+V6teta)
V6t = V6n + V6f

tneg=-5e-3:1e-5:0;
tt=[tneg,t];
Vsneg = Vs + 0*tneg
Vspos = sin(w*t)
V6neg = V6 + 0*tneg
V6pos = V6t
Vs_t = [Vsneg,Vspos];
V6_t = [V6neg,V6pos];

hf2 = figure (2);
plot (tt, Vs_t, "b");
hold on
plot (tt, V6_t, "r");

xlabel ("t[s]");
ylabel ("V_s(t) , V_6(t) [V]");
legend ('V_s(t)','V_6(t)','Location','Northeast')
print (hf2, "total.eps", "-depsc");

%% Frequency response
f = logspace(-1, 6, 50)

for i=1:1:50

w=2*pi*f(i)
A3 = [1, 0, 0, 0, 0, 0, 0 ; -G1, G1+G2+G3, -G2, -G3, 0, 0, 0; 0, Kb+G2, -G2, -Kb, 0, 0, 0 ; -G1, G1, 0, G4, 0, G6, 0 ; 0, 0, 0, 0, 0, -G6-G7, G7 ; 0, 0, 0, 1, 0, G6*Kd, -1 ; 0, -G3, 0, G3+G4+G5, -G5-(j*w*C), G6, j*w*C]
b3 = [-j; 0; 0; 0; 0; 0; 0]
V3 = A3\b3
V6(i)=V3(5)
V8(i)=V3(7)
Vs(i)=V3(1)
Vc(i)=V6(i)-V8(i)

end 

hf3 = figure (3);
plot (log10(f), 20*log10(abs(Vs)), "g");
hold on
plot (log10(f), 20*log10(abs(V6)), "r");
hold on
plot (log10(f), 20*log10(abs(Vc)), "b");

xlabel ("Frequency, in logarithmic scale [Hz]");
ylabel ("Magnitude [dB]");
legend ('Vs','V6', 'Vc','Location','Northeast')
print (hf3, "magnitude.eps", "-depsc");

hf4 = figure (4);
plot (log10(f), (180*angle(Vs))/pi, "g");
hold on
plot (log10(f), (180*angle(V6))/pi, "r");
hold on
plot (log10(f) , (180*angle(Vc))/pi, "b");

xlabel ("Frequency, in logarithmic scale [Hz]");
ylabel ("Phase [Degrees]");
legend ('Vs','V6', 'Vc','Location','Northeast')
print (hf4, "phase.eps", "-depsc");

















