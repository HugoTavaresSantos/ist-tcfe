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


G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;



N=[1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 Kd*G6 -1;0 -G3 0 G3+G4+G5 -G5 G6 0;];

x=[Vs;0;0;0;0;0;0];



solnodes=N\x;



printf("\n\nNodal Method\n");
printf("\n\nNodal Method\n");
V1 = solnodes(1,1) %V
V2 = solnodes(2,1) %V
V3 = solnodes(3,1) %V
V5 = solnodes(4,1) %V
V6 = solnodes(5,1) %V
V7 = solnodes(6,1) %V
V8 = solnodes(7,1) %V
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


G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;
G5=1/R5;
G6=1/R6;
G7=1/R7;

Vx=V6-V8

N=[1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 Kd*G6 -1;0 0 0 0 1 0 -1];

x=[0;0;0;0;0;0;Vx];


solnodes=N\x;



printf("\n\nNodal Method\n");

V1 = solnodes(1,1) %V
V2 = solnodes(2,1) %V
V3 = solnodes(3,1) %V
V5 = solnodes(4,1) %V
V6 = solnodes(5,1) %V
V7 = solnodes(6,1) %V
V8 = solnodes(7,1) %V

Ix = ((V6-V5)/R5) + ((V3-V2)/R2)
REq = abs(Vx/Ix)
TAU = REq*C

tab_file=fopen('theorcir2_TAB.tex', 'wt');
fprintf(tab_file, "$V_x$ & %f $V$\\\\ \\hline\n$I_x$ & %f $mA$\\\\ \\hline\n$R_{Eq}$ & %f $k\\Omega$\\\\ \\hline\n$\\tau$ & %f $s$\\\\ \\hline",Vx, Ix*1000, REq/1000,TAU);
fclose(tab_file);

fid=fopen('data.txt');
for k=1:9
fgetl(fid);
end
C = textscan(fid,'%s = %s');
fclose(fid);



fid=fopen('../sim/circuit1.txt', 'wt');

fprintf(fid, 'Vs 1 0 DC %s\n', C{2}{8});
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 9 8 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Cap 6 8 %sUF\n', C{2}{9});

fclose(fid);


Vx = V6 - V8;

fid=fopen('../sim/circuit2.txt','wt');

fprintf(fid, 'Vs 1 0 0V\n');
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 8 9 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Vx 6 8 %.15f', Vx);

fclose(fid);


fid=fopen('data.txt');
for k=1:9
fgetl(fid);
end
C = textscan(fid,'%s = %s');
fclose(fid);

fid=fopen('../sim/circuit3.txt','wt');

fprintf(fid, 'Vs 1 0 DC 0\n');
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 9 8 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Cap 6 8 %sUF\n', C{2}{9});

fclose(fid);

fid=fopen('../sim/circuit_ic.txt', 'wt');

fprintf(fid, '.IC v(6)=%fV v(8)=%fV\n',V6,V8);

fclose(fid);


fid=fopen('../sim/circuit4.txt','wt');

fprintf(fid, 'Vs 1 0 0.0 AC 1.0 -90 SIN(0 1 1000)\n');
fprintf(fid, 'R1 1 2 %sK\n', C{2}{1});
fprintf(fid, 'R2 2 3 %sK\n', C{2}{2});
fprintf(fid, 'R3 5 2 %sK\n', C{2}{3});
fprintf(fid, 'R4 5 0 %sK\n', C{2}{4});
fprintf(fid, 'R5 6 5 %sK\n', C{2}{5});
fprintf(fid, 'R6 0 7 %sK\n', C{2}{6});
fprintf(fid, 'R7 9 8 %sK\n', C{2}{7});
fprintf(fid, 'GIb 6 3 (2,5) %sm\n', C{2}{10});
fprintf(fid, 'VIc 7 9 0V\n');
fprintf(fid, 'HVd 5 8 VIc %sK\n', C{2}{11});
fprintf(fid, 'Cap 6 8 %sUF\n', C{2}{9});

fclose(fid);

syms t
syms A
syms wn
syms V6n(t)

V6n(t) = A*exp(wn*t)

A = Vx;
wn = - (1/TAU);

t=0:1e-5:20e-3;
V6n = A*exp(wn*t);

hf = figure (1);
plot (t*1000, V6n, "g");

xlabel ("t [ms]");
ylabel ("V 6_n [V]");
legend ('V 6_n(t)','Location','Northeast')
print (hf, "natural.eps", "-depsc");


data_circuit1

w = 2*pi*1000

N = [1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 G6*Kd -1;0 -G3 0 G3+G4+G5 -G5-(j*w*C) G6 j*w*C];

x = [-j; 0; 0; 0; 0; 0; 0];

solnodes = N\x;

V1r = abs(solnodes(1,1))
V1phase = angle(solnodes(1,1))
V2r = abs(solnodes(2,1))
V2phase = angle(solnodes(2,1))
V3r = abs(solnodes(3,1))
V3phase = angle(solnodes(3,1))
V5r = abs(solnodes(4,1))
V5phase = angle(solnodes(4,1))
V6r = abs(solnodes(5,1))
V6phase = angle(solnodes(5,1))
V7r = abs(solnodes(6,1))
V7phase = angle(solnodes(6,1))
V8r = abs(solnodes(7,1))
V8phase = angle(solnodes(7,1))

tab_file=fopen("forced_tab.tex","w");

fprintf(tab_file, "$V_{1}$ & %f & %f\\\\ \\hline\n$V_{2}$ & %f & %f\\\\ \\hline\n$V_{3}$ & %f & %f\\\\ \\hline\n$V_{5}$ & %f & %f\\\\ \\hline\n$V_{6}$ & %f & %f\\\\ \\hline\n$V_{7}$ & %f & %f\\\\ \\hline\n$V_{8}$ & %f & %f\\\\ \\hline", V1r,V1phase*180/pi,V2r,V2phase*180/pi,V3r,V3phase*180/pi,
V5r,V5phase*180/pi,V6r,V6phase*180/pi,V7r,V7phase*180/pi,V8r,V8phase*180/pi);

fclose(tab_file);



V6f = V6r*cos((w*t)+V6phase);
V6t = V6n + V6f;

negative_t=-5e-3:1e-5:0;
total_t=[negative_t,t];
Vs_negative = Vs + 0*negative_t;
Vs_positive = sin(w*t);
V6_negative = V6 + 0*negative_t;
V6_positive = V6t;
Vs_total = [Vs_negative,Vs_positive];
V6_total = [V6_negative,V6_positive];

hf = figure (2);
plot (total_t, Vs_total, "b");
hold on
plot (total_t, V6_total, "r");

xlabel ("t [s]");
ylabel ("V_s(t) , V_6(t) [V]");
legend ('V_s(t)','V_6(t)','Location','Northeast')
print (hf, "total.eps", "-depsc");


freq = logspace(-1, 6, 30);

for i=1:1:30

w=2*pi*freq(i);

N = [1 0 0 0 0 0 0;-G1 G1+G2+G3 -G2 -G3 0 0 0;0 Kb+G2 -G2 -Kb 0 0 0;-G1 G1 0 G4 0 G6 0;0 0 0 0 0 -G6-G7 G7;0 0 0 1 0 G6*Kd -1;0 -G3 0 G3+G4+G5 -G5-(j*w*C) G6 j*w*C];

x = [-j; 0; 0; 0; 0; 0; 0];

solnodes = N\x;
V6(i)=solnodes(5,1);
V8(i)=solnodes(7,1);
Vs(i)=solnodes(1,1);
Vc(i)=V6(i)-V8(i);

end 

hf = figure (3);
plot (log10(freq), 20*log10(abs(Vs)), "g");
hold on
plot (log10(freq), 20*log10(abs(V6)), "r");
hold on
plot (log10(freq), 20*log10(abs(Vc)), "b");

xlabel ("Frequency, in logarithmic scale [Hz]");
ylabel ("Magnitude [dB]");
legend ('Vs','V6', 'Vc','Location','Northeast')
print (hf, "magnitude.eps", "-depsc");

hf = figure (4);
plot (log10(freq), (180*angle(Vs))/pi, "g");
hold on
plot (log10(freq), (180*angle(V6))/pi, "r");
hold on
plot (log10(freq) , (180*angle(Vc))/pi, "b");

xlabel ("Frequecy, in logarithmic scale [Hz]");
ylabel ("Phase [Degrees]");
legend ('Vs','V6', 'Vc','Location','Northeast')
print (hf, "phase.eps", "-depsc");

