
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

