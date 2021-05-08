close all
clear all

pkg load symbolic;
format long;
%By Ricardo Rodrigues Vitor Liotti and Hugo Santos
%All OUR CHOSEN VALUES
C = 0.0005
R1 = 1300
R2 = 5300
f=50;
Vs = 230;
n = 16;
Vr = Vs/n;
num_diodes = 20;
von = 0.6;
I_s = 1e-14;
Vt = 0.025;
Req = 1/((1/R1)+(1/R2));
T = 1/(2*f); 
w = 2*pi*f;

%ESTIMATION FOR TOFF USING NEWTON RAPHSODY
t_OFF = (1/4)*T;
for i = 1:20
  f = Vr*C*w*sin(w*t_OFF) - (Vr/R1)*cos(w*t_OFF) - I_s*(exp(12/(Vt*num_diodes))-1);
  df = Vr*C*(w^2)*cos(w*t_OFF)+(Vr/R1)*w*sin(w*t_OFF);
  t_OFF = t_OFF - (f/df);
end


%ESTIMATION FOR TON USING NEWTON RAPHSODY
t_ON = (3/4)*T;
for i = 1:20
  g = Vr*cos(w*t_ON)+Vr*cos(w*t_OFF)*exp(-(1/(Req*C))*(t_ON-t_OFF));
  dg = -w*Vr*sin(w*t_ON)-Vr*cos(w*t_OFF)*(1/(Req*C))*exp(-(1/(Req*C))*(t_ON-t_OFF));
  t_ON = t_ON - (g/dg);
end

t = 0:(1e-6):0.4;

l = length(t);

vO = zeros(1,l);

for i = 1:l
  if t(i)<=t_OFF
    vO(i) = abs(Vr*cos(w*t(i)));
  else
    if t(i)<=t_ON
      vO(i) = Vr*abs(cos(w*t_OFF))*exp(-((t(i)-t_OFF)/(Req*C)));
    else
      t_OFF = t_OFF + T;
      t_ON = t_ON + T;
      if t(i)<=t_OFF
        vO(i) = abs(Vr*cos(w*t(i)));
      else
        if t(i)<=t_ON
          vO(i) = Vr*abs(cos(w*t_OFF))*exp(-((t(i)-t_OFF)/(Req*C)));
        end
      end
    end
  end    
end


dc_vO = mean(vO)
ripple_env = max(vO) - min(vO);
vO_aux = (ripple_env/2) + min(vO);

vr = zeros(1, length(t));
vr = Vr * cos(w*t);
for i=1:length(t)
	vr(i) = abs(vr(i));
endfor

rd = (Vt)/(I_s*exp((12/num_diodes)/(Vt))); %DIODES RESISTANCE
ac_vO_fin = ((num_diodes*rd)/(num_diodes*rd + R2))*(vO - dc_vO);

if vO_aux >= 12
    dc_vO_fin = 12;
else
    dc_vO_fin = vO_aux;
end

vO_fin = ac_vO_fin + dc_vO_fin;

ripple_env = max(vO) - min(vO)
average_env = mean(vO)

average_reg = mean(vO_fin)
ripple_reg = max(vO_fin)-min(vO_fin)




%PLOTS
	
%OUTPUTS AND GRAPHS OF VOLTAGE RECTIFIER, ENVELOPE DETECTOR AND REGULATOR TERMINALS
hf = figure(1);
plot(t*1000, vr,"k");
hold on
plot(t*1000,vO, "b");
hold on
plot(t*1000,vO_fin,"r");
xlabel ("t[ms]");
ylabel ("vO [V]");
legend('vr', 'vO envelope', 'vO regulator', 'Location', 'northeast');
print (hf, "output.eps", "-depsc");
	
%COMPONENTS OF AC
hf = figure(2);
plot (t*1000, vO-12, "b");
hold on
plot(t*1000, vO_fin-12,"r");
xlabel ("t[ms]");
ylabel ("vO [V]");
legend('vO envelope - 12V', 'vO regulator - 12V',  'Location', 'northeast');
print (hf, "ac_component.eps", "-depsc");



%TABLES

printf ("chosen_values_TAB\n");
printf ("R1 = %e \n", R1);
printf ("R2 = %e \n", R2);
printf ("C = %e \n", C);
printf ("chosen_values_END\n\n");

printf ("envelope_TAB\n");
printf ("enveloperipple = %e \n", ripple_env);
printf ("envelopeaverage = %e \n", average_env);
printf ("envelope_END\n\n");

printf ("regulator_TAB\n");
printf ("regulatorripple = %e \n", ripple_reg);
printf ("regulatoraverage = %e \n", average_reg);
printf ("regulator_END\n\n");


