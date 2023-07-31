clear all
clc
deposition = [ 1 1 1 1 1 1 ]
x=0.01;
   
height = length(deposition(:,1));
width = length(deposition(1,:));
convenience = zeros(height,2);
z = 0;
for i = height:-1:1
 z = z + 1;
 for j = 1:width
 if deposition(i,j) ~= 0
 convenience(z,1) = convenience(z,1) + 1;
 end
 end
end

m = length(convenience(:,1)); % m=number of layers in the object
%n = Number of individual filaments
n = 0;
n_new = 0;
for j = 1:m
 if m == 1
 n = convenience(1,1);
 else
 if convenience(j,2) <= 1
 n_new = n + convenience(j,1);
 n = n_new;
 end
 end
end
%%Problem definition and properties
time_step = 0.05; %Step time
erro = 0.001; %Convergence error
%Size of the variables
h(1,5)=0;
lambda(1,5)=0;
attributes(n,5)=0;
T(n,5)=0;
b_component(n,5)=0;
Q_component(n,5)=0;
b(1,n)=0;
Q(1,n)=0;
T_begin(1,n)=0;
dif(1,n)=0;
Temp_save(1,n)=0;
Temp_old(1,n)=0;
save_time(1,n)=0;
recaller = zeros(11,n);

%Process Variables
T_filament = 300; %Initial Temperature of the filament
T_Envi = 77; %Environment Temperature
v = 0.02; % Extrusion head
for lin = 1:n % SUpport material temperature
 T(lin,5) = T_Envi;
end
%Dimensions of filament
w = 0.00013; %Layer
L = 0.02; %Length
area = pi * (w/2)^2; %Area of the cross section
perimeter = pi * w; %Perimeter of the cross section of filament
% Material Properties
%Thermal conductivity (W/m.K)
conductivity = 0.18;

%Density (kg/m^3)
ro=  1150; 

%Specific heat (J/kg.K)
C = 2200;
% Heat transfer coefficient (lost of heat by natural convection)
h_convection = 45;
%Thermal contact conductances between
h(1,1) = 200; % left filament
h(1,2) = 200; % down filament
h(1,3) = 200; % right filament
h(1,4) = 200; % top filament
h(1,5) = 10; % support filament
%Fraction of perimeter contact between
lambda(1,1) = 0.3; % left filament
lambda(1,2) = 0.3; % down filament
lambda(1,3) = 0.3; % right filament
lambda(1,4) = 0.3;% top filament
lambda(1,5) = 0.3; % support filament

%Definition of the parameters influenced by the contacts
for col = 1:5
 for lin = 1:n
 vec_b(lin,col) = h(1,col)*lambda(1,col);
 vec_Q(lin,col) = vec_b(lin,col)*T(lin,col);
 end
end
z = 0;
number_filament = 0;
for i = height:-1:1
    for j = width:-1:1
         if deposition(i,j) ~= 0
            number_filament = number_filament + 1;
            Percent_contact(number_filament) = -perimeter/(ro*area*C);
            k(number_filament) = conductivity;
        end
 end
end
%calculation of temperatures at different points
for i = 1:(n+2)
 if odd(i) == 1
 temp_break(i,1) = (i*L-x)/v;
 temp_break(i,2) = (i*L+x)/v;
 else
 temp_break(i,1) = temp_break(i-1,2);
 temp_break(i,2) = ((i+1)*L-x)/v;
 end
end
for road = 1:n
 list = 0;
 for i = 0:time_step:temp_break(n,2)
 list = list + 1;
 temp(list,road) = T_filament;
 end
end
for layer = 1:m
 if layer == 1
 for num = 1:convenience(layer,1)
 if num == 1
 attributes(num,5) = 1; %Activation of the contact with support
 %Definition of the variables b and Q defined in equation Eq. 7
 b(num) = h_convection*(1-lambda*attributes(num,:)') + vec_b(num,:)*attributes(num,:)';
 Q(num) = (h_convection*(1-lambda*attributes(num,:)')*T_Envi + vec_Q(num,:)*attributes(num,:)')/b(num);
p = 0;
 for t = 0:time_step:temp_break(num,1)
 p = p+1;
 timer(p) = t;
 end
 for t = (temp_break(num,1)+time_step):time_step:temp_break(num,2)
 p = p+1; timer(p) = t;
 temp(p,num)=(T_filament-Q(num))*exp(Percent_contact(num)*b(num)*(t-temp_break(num,1)))
 +Q(num);
 end
 %Saving the last temperature of the period time of cooling down
 T_begin(num) = temp(p,num);
else
 %Activation of the contacts
 attributes(num-1,3) = 1;
 attributes(num,1) = 1;
 attributes(num,5) = 1;

 %Up-dating of the variable b
 for j = 1:num
 b(j) = h_convection*(1-lambda*attributes(j,:)') + vec_b(j,:)*attributes(j,:)';
 end
 if m == 1
 if num == convenience(layer,1)
 time_final = temp_break(num,2)+10;
 else
 time_final = temp_break(num,2);
 end
 else
 time_final = temp_break(num,2);
 end

 for t = (temp_break(num,1)+time_step):time_step:time_final
 p = p+1; timer(p) = t;
 last = p-1;
 for j = 1:num
 Temp_save(j) = temp(last,j);
 end

 %Iterative process untill steady state
 for q = temp_break(1,2):time_step:temp_break(4,2)
 %Storing previous temperatures
 for j = 1:num
 if j == 1
 T(j,3) = Temp_save(j+1);
 recaller(3,j) = j+1;
 end
 if j > 1 & j < num
 T(j,1) = Temp_save(j-1);
 recaller(1,j) = j-1;
 T(j,3) = Temp_save(j+1);
 recaller(3,j) = j+1;
 end
 if j == num
 T(j,1) = Temp_save(j-1);
 recaller(1,j) = j-1;
 end
 for k = 1:5
 if T(j,k) ~= 0 & k ~= 5
     vec_Q(j,k) = vec_b(j,k)*T(j,k);
 end
 end

 Q(j) = (h_convection*(1-lambda*attributes(j,:)')*T_Envi + vec_Q(j,:)*attributes(j,:)')/b(j);
 Temp_old(j) = Temp_save(j);
 end

 if num == 2
 Temp_save(1) = (T_begin(1)-Q(1))*exp(Percent_contact(1)*b(1)*(t-temp_break(1,2)))+Q(1);
 Temp_save(2) = (T_filament-Q(2))*exp(Percent_contact(2)*b(2)*(t-temp_break(1,1)))+Q(2);
 save_lim(1,1) = temp_break(num,2);
 save_lim(1,2) = temp_break(num,2);
 else
 for j=1:num-2
 Temp_save(j) = (T_begin(j)-Q(1))*exp(Percent_contact(j)*b(j)*(t-save_lim(1,j)))+Q(j);
 end
 Temp_save(num-1) = (T_begin(num-1)-Q(num-1))*exp(Percent_contact(num-1)*b(num-1)*(t-temp_break(num,1)))+Q(num-1);
 Temp_save(num) = (T_filament-Q(num))*exp(Percent_contact(num)*b(num)*(t-temp_break(num,1)))+ Q(num);
 save_lim(1,num-1) = temp_break(num,1);
 save_lim(1,num) = temp_break(num,1);
 end
 for j = 1:num
 dif(j) = abs(Temp_save(j)-Temp_old(j));
 end
 tru = 1;
 stop = 0;
 for j = 1:num
 if dif(tru) < erro
 tru = tru+1;
 end
if tru == num+1;
 stop = 1;
 end
 end
 if stop == 1
 for j = 1:num
 temp(p,j) = Temp_save(j);
 end
 end
 end
 end
 T_begin(num) = temp(p,num);
 plot(temp)
 legend('Filament1','Filament2','Filament3','Filament4','Filament5','Filament6')
 title('Temperature Profile For Filaments at x = ',x)
 xlabel('Time Steps, (1 time step = 0.05s)')
 ylabel('Temperature, Celsius')
 end
 end
 end
 
end

function tf = odd(x)
  if ~isreal(x)
    error('iseven:badinput','iseven requires real inputs');
  else
    tf = mod(x,2)==1;
  end
end
