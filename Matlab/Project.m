%% Project part 2
clear;
load iddata-13.mat
uid=id.u; yid=id.y;
uval=val.u; yval=val.y;

m=1:6;

%% IDENTIFICATION - MSE FORMULA
nk=1;
N=length(uid);
for na=1 %for efficiency run it for each na separately
    nb=na;
    for m=4
        d_id=delay(na,nb,nk,uid,yid);
        phi_id=regressor(d_id,m);
        theta=phi_id\yid; %finding theta (linear regresion)
%         ypred=phi_id*theta;
        ysim=simulation(na,nb,nk,theta,m,uid);
        MSEs1=0;MSEs2=0;
        for k=2:N
%             MSEs1=MSEs1+(ypred(k)-yid(k))^2;
%             MSEp=MSEs1*1/(N-1);
            MSEs2=MSEs2+(ysim(k)-yid(k))^2;
            MSEs=MSEs2*1/N;
        end
%         MSEid(1,m)=MSEp;   
        MSEsim(1,m)=MSEs;   
    end
end

%% IDENTIFICATION - PREDICTION
%by running the MSE IDENTIFICATION code for every na=nb (from 1 to 3), the following MSEs are computed
MSEidpred1=[0.221284654951598,0.00926482703457547,0.00903751551445065,0.00890539067035482,0.00878475717397580,0.00877474653615681];
MSEidpred2=[0.113357021241527,1.01377706266224e-05,3.20160492415993e-06,3.05283597763085e-06,4.01084217683975e-06,3.06546755160923e-06];
MSEidpred3=[0.0580278088236740,3.44116024434305e-06,1.72248860325734e-06,1.53617601464456e-06,2.12739037306050e-05,4.46671969106262e-05];
%MSEid_pred_min=1.53617601464456e-06; mse from identification dataset, prediction (na=nb=3 and m=4) MINIMAL VALUE

%Plot MSE vs na,nb,m
plot(m,MSEidpred1); hold on; plot(m,MSEidpred2); hold on; plot(m,MSEidpred3); axis([1 6 -0.05 0.25]); grid;
xlabel('m');ylabel('MSE'); legend('na=nb=1','na=nb=2','na=nb=3');title('MSE vs na,nb,m');subtitle('Prediction values - identification');

%% IDENTIFICATION - SIMULATION
%by running the MSE IDENTIFICATION code for every na=nb (from 1 to 3), the following MSEs are computed
MSEidsim1=[0.284098901543605,0.0823452160067776,0.0817209373711655,0.0821901408262724,NaN,NaN];
MSEidsim2=[0.131859897279292,0.000344812542741304,0.000166056480583433,0.000165101816824158,NaN,NaN];
MSEidsim3=[0.0651128619798865,0.000363565061460954,NaN,NaN,NaN,NaN];
%MSEid_sim_min=1.651018168241580e-04; mse from identification dataset, simulation (na=nb=2 and m=4) MINIMAL VALUE

%Plot MSE vs na,nb,m
plot(m,MSEidsim1); hold on; plot(m,MSEidsim2); hold on; plot(m,MSEidsim3); axis([1 4 -0.05 0.3]); grid;
xlabel('m');ylabel('MSE'); legend('na=nb=1','na=nb=2','na=nb=3');title('MSE vs na,nb,m');subtitle('Simulation values - identification');

%% VALIDATION - MSE FORMULA
nk=1;
N=length(uval);
for na=3 %for efficiency run it for each na separately
    nb=na;
    for m=1:6
        d_id=delay(na,nb,nk,uid,yid);
        phi_id=regressor(d_id,m);
        theta=phi_id\yid;
        d_val=delay(na,nb,nk,uval,yval);
        phi_val=regressor(d_val,m);
        ypred=phi_val*theta;
%         ysim=simulation(na,nb,nk,theta,m,uval);
        MSEs1=0;MSEs2=0;
        for k=4:N
            MSEs1=MSEs1+(ypred(k)-yval(k))^2;
            MSEp=MSEs1*1/(N-3);
%             MSEs2=MSEs2+(ysim(k)-yval(k))^2;
%             MSEs=MSEs2*1/N;
        end
        MSEpred(1,m)=MSEp;   
%         MSEsim(1,m)=MSEs;   
    end
end

%% VALIDATION - PREDICTION WITHOUT CORRECTIONS
%by running the MSE VALIDATION code for every na=nb (from 1 to 3), the following MSEs are computed
MSEvalpred1=[0.452903587190824,0.0315696350258625,0.0379158833555720,0.0716874724790969,130.077176696005,130.172062760374];
MSEvalpred2=[0.275993524728497,49.1571564499288,690642.853439489,27676428.5292570,4186445189.08069,41143542.9418356];
MSEvalpred3=[0.185951626349582,21.0712300361415,963996850.511291,425359040.823351,457863770180.700,606941601.416999];
%MSEval_prederr_min=0.0315696350258625; mse from validation dataset, prediction with errors (na=nb=1 and m=2) MINIMAL VALUE

%Plot MSE vs na,nb,m
plot(m,MSEvalpred1); hold on; plot(m,MSEvalpred2); hold on; plot(m,MSEvalpred3); axis([1 6 -1e11 5e11]); grid;
xlabel('m');ylabel('MSE'); legend('na=nb=1','na=nb=2','na=nb=3');title('MSE vs na,nb,m');subtitle('Prediction values - validation');

%% VALIDATION - PREDICTION WITH CORRECTIONS
%by running the MSE VALIDATION code for every na=nb (from 1 to 3), the following MSEs are computed
MSEvalpred_withcorr1=[0.439121342465439,0.0258523703366705,0.0322092042973419,0.0660259790277956,0.187013632784712,0.282025953360664];
MSEvalpred_withcorr2=[0.263778629965884,7.04824994773354e-05,2.83374798784107e-05,3.05173387488735e-05,0.000122050758314170,8.85727661178942e-05];
MSEvalpred_withcorr3=[0.173510043222695,2.06436808046215e-05,2.81830143428014e-05,5.88169687184980e-05,0.00329919059365767,0.00312230218755499];
%MSEval_predwcorr_min=2.83374798784107e-05; mse from validation dataset, prediction without errors (na=nb=2 and m=3) MINIMAL VALUE
%correcting = ignoring the values (we had maximum 3) that are exeeding the range and giving large errors

plot(m,MSEvalpred_withcorr1); hold on; plot(m,MSEvalpred_withcorr2); hold on; plot(m,MSEvalpred_withcorr3); axis([1 6 -0.1 0.5]); grid;
xlabel('m');ylabel('MSE'); legend('na=nb=1','na=nb=2');title('MSE vs na,nb,m');subtitle('Prediction values (with corrections) - validation');

%% VALIDATION - SIMULATION
%by running the MSE VALIDATION code for every na=nb (from 1 to 3), the following MSEs are computed 
MSEvalsim1=[0.571773420112930,0.269172366333928,0.825941327148165,NaN,NaN,NaN];
MSEvalsim2=[0.320825680621420,NaN,NaN,NaN,NaN,NaN];
MSEvalsim3=[0.210433860206673,NaN,NaN,NaN,NaN,NaN];
%MSEval_sim_min=0.210433860206673; mse from validation dataset, simulation (na=nb=3 and m=1) MINIMAL VALUE

%TABLE

%% Representative plots identification
%% Smallest MSE prediction - for na=nb=3, m=4
na=3;nb=na;nk=1;m=4;
d_id=delay(na,nb,nk,uid,yid);
phi_id=regressor(d_id,m);
theta=phi_id\yid;
ypred=phi_id*theta;
plot(ypred); hold on; plot(yid); legend('Prediction','Identification'); title('Yprediction vs Yidentification (best fitting)');
subtitle('na=nb=3 and m=4'); grid;
%% Smallest MSE simultaion - for na=nb=2, m=4
na=2;nb=na;nk=1;m=4;
d_id=delay(na,nb,nk,uid,yid);
phi_id=regressor(d_id,m);
theta=phi_id\yid;
ypred=phi_id*theta;
ysim=simulation(na,nb,nk,theta,m,uid);
figure;
plot(ysim); hold on; plot(yid); legend('Simulation','Identification'); title('Ysimulation vs Yidentification (best fitting)');
subtitle('na=nb=2 and m=4'); grid;

%% Representative plots validation
%% Smallest MSE prediction without corrections - for na=nb=1, m=2
na=1;nb=na;nk=1;m=2;
d_id=delay(na,nb,nk,uid,yid);
d_val=delay(na,nb,nk,uval,yval);
phi_id=regressor(d_id,m);
theta=phi_id\yid;
phi_val=regressor(d_val,m);
ypred=phi_val*theta;
plot(ypred); hold on; plot(yval); legend('Prediction','Validation'); title('Yprediction vs Yvalidation (best fitting)');
subtitle('na=nb=1 and m=2'); grid;
%% Smallest MSE prediction with corrections - for na=nb=3, m=2
na=3;nb=na;nk=1;m=2;
d_id=delay(na,nb,nk,uid,yid);
d_val=delay(na,nb,nk,uval,yval);
phi_id=regressor(d_id,m);
theta=phi_id\yid;
phi_val=regressor(d_val,m);
ypred=phi_val*theta;
figure;
plot(ypred(4:end)); hold on; plot(yval(4:end)); legend('Prediction','Validation'); title('Yprediction vs Yvalidation (best fitting) - with corrections');
subtitle('na=nb=3 and m=2'); grid;
%% Smallest MSE simulation - for na=nb=2, m=3
na=3;nb=na;nk=1;m=1;
d_id=delay(na,nb,nk,uid,yid);
phi_id=regressor(d_id,m);
theta=phi_id\yid;
ysim=simulation(na,nb,nk,theta,m,uval);
figure;
plot(ysim); hold on; plot(yval); legend('Simulation','Validation'); title('Ysimulation vs Yvalidation (best fitting)');
subtitle('na=nb=3 and m=1'); grid;

%% Functions

%Delay function
%creating the delay vectors for each step k
function d=delay(na,nb,nk,u_data,y_data) 
for k=1:length(y_data)
    for n1=1:na
        if k-n1>0 %for 0 or negative indexes the value will be taken 0
            y(n1,k)=[y_data(k-n1)]; %applying the formula (recursion)--> true previous outputs
        else
            y(n1,k)=[zeros(1,1)];
        end
    end
    for n2=1:nb
        if k-nk-n2+1>0 %for 0 or negative indexes the value will be taken 0 
            u(n2,k)=[u_data(k-nk-n2+1)]; %applying the formula + considering the delay --> previous inputs
        else
            u(n2,k)=[zeros(1,1)];
        end
    end
end
de=[y;u]; %concatenating the two vectors
d=de;
end


%powers function
%creates the matrix of powers for out polynomial
function find_power=powers(delay,m)
no_cl=length(delay(:,1));
v=0:m; %vector of chosen powers (for m=2 => [0 1 2])
new_v=[];
for i=1:no_cl %algorithm to copy the given vector 'number of colums' times (ensures repetition)
    new_v=[new_v,v];
end
matrix=nchoosek(new_v,no_cl); %creates a matrix with all possible combinations (without repetition) 
%from a given vector, specifying the size of the columns
matrix=unique(matrix,'rows'); %removes the repeating rows from the matrix

i=1;
while i<length(matrix(:,1))+1 %the sum of the elements from each row must be equal to m
    if sum(matrix(i,:))~=m
        matrix(i,:)=[];
    else
        i=i+1;
    end
end

find_power=matrix;
end

%reg function
%computes the polynomial for a specific m
function reg_power=reg(delay,m)
matrix=powers(delay,m); %calling powers for constructing the powers matrix

for k=1:length(delay(1,:)) 
    for i=1:length(matrix(:,1))
        p=prod((delay(:,k))'.^matrix(i,:)); %raises each column (transposed) 
        row(k,i)=p;                         %from the delay vector at each row of the power matrix
    end
end
reg_power=row;
end

%regressor function
%computes the overall polynomial for each m
function my_regressor=regressor(delay,m)
one=ones(length(delay(1,:)),1);
g=[];
for i=1:m
    g=[g,reg(delay,i)]; %calls reg for each power m
end
g=[one,g];
my_regressor=g;
end

%simulation function
%computes ysimulation
function ysimulation=simulation(na,nb,nk,theta,m,u_data)
for k=1:length(u_data)
    for n1=1:na
        if k-n1>0 %for 0 or negative indexes the value will be taken 0 --> zero initial conditions
            y(n1,k)=[ysim(k-n1)]; %calculates y using the previous simulated output
        else
            y(n1,k)=[zeros(1,1)];
        end
    end
    for n2=1:nb
        if k-nk-n2+1>0 %for 0 or negative indexes the value will be taken 0
            u(n2,k)=[u_data(k-nk-n2+1)]; %applying the formula + considering the delay --> previous inputs
        else
            u(n2,k)=[zeros(1,1)];
        end
    end
    matrix_val(:,k)=[y(:,k);u(:,k)]; %concatenating  y and u (taking rows for each k)
    ysim(k)=regressor(matrix_val(:,k),m)*theta; %recomputing ysim for each k
end
ysim=ysim';
ysimulation=ysim;
end
