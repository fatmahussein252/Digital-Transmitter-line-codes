%% Clearing workspace
clear all;
close all;

A=4;
ensemble=zeros(500,700);
st_mean=zeros(1,700);
t_mean=zeros(1,500);

Data=randi([0 1], 500, 101);
%%%% Generation of polar NRZ line code %%%
for i = 1:500
Tx=((2*Data(i,:))-1)*A; % maping for 0 to be –A, 1 to be A
Tx2=repmat(Tx,7,1);
Tx_out=reshape(Tx2,size(Tx2,1)* size(Tx2,2),1);
td=randi([1 7],1,1);
Delayed_Tx=Tx_out(td:(size(Tx_out))-(7-td)-1);
ensemble(i,:)=Delayed_Tx;
end

%%%% plot 4 realizations of polar NRZ line code %%%
 plot_line_codes(ensemble,('polar NRZ line code realisation '));
 
%%%% calculate statistical and time means of polar NRZ line code %%%
 statistical_mean(ensemble, ("statistical mean - polar NRZ line code"));
 time_mean(ensemble, ("time mean - polar NRZ line code"));

%%%% Generation of Unipolar line code %%%
for i = 1:500
Tx=A*Data(i,:); % maping for 0 to be 0, 1 to be A
Tx2=repmat(Tx,7,1);
Tx_out=reshape(Tx2,size(Tx2,1)* size(Tx2,2),1);
td=randi([1 7],1,1);
Delayed_Tx=Tx_out(td:(size(Tx_out))-(7-td)-1);
ensemble(i,:)=Delayed_Tx;
end

%%%% plot 4 realizations of Unipolar line code %%%
 plot_line_codes(ensemble,('Unipolar line code realisation '));
 
%%%% calculate statistical and time means of Unipolar line code %%%
 statistical_mean(ensemble,  ("statistical mean - Unipolar line code"));
 time_mean(ensemble,   ("time mean - Unipolar line code"));
   
   
%%%% Generation of polar RZ line code %%%
for i = 1:500
Tx=((2*Data(i,:))-1)*A; % maping for 0 to be –A, 1 to be A
Tx2=repmat(Tx,7,1);
Tx2(6:7,:)=0;
Tx_out=reshape(Tx2,size(Tx2,1)* size(Tx2,2),1);
td=randi([1 7],1,1);
Delayed_Tx=Tx_out(td:(size(Tx_out))-(7-td)-1);
ensemble(i,:)=Delayed_Tx;
end

%%%% plot 4 realizations of polar RZ line code %%%
 plot_line_codes(ensemble,('polar RZ line code realisation '));
 
%%%% calculate statistical and time means of polar RZ line code %%%
 statistical_mean(ensemble, ("statistical mean - polar RZ line code"));
 time_mean(ensemble, ("time mean - polar RZ line code"));

function plot_line_codes(ensemble,str)
    figure;
    for i=1:4
    subplot(4,1,i);
    plot(ensemble(i,:));
    title(str, num2str(i));
    ylim([-5 5]);
    end
end
   
function statistical_mean(ensemble,str)
    for i=1:700
     st_mean(1,i)=sum(ensemble(:,i))/500;
    end
    figure;
    subplot(2,1,1);
    plot(st_mean);
    title(str);
    ylim([-5 5]);
end

function time_mean(ensemble,str)
    
    for i=1:500
     t_mean(1,i)=sum(ensemble(i,:))/700;
    end
    subplot(2,1,2);
    plot(t_mean);
    title(str);
    ylim([-5 5]);

end
