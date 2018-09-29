function [x]=S1979(s,fs);
L=0.032*fs;
w=hamming(L);
overlap=L/2;
len=length(s);
number=ceil((len-overlap)/(L-overlap));
buqi=number*L-((number-1)*overlap)-len;
ss=[s;zeros(buqi,1)];
Energy=[];
sss3=[];
for i=0:number-1
    sss=ss(i*overlap+1:i*overlap+L);
    sss1=sss.*w;
    sss2=fft(sss1);
    FrameEnergy=sum(sss1.^2);
    Energy=[Energy,FrameEnergy];
    sss3=[sss3,sss2];
end
sss4=abs(sss3);
m=1:5;
ll=length(m);
Ndash=[];
sss6=[];
for k=1:L
    sum1=0;
    for i=1:ll
        sss5=sss4(:,i);
        sss6=[sss6,sss5];
        sum1=sum1+sss5(k);
    end
    sum2=sum1/ll;
    Ndash=[Ndash,sum2];
end
sss6=sss6(:,1:ll);
Ndash1=Ndash';
one1=ones(L,1);
HR=[];
for i=1:number
    sss5=sss4(:,i);
    H=one1-Ndash1./sss5;
    Hr=(H+abs(H))/2;
    HR=[HR,Hr];
end
sestimate=HR.*sss3;
Nr=max(sss6,[],2);
for i=1:number
    t1=0;
    T=0;
    for k=1:L
        if sestimate(k,i)>Nr(k)
            continue;
        else
            if i==1
                sestimate(k,i)=min(sestimate(k,i),sestimate(k,i+1));
            elseif i==number
                sestimate(k,i)=min(sestimate(k,i),sestimate(k,i-1));
            else
                sestimate(k,i)=min(sestimate(k,i),sestimate(k,i-1));
                sestimate(k,i)=min(sestimate(k,i),sestimate(k,i+1));
            end
        end     
        t=abs(sestimate(k,i)/Ndash(k));
        t1=t1+t;
    end
    t1=t1/L;
    T=20*log10(t1);
    if T<=-40
        sestimate(:,i)=0.0316*sss3(:,i);
    end
end
[m,n]=size(ss);
x=zeros(m,number);
for i=1:number
    xx=ifft(sestimate(:,i));
    x((i-1)*overlap+1:(i-1)*overlap+L,i)=xx;
end
x=sum(x,2);
x=x(1:len);
b=[1,0.9];
x=filter(b,1,x); % smooth
% str=['enhanced babble 5db'];
% figure;
% subplot(2,1,1);plot(s);
% subplot(2,1,1);plot(x);title(str);
% figure;
% subplot(2,1,1);spectrogram(s);
subplot(2,1,2);spectrogram(x);title(str);
end