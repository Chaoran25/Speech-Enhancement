function [x]=martin1994(s,fs)
gama=0.88;%¦Ãin eq.(2£©
alpha=0.92;%¦Áin eq.(4)
omin=1.5;% omin in eq.(5)
L=length(s);%length of the signal
WFft=256;%256 points FFT,also the length of window
R=64;% R in paper
nf=ceil((L-(WFft-R))/R); % number of frames
plus=nf*R+WFft-L;
s=[s;zeros(plus,1)];
w=hamming(WFft);
% w=ones(WFft,1);%window
D=80;% same with D in paper
subf=0.04;% spectral floor constant
px=zeros(WFft/2+1,1);%power of the signal(with noise) in frame
pxn=px;
os=px;
PSFft=[];
Pactual=[];%actual Px
XX=[];
for i=0:nf-1
    ss=s(i*R+1:i*R+WFft).*w;%windowed signal
    ssfft=fft(ss,WFft);%take 256 points fft;
    ssfft=ssfft(1:end/2+1);% the former half
    XX=[XX,ssfft];
    px=ssfft.*conj(ssfft);%power of the windowed signal Pasheval theorem
    if i==0
        os=gama.*px;
        pxn=alpha.*px;
    else
        os=gama.*os+(1-gama).*px;
        pxn=alpha.*pxn+(1-alpha).*px;%eq.(4) in paper
    end
    PSFft=[PSFft,pxn];
    Pactual=[Pactual,os];
end
XX=XX';
PSFFt=PSFft';
Pactual=Pactual';
pmin=[];
Num=ceil(nf/D);
plusfra=Num*D-nf;
plusmatrix=100*ones(plusfra,WFft/2+1);
plusmatrix1=zeros(plusfra,WFft/2+1);
PSFFt=[PSFFt;plusmatrix];
Pactual=[Pactual;plusmatrix1];
XX=[XX;plusmatrix1];
Y=[];
for i=0:Num-1
    YK=[];
    pmina=PSFFt(1+i*D:(i+1)*D,:);
    Xk=Pactual(1+i*D:(i+1)*D,:);
    XX1=XX(1+i*D:(i+1)*D,:);
    pmin=min(pmina);
    Pn=omin.*pmin;
    for j=1:D
        Yk=[];
        Q1=[];
        for n=1:129
            NSR=10*log10(Pn(n)./Xk(j,n));
            if NSR>=-7;
                osub=1;
            else
                osub=0.2;
            end
            Q=1-sqrt(osub.*Pn(n)./Xk(j,n));
            Q1=[Q1,Q];
        end
        compare1=(XX1(j,:)).*Q1;
        compare2=sqrt(subf.*Pn(n));
        Yk=max(compare1,compare2);
        YK=[YK;Yk];
    end
    Y=[Y;YK];
end
Y=Y(1:nf,:);
x1=zeros(length(s),nf);
x2=[];
for i=0:nf-1
    Y1=Y(i+1,:);
    Y2=flipud(Y1(2:128));
    Y3=[Y1';Y2'];
    X1=ifft(Y3);
    X1=real(X1);
    x1(i*R+1:i*R+WFft,i+1)=X1;
end
x2=sum(x1,2);
x=x2(1:L);
b=fir1(7,0.5);
x=filter(b,1,x);
x=filter(b,1,x);
% figure;
% subplot(2,1,1);plot(s); title('orignal signal');
% str=['enhanced car 0db'];
% subplot(2,1,1);plot(x);title(str);
% % figure;
% % subplot(2,1,1);spectrogram(s);
% subplot(2,1,2);spectrogram(x);title(str);
end