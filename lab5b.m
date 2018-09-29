function x=lab5b(s,fs)
a=1;N=4;
fc=[50 150 250 350 450 570 700 840 1000 1170 1370 1600 1850 2150 2500 2900 3400];
t=1:160;
t1=160-t;
t=t/fs;
t1=t1/fs;
pt=[];
gt=[];
ERB=24.7+0.108*fc;
b=2*pi*ERB*1.019*(-1);
h=[];
for i =1:17
    pti= a.*t.^(N-1).*exp(b(i).*t).*cos(2.*pi.*fc(i).*t);
    gti= a.*t1.^(N-1).*exp(b(i).*t1).*cos(2.*pi.*fc(i).*t1);
    pt=[pt;pti];
    gt=[gt;gti];
end
Fpt=[];Fgt=[];Fqt=[];qt=[];
for i=1:17
    Fpti=fft(pt(i,:),8000);Fgti=fft(gt(i,:),8000);
    Fpt=[Fpt;Fpti];Fgt=[Fgt;Fgti];
end
for i=1:17
    Fqti=Fpt(i,:).*Fgt(i,:);
    Fqt=[Fqt;Fqti];
end
for i=1:17
    qti=ifft(Fqt(i,:));
    qti=qti./max(qti);
    qt=[qt;qti];
end
qt=qt';
L=length(s);
number=ceil(L/160);
ss=[s',zeros(1,number*160-L)]';
xm3=[];
u=100;
for i=1:17
    zf=[0];
    xmi=[];
    En=0;
    for j=0:4
        if zf==[0];
            [xmin,zf]=filter(qt(i,:),1,ss((j*160+1):((j+1)*160)));
        else
            [xmin,zf]=filter(qt(i,:),1,ss((j*160+1):((j+1)*160)),zf);          
        end   
        sum11=xmin.^2;
        sum11=sum(sum11);
        En=En+sum11;
    end
    zf=[0];
    En=En/5;
    for n=0:(number-1)
        if zf==[0];
            [xmin,zf]=filter(qt(i,:),1,ss((n*160+1):((n+1)*160)));
        else
            [xmin,zf]=filter(qt(i,:),1,ss((n*160+1):((n+1)*160)),zf);          
        end
        Esn=xmin.^2;
        Es=sum(Esn)-En;
        k=Es./(Es+u*En);
        xmin1=xmin.*k;
        xmi=[xmi;xmin1];
    end
    xm3=[xm3,xmi];
end
x=sum(xm3,2);
norm=max(x);
x=x./norm;
% str=['enhanced babble 10db'];
% figure;
% subplot(2,1,1);plot(s);
% subplot(2,1,1);plot(x);title(str);
% figure;
% subplot(2,1,1);spectrogram(s);
% subplot(2,1,2);spectrogram(x);title(str);
end