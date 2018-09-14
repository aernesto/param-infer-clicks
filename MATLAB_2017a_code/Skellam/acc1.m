% check values of S
clear
l1=[1,10:80]; % lambda high
minS=0;
maxS=0;


for i=1:length(l1)
    lh=l1(i);
    nl2=i-1;
    for j=1:nl2
        l2=l1(j); % low rate
        s=S(l2,lh);
        if s<minS
            minS=s;
        end
        if s>maxS
            maxS=s;
        end
    end
end
minS
maxS







function s=S(ll,lh)
% computes S from 2 rates
    s=(lh-ll)/sqrt(ll+lh);
end
