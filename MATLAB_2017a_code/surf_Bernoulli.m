q1=0:0.1:1; num_q=length(q1);
q2=q1;
[X,Y]=meshgrid(q1,q2);
S=zeros(num_q);

for i=1:num_q
    for j=1:num_q
        S(i,j)=q1(i)*q2(j)+(1-q1(i))*(1-q2(j));
    end
end
surf(X,Y,S)