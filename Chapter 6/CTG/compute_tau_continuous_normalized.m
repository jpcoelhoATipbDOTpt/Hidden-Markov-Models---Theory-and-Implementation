function tau=compute_tau_continuous_normalized(Alfa,Beta,A,B)

[m,N]=size(B);
tau=zeros(m,m);

for k=1:N-1,
    tmp1=A.*(Alfa(k,:).'*Beta(k+1,:)).*(B(:,k+1)*ones(1,m)).';
    tmp2=ones(1,m)*tmp1*ones(m,1);
    tau=tau+tmp1/tmp2;
end

