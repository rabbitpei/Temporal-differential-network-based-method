clear;
clc;
close all;

load minmax_sd;
min_sd=minmax_sd(1);
max_sd=minmax_sd(2);
e=exp(1);
for i=1:8
    X(i,1)=rand()*(5-1)+1;
end
delta_t=1;
L=2100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.5;

p=zeros(1,25);
p(1)=-0.5;
p(2)=-0.4;
p(3)=-0.3;
p(4)=-0.2;
p(5)=-0.1;
p(6)=0.002;
%p(2)=0.03;
for i=7:25
    p(i)=(i-3)/40;
end
for i=1:25
    pr(i)=p(26-i);
end

pp(1)=0.2;
pp(2)=0.1;
pp(3)=0.005;
pp(4)=-0.1;

sample_num=10;
D= [-2 1 0 0 0 0 0 0 ;...
    -1 -1 0 0 0 0 0 0 ;...
    1 0 -1 0 0 0 0 0 ;...
    1 0 0 -1 0 0 0 0 ;...
    1 0 0 0 -1 0 0 0 ;...
    0 0 0 0 1 -1 0 1 ;...
    0 0 0 0 0 0 -1 1 ;...
    0 0 0 0 0 0 -1 -1];
T=1;
TT=50;
signif_edges=zeros(24,T);

init_trans=[0.5,0.5;0.5,0.5];
for i=1:6
    init_emiss(1,i)=1/6;
    init_emiss(2,i)=1/6;
end

%%
adjacent_network= [1 2 3 4 5 6;...
    2 1 3 4 5 6;...
    3 1 2 0 0 0;...
    4 1 2 0 0 0;...
    5 1 2 6 0 0;...
    6 1 2 5 7 8;...
    7 6 8 0 0 0;...
    8 6 7 0 0 0];
%%
pcc_thresh=0.6;
Inconsis_I=zeros(25-2,1);
for ss=1:TT
    new_seq=zeros(8,25-1);
    CC=zeros(8,4,sample_num);
    for l=1:25
        q(l)=0.96^(1/abs(pr(l)));
        E=[-2/5*q(l) 0 0 0 0 0 0 0 ;...
            0 -4/5 0 0 0 0 0 0 ;...
            0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);
        
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:8
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,2000);
        end
        TC=reshape(CC(:,l,:),8,sample_num);
        
        for na=1:8
            edges_list=[];
            neighbour_list=adjacent_network(na,2:6);
            e=0;
            for n1=2:6
                nei1=adjacent_network(na,n1);
                if nei1==0
                    break;
                end
                e=e+1;
                edges_list(e,:)=[na n1];
                for n2=n1+1:6
                    nei2=adjacent_network(na,n2);
                    if nei2==0
                        break;
                    end
                    if isempty(find(adjacent_network(nei1,:)==nei2))==0
                        e=e+1;
                        edges_list(e,:)=[n1 n2];
                    end
                end
            end
            
            
            if l>=2
                diff_edges_num=0;
                for i=1:e
                    pre_pcc=abs(corr(pre_TC(edges_list(i,1),:)',pre_TC(edges_list(i,2),:)'));
                    post_pcc=abs(corr(TC(edges_list(i,1),:)',TC(edges_list(i,2),:)'));
                    %test_edges_pcc(l,i)=post_pcc;
                    if (pre_pcc-pcc_thresh)*(post_pcc-pcc_thresh)<0
                        diff_edges_num=diff_edges_num+1;
                    end
                end
                
                for i=1:8
                    pre_sd=std(pre_TC(i,:));
                    post_sd=std(TC(i,:));
                    pre_scaled_sd=(pre_sd-min_sd)/(max_sd-min_sd);
                    post_scaled_sd=(post_sd-min_sd)/(max_sd-min_sd);
                    if (pre_scaled_sd-pcc_thresh)*(post_scaled_sd-pcc_thresh)<0
                        diff_edges_num=diff_edges_num+ std(TC(i,:));
                    end
                end
                obs_seq(l-1)=diff_edges_num;
                new_seq(na,1:l-1)=myquantile_6p(obs_seq(1:l-1),e);
            end
        end
        pre_TC=TC;
        
        if l>=3
            cell_new_seq=mat2cell(new_seq(:,1:l-2),ones(1,8))';
            [state_transi,emission]=hmmtrain(cell_new_seq,init_trans,...
                init_emiss,'ALGORITHM','Viterbi');
            pi=[0.5,0.5];
            aver_pt=0;
            for i=1:8
                [beta,pt]=pr_hmm2(new_seq(i,l-2:l-1),state_transi,emission,pi);
                p_at_t=beta(1,2);
                aver_pt=aver_pt+p_at_t/8;
            end
            Inconsis_I(l-2)=Inconsis_I(l-2)+(1-aver_pt)/TT;
            l, aver_pt
        end
    end
    ss
end

%%
plot([1:23],Inconsis_I);






