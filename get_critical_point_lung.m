clc;
clear;
close all;

load minmax_sd;
min_sd=minmax_sd(1);
max_sd=minmax_sd(2);

fid=fopen('adj_edges.txt');
adjacent_network={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_network{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
total_node_num=j;
%%
fid=fopen('adj_idx_edges.txt');
adjacent_idx={};
j=0;
while ~feof(fid)
    tline=fgetl(fid);
    j=j+1;
    adjacent_idx{j}=regexp(tline, '\t', 'split');
end
fclose(fid);
%%
fpi=fopen('lung_injury_control.txt');
hline1 = textscan(fpi, '%s', 1, 'delimiter', '\n');
field1=textscan(hline1{1}{1},'%s');
clear format;
format='%s';
for i=2:57
    format=[format,' %f'];
end
plines =textscan(fpi, format,1000000,'delimiter', '\t');
pipi=plines{1};
pprofile = [];
for i = 2 : 57
    pprofile = [pprofile, plines{i}];
end
fclose(fpi);

fpi=fopen('lung_injury_case.txt');
hline2 = textscan(fpi, '%s', 1, 'delimiter', '\n');
field2=textscan(hline2{1}{1},'%s');
clear format;
format='%s';
for i=2:55
    format=[format,' %f'];
end
mlines =textscan(fpi, format,1000000,'delimiter', '\t');
mipi=mlines{1};
mprofile = [];
for i = 2 : 55
    mprofile = [mprofile, mlines{i}];
end
fclose(fpi);

allprofile(:,1:56)=pprofile;
allprofile(:,57:110)=mprofile;
allprofile=zscore(allprofile');
allprofile=allprofile';
controlprofile=allprofile(:,1:56);
caseprofile=allprofile(:,57:110);
for i=1:9
    tempcontrol(:,i,:)=controlprofile(:,6*i-5:6*i);
    tempcase(:,i,:)=caseprofile(:,6*i-5:6*i);
end
tempcase(:,2:10,:)=tempcase;
tempcase(:,1,:)=tempcontrol(:,1,:);
psize=size(tempcase);
%%
pcc_thresh=0.9;
Inconsis_I=zeros(10-2,1);
aver_change=zeros(10-1,1);
test_pcc=zeros(10,1);

for l=1:10
    for na=1:total_node_num
        edges_list=[];
        %neighbour_list=adjacent_network{na};
        center=adjacent_network{na}{1};
        e=0;
        for n1=2:length(adjacent_network{na})
            nei1=adjacent_network{na}{n1};
            idx1=adjacent_idx{na}{n1};
            e=e+1;
            edges_list(e,:)=[str2num(center) str2num(nei1)];
            for n2=n1+1:length(adjacent_network{na})
                nei2=adjacent_network{na}{n2};
                idx2=adjacent_idx{na}{n2};
                if isempty(find(str2num(char(adjacent_idx{str2num(idx1)}))==str2num(idx2)))==0
                    e=e+1;
                    edges_list(e,:)=[str2num(nei1) str2num(nei2)];
                end
            end
        end
        if l>=2
            diff_edges_num=0;
            for i=1:e
                pre_pcc=abs(corr(reshape(tempcase(edges_list(i,1),l-1,:),psize(3),1),...
                    reshape(tempcase(edges_list(i,2),l-1,:),psize(3),1)));
                post_pcc=abs(corr(reshape(tempcase(edges_list(i,1),l,:),psize(3),1),...
                    reshape(tempcase(edges_list(i,2),l,:),psize(3),1)));
                test_edges_pcc(l,i)=post_pcc;
                if (pre_pcc-pcc_thresh)*(post_pcc-pcc_thresh)<0
                    diff_edges_num=diff_edges_num+1;
                end
            end
            test_pcc(l)=test_pcc(l)+mean(test_edges_pcc(l,1:e));
            
            obs_seq(l-1)=diff_edges_num;
            new_seq(na,1:l-1)=myquantile_6p(obs_seq(1:l-1),e);
            aver_change(l-1)=aver_change(l-1)+new_seq(na,l-1);
        end
    end
    test_pcc(l)=test_pcc(l)/na;
    if l>=3
        cell_new_seq=mat2cell(new_seq(:,1:l-2),ones(1,total_node_num))';
        init_trans=[0.5,0.5;0.5,0.5];
        for i=1:6
            init_emiss(1,i)=1/6;
            init_emiss(2,i)=1/6;
        end
        [state_transi,emission]=hmmtrain(cell_new_seq,init_trans,...
            init_emiss,'ALGORITHM','Viterbi');
        pi=[0.5,0.5];
        aver_pt=0;
        for i=1:total_node_num
            [beta,pt]=pr_hmm2(new_seq(i,l-2:l-1),state_transi,emission,pi);
            p_at_t=beta(1,2);
            aver_pt=aver_pt+p_at_t/total_node_num;
        end
        Inconsis_I(l-2)=Inconsis_I(l-2)+(1-aver_pt);
        l, aver_pt
    end
    l
end
plot([2:9],Inconsis_I);
test_pcc