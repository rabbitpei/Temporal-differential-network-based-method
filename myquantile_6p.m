function [quantiled_seq] = myquantile_6p( obs_seq,total_edges)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes her
   quantiled_seq=[];
   tmp_seq=[1:total_edges];
   fields=quantile(tmp_seq,[.167,.333,.5,.667,.833]);
   for i=1:length(obs_seq)
       if obs_seq(i)<fields(1)
           quantiled_seq(i)=1;
       elseif obs_seq(i)>=fields(1)&&obs_seq(i)<fields(2)
           quantiled_seq(i)=2;
       elseif obs_seq(i)>=fields(2)&&obs_seq(i)<fields(3)
           quantiled_seq(i)=3;
       elseif obs_seq(i)>=fields(3)&&obs_seq(i)<fields(4)
           quantiled_seq(i)=4;
       elseif obs_seq(i)>=fields(4)&&obs_seq(i)<fields(5)
           quantiled_seq(i)=5;
       else 
           quantiled_seq(i)=6;
       end
   end
end

