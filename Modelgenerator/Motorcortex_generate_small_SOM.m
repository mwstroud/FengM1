Fclear all;close all;clc;
seed= 12345789; %(random1) % control the sequencies of generation the random values.
%seed= 3344;%(random2)
%seed=271178;%random5 5801 %random4 96045 %random3 10858;
rng(seed);

tic

%INSERT INTERNEURONS EVENLY TO FILL VOLUME, 120 SPOTS CURRENTLY
AP=0.6;%%1.4;
LAT=0.6;%1.4;
DV=0.6;%1.4;
AP_I=0.6;%1.4;   %%%%inner dimension
LAT_I=0.6;%1.4;
DV_I=0.6;%1.4;


% SET NUERON DENSITY.
bla_volume = 610;       % mm^3 from human being Rubinow et al 2016
nrn_num = 12.68e6;  % from human being Rubinow et al 2016
dens = nrn_num/bla_volume;  % cell density mm^-3.
avg_dist = round((bla_volume*1e9/nrn_num)^(1/3));   % grid length = 36 um.
dmin = 0.025;   % minimum distance between soma.


upscale=1;
PCELL=800*upscale;

INTCELL=100*upscale;   %%according to Drew, only 20% is FS.

SOMCELL=100*upscale; %%%according to Drew
                       %%%%

NCELL=PCELL+INTCELL+SOMCELL;

% TypeA=640*upscale;Anone=110*upscale;ADA=190*upscale;ANE=120*upscale; ADANE=220*upscale;    %%70,110,79,141  %%%%%according to Drew data, A took up ~65% of Pyr
% TypeB=0;%TypeB=240; Bnone=40;BDA=65;BNE=48; BDANE=87;
% TypeC=260*upscale; Cnone=57*upscale;CDA=64*upscale;CNE=50*upscale; CDANE=89*upscale;    %%27,44,32,57

% volume = NCELL/dens;
% scale = (volume/(AP*LAT*DV))^(1/3);
% scale = 1; %%%%for large network, no need to scale
% AP = AP*scale;
% LAT = LAT*scale;
% DV = DV*scale;


n=0;
m=0;
xI=(0.05:0.05:AP_I)';
yI=(0.05:0.05:LAT_I)';
zI=(0.05:0.05:DV_I)';
loca=[];loc=[];
for z=dmin:dmin:DV;
    for y=dmin:dmin:LAT;              % 0.05 seems the minimal we could make, if less than this, would be too small for placing one neuron
        for x=dmin:dmin:AP;
            if any(abs(xI-x)<=dmin/2) && any(abs(yI-y)<=dmin/2) && any(abs(zI-z)<=dmin/2)
            %if ((any(uint16(xI.*100)==uint16(x.*100)))&&(any(uint16(yI.*100)==uint16(y.*100)))&&(any(uint16(zI.*100)==uint16(z.*100))))
                m=m+1;
                loca(m,1)=x;           % for put interneuron
                loca(m,2)=y;
                loca(m,3)=z;
            else
                n=n+1;
                loc(n,1)=x;              % for put primary cell
                loc(n,2)=y;
                loc(n,3)=z;
            end
        end
    end
end

p = randperm(n,PCELL);
p=p';
% NM=zeros(NCELL,1);

location=zeros(NCELL,3);
type=zeros(NCELL,1);

index=[1:PCELL];
location(index,1:3)=loc(p(index),1:3);
type(index,1)=1;



%put in Interneuorn
p=randperm(m,INTCELL+SOMCELL);
p=p';


index=[1:INTCELL];
location(index+PCELL,1:3)=loca(p(index),1:3);
type(index+PCELL,1)=100;

%%%%SOM CELL
index=[INTCELL+1:SOMCELL+INTCELL];
location(index+PCELL,1:3)=loca(p(index),1:3);
type(index+PCELL,1)=200;


figure
x_p=location(1:PCELL,1);y_p=location(1:PCELL,2);z_p=location(1:PCELL,3);
x_I=location(PCELL+1:PCELL+INTCELL,1);y_I=location(PCELL+1:PCELL+INTCELL,2);z_I=location(PCELL+1:PCELL+INTCELL,3);
x_I_SOM=location(PCELL+INTCELL+1:NCELL,1);y_I_SOM=location(PCELL+INTCELL+1:NCELL,2);z_I_SOM=location(PCELL+INTCELL+1:NCELL,3);
scatter3(x_p,y_p,z_p,'^','filled','MarkerFaceColor',[150/255 150/255 150/255],'MarkerEdgeColor',[150/255 150/255 150/255]);hold on;
scatter3(x_I,y_I,z_I,'filled','MarkerFaceColor',[76/255 170/255 255/255],'MarkerEdgeColor',[76/255 170/255 255/255]);hold on;
scatter3(x_I_SOM,y_I_SOM,z_I_SOM,'filled','MarkerFaceColor',[239/255 62/255 61/255],'MarkerEdgeColor',[239/255 62/255 61/255]);hold on;


%%%%for PV-PV connections%%%%%%%%
                           %34% of pairs are unidirectional, 43% are reciprocal
 preID=[PCELL+1:PCELL+INTCELL];postID=[PCELL+1:PCELL+INTCELL]; 
 uni_conn=0.34; reci_conn=0.43;
 [IandIpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'balanced'); %%%balance means pre and post from same group
[unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,IandIpairs,NCELL);

%%%%for PN-PN connections%%%%%%%%
                           %10% of pairs are unidirectional, out of them 30% are reciprocal
 preID=[1:PCELL];postID=[1:PCELL]; 
 uni_conn=0.10; reci_conn=0.30;
 [PNandPNpairs]=createconnpair_coupled(preID,postID,uni_conn,reci_conn,'balanced'); %%%coupled function meaning reci is part of uni
[unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,PNandPNpairs,NCELL);

%%%%for PV-CP connections%%%%%%%%
                           %20% of pairs are unidirectional, 30% are
                           %reciprocal, rely on reciprocal to build CP->PV
 preID=[];preID=[PCELL+1:PCELL+INTCELL];postID=[];postID=[1:PCELL]; 
 uni_conn=0.2; reci_conn=0.3;
 [FSItoPpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'unbalanced'); %%%unbalance means pre and post from diff group
 [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,FSItoPpairs,NCELL) 

%%%%for SOM-CP connections%%%%%%%%
                           %26/74 of pairs are unidirectional,Beierlein et
                           %al 2003
                           %13/74 of pairs are bidirectional,Beierlein et
                           %al 2003
 preID=[];preID=[PCELL+INTCELL+1:NCELL];postID=[];postID=[1:PCELL]; 
 uni_conn=26/74; reci_conn=13/74;
 [SOMtoPpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'unbalanced'); %%%unbalance means pre and post from diff group
 [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,SOMtoPpairs,NCELL) 

%%%%for CP-SOM connections%%%%%%%%
                           %36/63 of pairs are unidirectional, Beierlein et
                           %al 2003
                           %13/63 of pairs are bidirectional,Beierlein et
                           %al 2003
 preID=[];preID=[1:PCELL];postID=[];postID=[PCELL+INTCELL+1:NCELL]; 
 uni_conn=36/63; reci_conn=0%13/63;
 [PtoSOMpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'unbalanced'); %%%unbalance means pre and post from diff group
 [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,PtoSOMpairs,NCELL) 
 
%%%%just verify reciproal CP and SOM%%%%%%%%%%%
[unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,[SOMtoPpairs;PtoSOMpairs],NCELL) 



%%%%for SOM-PV connections%%%%%%%%
                           %17/32 of pairs are unidirectional,Gibson et
                           %al 1999
                           %7/32 of pairs are bidirectional,Gibson et
                           %al 1999
 preID=[];preID=[PCELL+INTCELL+1:NCELL];postID=[];postID=[PCELL+1:PCELL+INTCELL]; 
 uni_conn=17/32; reci_conn=7/32;
 [SOMtoPVpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'unbalanced'); %%%unbalance means pre and post from diff group
 [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,SOMtoPVpairs,NCELL) 

%%%%for PV-SOM connections%%%%%%%%
                           %11/32 of pairs are unidirectional,Gibson et
                           %al 1999
                           %7/32 of pairs are bidirectional,Gibson et
                           %al 1999
 preID=[];preID=[PCELL+1:PCELL+INTCELL];postID=[];postID=[PCELL+INTCELL+1:NCELL]; 
 uni_conn=11/32; reci_conn=0%7/32;
 [PVtoSOMpairs]=createconnpair(preID,postID,uni_conn,reci_conn,'unbalanced'); %%%unbalance means pre and post from diff group
 [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,PVtoSOMpairs,NCELL) 
 
%%%%just verify reciproal PV and SOM%%%%%%%%%%%
[unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,[PVtoSOMpairs;SOMtoPVpairs],NCELL) 


entire_conn_pair=[PNandPNpairs;IandIpairs;FSItoPpairs;SOMtoPpairs;PtoSOMpairs;SOMtoPVpairs;PVtoSOMpairs];
entire_conn_pair=entire_conn_pair-1; %%%convert to NEURONID
fid=fopen('active_syn_op','w');  %%%%for saving all pre (both INH and EXC) cells for PNs&ITNs 
k=0;
for post=0:NCELL-1;
    pre=[];pre=entire_conn_pair(entire_conn_pair(:,2)==post,1);pre=sort(pre);pre=unique(pre);
%     if numel(pre)~=numel(unique(pre))
%         k=k+1;
%         rep(k,1)=post;
%     end
if numel(pre)>=1;
clear fmt; fmt=[repmat('%d\t',1,numel(pre)) '\n'];   %%%%format
fprintf(fid,fmt,pre);
end
prenum(post+1,1)=numel(pre);  %%%to store pre num for each post
end

ind_num=cumsum(prenum);ind_num=[0;ind_num];
fid_ind=fopen('active_syn_ind','w'); %%%%%for saving index for segerate fid file for each neuron
fmt=[repmat('%d\t',1,3) '\n'];
for i=1:INTCELL+PCELL+SOMCELL;
    i
        temp=[];temp=[i-1,ind_num(i),ind_num(i+1)-1];
 fprintf(fid_ind,fmt,temp);
end

%%%%%to check connection numbers%%%%%%%%%%%%%%%%%%%%%%%%%%%
for post=0:NCELL-1;
    pre=[];pre=entire_conn_pair(entire_conn_pair(:,2)==post,1);pre=sort(pre);pre=unique(pre);
pre_num(post+1,1)=numel(pre(pre<PCELL));  %%%to store EXC pre num for each post
pre_num(post+1,2)=numel(pre(pre>=PCELL&pre<PCELL+INTCELL));  %%%to store INH from PV pre num for each post
pre_num(post+1,3)=numel(pre(pre>=PCELL+INTCELL));  %%%to store INH from PV pre num for each post
end
disp(['for CP ave. #CP received:',num2str(mean(pre_num(1:PCELL,1))),'/',num2str(std(pre_num(1:PCELL,1)))]);
    disp(['for CP ave. #PV received:',num2str(mean(pre_num(1:PCELL,2))),'/',num2str(std(pre_num(1:PCELL,2)))]);
        disp(['for CP ave. #SOM received:',num2str(mean(pre_num(1:PCELL,3))),'/',num2str(std(pre_num(1:PCELL,3)))]);
            disp(['for PV ave. #CP received:',num2str(mean(pre_num(PCELL+1:PCELL+INTCELL,1))),'/',num2str(std(pre_num(PCELL+1:PCELL+INTCELL,1)))]);
                disp(['for PV ave. #PV received:',num2str(mean(pre_num(PCELL+1:PCELL+INTCELL,2))),'/',num2str(std(pre_num(PCELL+1:PCELL+INTCELL,2)))]);
                   disp(['for PV ave. #SOM received:',num2str(mean(pre_num(PCELL+1:PCELL+INTCELL,3))),'/',num2str(std(pre_num(PCELL+1:PCELL+INTCELL,3)))]);
                        disp(['for SOM ave. #CP received:',num2str(mean(pre_num(PCELL+INTCELL+1:end,1))),'/',num2str(std(pre_num(PCELL+INTCELL+1:end,1)))]); 
                            disp(['for SOM ave. #PV received:',num2str(mean(pre_num(PCELL+INTCELL+1:end,2))),'/',num2str(std(pre_num(PCELL+INTCELL+1:end,2)))]);
                                disp(['for SOM ave. #SOM received:',num2str(mean(pre_num(PCELL+INTCELL+1:end,3))),'/',num2str(std(pre_num(PCELL+INTCELL+1:end,3)))]); 


%%%%%%%oritation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=1000000;                       %%%method is from website http://mathworld.wolfram.com/SpherePointPicking.html%%%%
X=randn(N,1);
X2=X.^2;
Y=randn(N,1);
Y2=Y.^2;
Z=randn(N,1);
Z2=Z.^2;
square=1./(sqrt(X2+Y2+Z2));
[X_n]=square.*X;
[Y_n]=square.*Y;
[Z_n]=square.*Z;
randompo=[X_n,Y_n,Z_n];
nooritation = repmat([0 0 1],NCELL,1); 

P=1; %%%%portions of cells that are oritated. 0.3 is the default one

p_oritation=nooritation;
oritation_seq = randperm(N,NCELL*P);
oritation_seq=oritation_seq';
randompo=randompo(oritation_seq,:);
rand_cell_oritated=randperm(NCELL,NCELL*P);
rand_cell_oritated=rand_cell_oritated';
p_oritation(rand_cell_oritated,:)=randompo;


%%%%%%save various files%%%%%%%%%

fid=fopen('oritation.txt','w');

matrix=p_oritation;    %mm saved%*1000; %%%%%convert to um, for NEURON use.
[m,n]=size(matrix);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fid,'%g\n',matrix(i,j));
        else
            fprintf(fid,'%g\t',matrix(i,j));
        end
    end
end
fclose(fid);


fid=fopen('location.txt','w');

matrix=location;    %mm saved%*1000; %%%%%convert to um, for NEURON use.
[m,n]=size(matrix);
for i=1:m
    for j=1:n
        if j==n
            fprintf(fid,'%g\n',matrix(i,j));
        else
            fprintf(fid,'%g\t',matrix(i,j));
        end
    end
end
fclose(fid);

GG_type = ['Cell_type','.txt'];
dlmwrite(GG_type,type,'delimiter','\t','precision', '%d');


function [pairs]=createconnpair(preID,postID,uni_conn,reci_conn,flag);
%%%flag indicate whether pre and post group are from the same group.
if strcmp(flag,'balanced')
[b2,b1] = ndgrid(preID,postID); 
q = b1~=b2 & b1<b2;
B = [b2(q) b1(q)];

B = B(randperm(length(B)),:);
reci_pairs=B(1:round(length(B)*reci_conn),:);
reci_pairs=[reci_pairs;flip(reci_pairs,2)]

uni_pairs_temp=B(round(length(B)*reci_conn)+1:round(length(B)*reci_conn)+round(length(B)*uni_conn),:);

uni_pairs_temp=uni_pairs_temp(randperm(length(uni_pairs_temp)),:);

uni_pairs=[flip(uni_pairs_temp(1:round(0.5*length(uni_pairs_temp)),:),2);uni_pairs_temp(round(0.5*length(uni_pairs_temp))+1:end,:)];
pairs=[reci_pairs;uni_pairs]; %%%in pairs first column is for pre, second is for post


else
    
[b2,b1] = ndgrid(preID,postID);
q = b1~=b2;
B = [b2(q) b1(q)];

B = B(randperm(length(B)),:);
reci_pairs=B(1:round(length(B)*reci_conn),:);
reci_pairs=[reci_pairs;flip(reci_pairs,2)]

uni_pairs=B(round(length(B)*reci_conn)+1:round(length(B)*reci_conn)+round(length(B)*uni_conn),:);

% uni_pairs=[flip(uni_pairs_temp(1:round(0.5*length(uni_pairs_temp)),:),2);uni_pairs_temp(round(0.5*length(uni_pairs_temp))+1:end,:)];
pairs=[reci_pairs;uni_pairs]; %%%in pairs first column is for pre, second is for post
  
    
    
end
end

function [pairs]=createconnpair_coupled(preID,postID,uni_conn,reci_conn,flag);
%%%flag indicate whether pre and post group are from the same group.
if strcmp(flag,'balanced')
[b2,b1] = ndgrid(preID,postID); 
q = b1~=b2 & b1<b2;
B = [b2(q) b1(q)];

B = B(randperm(length(B)),:);
uni_pairs=B(1:round(length(B)*uni_conn),:);
uni_pairs=uni_pairs(randperm(length(uni_pairs)),:);%%randomize sequence
uni_pairs=[flip(uni_pairs(1:round(0.5*length(uni_pairs)),:),2);uni_pairs(round(0.5*length(uni_pairs))+1:end,:)];

reci_pairs=uni_pairs(1:round(length(uni_pairs)*reci_conn),:);
reci_pairs=flip(reci_pairs,2);

% reci_pairs_temp=B(round(length(B)*reci_conn)+1:round(length(B)*reci_conn)+round(length(B)*uni_conn),:);
% 
% uni_pairs_temp=uni_pairs_temp(randperm(length(uni_pairs_temp)),:);
% 
% uni_pairs=[flip(uni_pairs_temp(1:round(0.5*length(uni_pairs_temp)),:),2);uni_pairs_temp(round(0.5*length(uni_pairs_temp))+1:end,:)];
pairs=[reci_pairs;uni_pairs]; %%%in pairs first column is for pre, second is for post


else
    
[b2,b1] = ndgrid(preID,postID);
q = b1~=b2;
B = [b2(q) b1(q)];

B = B(randperm(length(B)),:);
uni_pairs=B(1:round(length(B)*uni_conn),:);
uni_pairs=uni_pairs(randperm(length(uni_pairs)),:);%%randomize sequence
% uni_pairs=[flip(uni_pairs(1:round(0.5*length(uni_pairs)),:),2);uni_pairs(round(0.5*length(uni_pairs))+1:end,:)];

reci_pairs=uni_pairs(1:round(length(uni_pairs)*reci_conn),:);
reci_pairs=flip(reci_pairs,2);

% reci_pairs_temp=B(round(length(B)*reci_conn)+1:round(length(B)*reci_conn)+round(length(B)*uni_conn),:);
% 
% uni_pairs_temp=uni_pairs_temp(randperm(length(uni_pairs_temp)),:);
% 
% uni_pairs=[flip(uni_pairs_temp(1:round(0.5*length(uni_pairs_temp)),:),2);uni_pairs_temp(round(0.5*length(uni_pairs_temp))+1:end,:)];
pairs=[reci_pairs;uni_pairs]; %%%in pairs first column is for pre, second is for post

  
    
    
end
end


function [unidirection_ratio,reciprocal_ratio]=verifyconn(preID,postID,pairs,NCELL);
  maxtrix=zeros(NCELL,NCELL);
  reciprocal=0;unidirection=0;
  
 for i=1:length(pairs)  %%%fill in matrix
     maxtrix(pairs(i,1),pairs(i,2))=1;
 end
 
 if   isequal(preID,postID)
     
for pre=1:numel(preID);%%%%pre
  for post=pre+1:numel(preID);%%%%post
        
      if maxtrix(preID(pre),postID(post))==1&&maxtrix(postID(post),preID(pre))==1;
          reciprocal=reciprocal+1;
      elseif maxtrix(preID(pre),postID(post))+maxtrix(postID(post),preID(pre))==1;
              unidirection=unidirection+1;
      end
  end
end
        Npair_total=length(preID)*(length(postID)-1)/2;     
      
 else
     
 for pre=1:numel(preID);%%%%pre
  for post=1:numel(postID);%%%%post
        
      if maxtrix(preID(pre),postID(post))==1&&maxtrix(postID(post),preID(pre))==1;
          reciprocal=reciprocal+1;
      elseif maxtrix(preID(pre),postID(post))+maxtrix(postID(post),preID(pre))==1;
              unidirection=unidirection+1;
      end
  end
end    
       Npair_total=length(preID)*length(postID);     
      end
          
reciprocal_ratio=reciprocal/Npair_total;
unidirection_ratio=unidirection/Npair_total;
 end

