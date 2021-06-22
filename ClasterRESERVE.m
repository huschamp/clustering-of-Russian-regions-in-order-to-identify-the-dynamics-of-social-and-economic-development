clear;
clc;
% X = mvnrnd([0;0],[1 .9;.9 1],1000); 
 % Y = [1 1;1 -1;-1 1;-1 -1]; 
 % d2_mahal = mahal(X(:,1),X(:,2)) ;
regions=caseread('C:\big\four\regionst.txt');
p1=dlmread('C:\big\four\doh.txt');
p2=dlmread('C:\big\four\ves.txt');
p3=dlmread('C:\big\four\prir.txt');
p4=dlmread('C:\big\four\bezr.txt');
p5=dlmread('C:\big\four\inter.txt');
table=[p1 p2 p3 p4 p5];
X=table;
[t1 ps]=mapstd(X);%%\
%%X=t1;%%%%%%%%%%/ normirovka neverna (mogno estestvennuy)
%X=linkage(X,'average','mahalanobis')
X1=X;
for j=1:39
    for i=1:5
       minn=min(X(j,:));
       maxx=max(X(j,:));
       sr=sum(X(j,:))/5;
       X1(j,i)=(X(j,i)-minn)/(maxx-minn);
    end;
end;

% X1=X
% for i=1:5
%     for j=1:39
%        minn=min(X(:,i));
%        maxx=max(X(:,i));
%        sr=sum(X(:,i))/39;
%        %X1(j,i)=(X(j,i)-minn)/(maxx-minn);
%        X1(j,i)=X(j,i)/maxx;
%     end;
% end;
X=X1;
d1_mahal = mahal(X(:,1),X(:,1)) ;
d2_mahal = mahal(X(:,2),X(:,2)) ;
d3_mahal = mahal(X(:,3),X(:,3)) ;
d4_mahal = mahal(X(:,4),X(:,4)) ;
d5_mahal = mahal(X(:,5),X(:,5)) ;
d12_mahal=mahal(X(:,1),X(:,2));
d_mahal_table=[d1_mahal d2_mahal d3_mahal d4_mahal d5_mahal];
%X=d_mahal_table;
result = d2_mahal.*X(:,2);
X=table;
[idx,C] = kmeans(X,6); 
% figure
%  for i=1:5
%      for j=1:5
          figure
         gscatter(X(:,1),X(:,2),idx, 'bgmkrc' )
%      end;
%  end;
 %gscatter(d1_mahal,d2_mahal,idx, 'bgm' )
 hold on
 plot(C(:,1),C(:,2), 'kx' )
 legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3', 'Cluster Centroid' ) 
 %%%
%Z = linkage(X, 'weighted' , 'mahalanobis' ); 
%%%Z = linkage(X, 'average' , 'mahalanobis' );
Z = linkage(X);
%Z = linkage(X, 'average' , 'euclidean' ); 
%Z=linkage(X);
T = cluster(Z, 'maxclust' ,3);
q=Z(end-2,3);
cutoff = median([Z(end-5,3) Z(end-4,3)]);
%dendrogram(Z, 'ColorThreshold' ,cutoff) 
dendrogram(Z,0) 
%%%%%%%%%%%%%%%%%%%%%%%metod blish soseds
idxb=[1];
X=table;
while length(idxb)~=39
    idxb=[idxb 0];
end;
n=0;
for j=2:39
    rass=[];
    %k=[X(1,:);X(2,:)]
    %c=pdist(k)
    for i=1:39
        if i==j
            continue;
        end;
        k=[X(j,:);X(i,:)]; %%%просто двоеточие
        rass=[rass pdist(k)];
    end;
    [A m]=min(rass);
    if m>j
        m=m+1;
    end;
    while idxb(m)==0
        [A m]=min(rass);
        if idxb(m)==0  %%%%%%%%%%%во всех idbx - m+1
            rass(m)=[];
        end;
    end;
    if A<=3600 %%2000
        if idxb(m)==1
            idxb(j)=1;
        end;
        if idxb(m)==2
            idxb(j)=2;
        end;
        if idxb(m)==3
            idxb(j)=3;
        end;
        if idxb(m)==4
            idxb(j)=4;
        end;
         if idxb(m)==5
            idxb(j)=5;
        end;
         if idxb(m)==6
          idxb(j)=6;
        end;
      %  idxb(m+1)=1;
    else
      if n==3 && A >20000
         idxb(j)=4;
      end;
         if n==1
             idxb(j)=3;
             n=n+1;
             continue;
         end;
        if n==0
            idxb(j)=2;
            n=n+1;
        end;
         if n==2 && A>10000
            idxb(j)=4;
            n=n+1;
            continue;
        end;
        if n==3 && A>5000
            idxb(j)=5;
            n=n+1;
            continue;
        end;
         if n==4
            idxb(j)=6;
            %n=n+1;
            continue;
        end;
        %n=n+1;
    end;
end;
 figure
 gscatter(X(:,1),X(:,2),idxb', 'bgmkrc' )
 %gscatter(d1_mahal,d2_mahal,idx, 'bgm' )
 %hold on
 %plot(C(:,1),C(:,2), 'kx' )
 legend( 'Cluster 1' , 'Cluster 2', 'Cluster 3', 'Cluster 4','Cluster 5','Cluster 6' ) 
 %[idxbb,C2] = kmeans(X,6); 
 %figure           %%%%%%%%% FOR PROVERKA
 %gscatter(X(:,1),X(:,2),idxbb, 'bgmkrc' )
%%%%%%%%%%%%%%%%%%%%%%metod  sred svyzi
 
 
 
 
 
 flex=[1 2 1 3 5 4 1 3 6 7];
 y=find(flex == 1);
 %%%%%%%%%%%%%%%%%%%metod  dalnego soseda
 idxd=1:39;
 id=idxd([1 3 6])
 index=[0 0];
% u=[];
 u2=[];
 for j=1:39
     u=[];
     rass=[];
    for i=1:39
        if ismember(i,u)==1
            continue;
        end;
        if i==j
            continue;
        end;
        if length(find(idxd==i))>1
            u2=find(idxd==i)
            u=[u u2];
            k=[X(j,:)];
            for h=1:length(u2)
                k=[k;X(u2(h),:)];
            end;
            rassd=pdist(k);
            [B z]=max(rassd);
            rass=[rass B];
            continue;
        end;
        k=[X(j,:);X(i,:)]; %%%просто двоеточие
        rass=[rass pdist(k)];
    end;
    [A m]=min(rass);
     if m>j
         m=m+1;
     end;
     if ismember(m,index(:,2))==1 || ismember(m,index(:,1))==1
         idxd(j)=idxd(m);
     else
         idxd(m)=idxd(j);
     end;
    index=[index;j m];
    if j==2
        index(1,:)=[];
    end;
    clast1=idxd([j m]);
 end;
 kolvo=length(unique(idxd));
 uuu=unique(u);
  figure
 gscatter(X(:,1),X(:,2),idxd', 'bgmkrcy' )
 %gscatter(d1_mahal,d2_mahal,idx, 'bgm' )
 %hold on
 %plot(C(:,1),C(:,2), 'kx' )
 legend( 'Cluster 1' , 'Cluster 2', 'Cluster 3', 'Cluster 4','Cluster 5','Cluster 6' ) 
