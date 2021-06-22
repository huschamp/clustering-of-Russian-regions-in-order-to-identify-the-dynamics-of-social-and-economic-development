clear;
clc;
regions=caseread('C:\big\four\regionst.txt');
p1=dlmread('C:\big\four\doh.txt');
p2=dlmread('C:\big\four\ves.txt');
p3=dlmread('C:\big\four\prir.txt');
p4=dlmread('C:\big\four\bezr.txt');
p5=dlmread('C:\big\four\inter.txt');
table=[p1 p2 p3 p4 p5];
X=table;
X(24,:)=[];
X(20,:)=[];
X(18,:)=[];
X(10,:)=[];
regions(24,:)=[];
regions(20,:)=[];
regions(18,:)=[];
regions(10,:)=[];
%%%X(27,:)=[];
%%%X(24,:)=[];
%X(21,:)=[];
%X(20,:)=[];
%X(19,:)=[];
%X(18,:)=[];
Xc=X;
table=X;
X1=X;
Coef=[];
 for j=1:length(X);
     for i=1:5
        minn=min(X(j,:));
        maxx=max(X(j,:));
        sr=sum(X(j,:))/5;
    %    X1(j,i)=(X(j,i)-sr)/(maxx-minn);
      % X1(j,i)=(X(j,i))/(maxx);
      % X1(j,i)=(X1(j,i)-minn)/(maxx-minn);
     end;
 end;
 for i=1:5
     maxx=max(X(:,i));
     X(:,i)=X(:,i)/maxx;
 end;
%X=X1;
% X=Xc;
% % %  ma=max(max(X));
% % %  mi=min(min(X));
% % %  Xnorm=norm(X,1);
% % %  X(:,:)=X/Xnorm;
 X1=X;
% % %  X1(:,2)=X1(:,2)*-10000;
% % %  X1(:,3)=X1(:,3)*10000;
% % %  X1(:,4)=X1(:,4)*1000;
% % % X1(:,5)=X1(:,5)*100;
X=X1;
 XXX=[];
 for i=1:length(X)
 XXX=[XXX sum(X1(i,:))];
 end;
 XXX=XXX';
 [idkm,C1] = kmeans(X,5);
%for i=1:5
%    for j=1:5
   figure
gscatter(X(:,3),X(:,2),idkm, 'rkbgm' )%%1 4!! %% 1 3 norm %% 3%11 %% 4 16 17 18
%    end;
%end
hold on;
plot(C1(:,3),C1(:,2), 'kx' )
legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3','Cluster 4','Cluster 5' , 'Cluster Centroid' )
for cl=1:length(X)
[idkm,C1] = kmeans(X,cl); 
%figure
%gscatter(X(:,4),X(:,5),idx, 'rkb' )
%hold on
%plot(C(:,1),C(:,2), 'kx' )
%legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3' , 'Cluster Centroid' ) 
%%%%%%%%%%%%%%%%
%X=table;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Zkm = linkage(X, 'single', 'mahalanobis');
%Zkm = linkage(X, 'average');
%cutoff = median([Z(end-2,3) Z(end-1,3)]);
%dendrogram(Z,'ColorThreshold',cutoff)
%X=X1;

%X=X1;
%%%%pb=squareform(pdist(X,'mahalanobis'));
psr=squareform(pdist(X,'mahalanobis'));
%figure
Zsr = linkage(psr, 'average');
Tsr = cluster(Zsr,'maxclust',cl);
%X=X1;
%gscatter(X(:,2),X(:,3),T, 'rkb' )
%hold on
%legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3')
%cutoff = median([Zsr(end-2,3) Zsr(end-1,3)]);
%dendrogram(Zsr,0,'ColorThreshold',cutoff)
%dendrogram(Zsr,0)
%%%%%X=table;
%%%figure
%%%H2=dendrogram(Zsr,0)
%%figure 
%%Tsrs = cluster(Zsr,'maxclust',3);
%%gscatter(X(:,2),X(:,3),Tsrs, 'rkb' )

pb=squareform(pdist(X,'mahalanobis'));
%%%figure
Zb = linkage(pb, 'single');
%%%%H1=dendrogram(Zb,0)
%ZZ2 = linkage(X, 'single');
Tb = cluster(Zb,'maxclust',cl);
%T = clusterdata(X, 'linkage' , 'single', 'mahalanobis','maxclust',6);
%gscatter(X(:,2),X(:,3),T, 'bgm' )
%hold on
%legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3')
%cutoff = median([Zb(end-2,3) Zb(end-1,3)]);
%dendrogram(Zb,'ColorThreshold',cutoff)
%dendrogram(Zb,0)

%figure
pp=squareform(pdist(X,'mahalanobis'));
Zp = linkage(pp, 'complete');
Tp = cluster(Zp,'maxclust',cl);
% % % % figure
%%%%H3=dendrogram(Zp,0)
%for i=1:5
 %     for j=1:5
%         gscatter(X(:,2),X(:,3),T, 'rkb' )
   %   end;
 % end;
%hold on
%legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3')
%cutoff = median([ZZ(end-2,3) ZZ(end-1,3)]);
%dendrogram(ZZ,'ColorThreshold',cutoff)
%%%%%%%%%%%%%%%TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
%cl=5;
%ZZ = linkage(pb, 'complete');
%TU = cluster(ZZ,'maxclust',cl);
for metod=1:4
%%%%
if metod==1
    TU=idkm;
end;
if metod==2
    TU=Tsr;
end;
if metod==3
    TU=Tb;
end;
if metod==4
    TU=Tp;
end;
BBs=[];
srvb=[];
for i=1:length(X)
    srvb=[srvb sum(X(i,:))];
end;
srvb=sum(srvb)/length(X);
for i=1:cl
    n=find(TU==i);
    sr=[];
        for j=1:length(n)
            sr=[sr sum(X(n(j),:))];
        end;
        sr=sum(sr);
        sr=sr/length(n);
        s=[sr;srvb];
    BBs=[BBs length(n)*(pdist(s))^2]%(sr-srvb)^2)];
end;
BB=sum(BBs);
%B=sum(B);
L=[];
for i=1:length(X)
    ss=[sum(X(i,:));srvb];
    L=[L pdist(ss)^2];
end;
V=sum(L);
B=BB;
W=V-B;
%t=1-((V-B)/V)=1-(1-B/V)=B/V
T=1-(W/V);
%TT=B/V;
if cl>1 && metod==1
   %Coef=[Coef; T]; 
end;
Coef=[Coef T];

end;

end;
Coef1=[];
Coef2=[];
Coef3=[];
Coef4=[];
i=1;
while i<=length(Coef)
    if i==1
    Coef1=[Coef1 Coef(i)];
    i=i+4;
    else
    Coef1=[Coef1 Coef(i)];
    i=i+4;
    end;
end;
i=2;
while i<=length(Coef)
    if i==1
    Coef2=[Coef2 Coef(i)];
    i=i+4;
    else
    Coef2=[Coef2 Coef(i)];
    i=i+4;
    end;
end;
i=3;
while i<=length(Coef)
    if i==1
    Coef3=[Coef3 Coef(i)];
    i=i+4;
    else
    Coef3=[Coef3 Coef(i)];
    i=i+4;
    end;
end;
i=4;
while i<=length(Coef)
    if i==1
    Coef4=[Coef4 Coef(i)];
    i=i+4;
    else
    Coef4=[Coef4 Coef(i)];
    i=i+4;
    end;
end;
Coef=[Coef1;Coef2;Coef3;Coef4];
Coef=Coef';
for i=1:length(X)
    for j=1:4
        if Coef(j,i)<0.8
            Coef(j,i)=2;
        end
    end
end
Coef=Coef';
[A1 id1]=min(Coef(:,1));
[A2 id2]=min(Coef(:,2));
[A3 id3]=min(Coef(:,3));
[A4 id4]=min(Coef(:,4));
id=[id1 id2 id3 id4];
A=[A1 A2 A3 A4];
id=min(min(id1,id2),min(id3,id4));
cl=id;
figure
psr=squareform(pdist(X,'mahalanobis'));
Zsr = linkage(psr, 'average');
Tsr = cluster(Zsr,'maxclust',cl);
zzz1=Zsr(end-1,3)
zzz2=Zsr(end-2,3)
zzz3=Zsr(end-3,3)
zzz4=Zsr(end-4,3)
cutoff = median([ Zsr(end-3,3) Zsr(end-2,3) Zsr(end-1,3)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dendrogram(Zsr,0)%,'ColorThreshold',cutoff)
% figure
%  gscatter(X(:,1),X(:,4),Tsr, 'rkbgm' ) %%%%% 4;5
% figure
% pb=squareform(pdist(X,'mahalanobis'));
% Zb = linkage(pb, 'single');
% Tb = cluster(Zb,'maxclust',cl);
% cutoff = median([Zb(end-2,3) Zb(end-1,3)]);
% dendrogram(Zb,0)%,'ColorThreshold',cutoff)
% figure
% pp=squareform(pdist(X,'mahalanobis'));
% Zp = linkage(pp, 'complete');
% Tp = cluster(Zp,'maxclust',cl);
% cutoff = median([Zp(end-2,3) Zp(end-1,3)]);
% dendrogram(Zp,0)%,'ColorThreshold',cutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k=0;
% if id(1)~=min(id);
%     k=1;
% end;
% for i=4:-1:2
%     if id(i)>id(i-1)
%         id(i)=-1;
%     end;
% end;
% if k==1
%     id(1)=-1;
% end;
% AI=[A;id];
% AI=sum(AI());
% [S idx]=max(AI);
% cl=max(id);
cl=id;
if id==id1
    cl=4;
    [idkm,C] = kmeans(X,cl); 
    figure
   % for i=1:5
   %     for j=1:5
   %         figure
         gscatter(X(:,3),X(:,2),idkm, 'rkbgm' ) %%%%% 4;5
   %     end;
   % end;
    hold on
    plot(C(:,3),C(:,2), 'kx' )
    legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3', 'Cluster 4', 'Cluster Centroid' ) 
    %%%%%%%%%%%%%%%%
    %X=table;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Zkm = linkage(X, 'single', 'mahalanobis');
  %  Z = linkage(X, 'single','mahalanobis');
  %  cutoff = median([Z(end-2,3) Z(end-1,3)]);
  %  dendrogram(Z,'ColorThreshold',cutoff)
end;
if id==id2
    psr=squareform(pdist(X,'mahalanobis'));
    figure
    ZZ = linkage(psr, 'average');
    T = cluster(ZZ,'maxclust',cl);
    C=[];
    for i=1:cl
       Cc=[];
       Cc=[Cc;X(T==i,:)];
       if i==1
       C=[C sum(Cc)/length(Cc)];
       else
        if length(sum(Cc))==1
           C=[C; Cc]
           continue;
        end;
        C=[C;sum(Cc)/length(Cc)];   
       end;
    end;
    X=X1;
    gscatter(X(:,1),X(:,4),T, 'rkbgm' )
    hold on
    plot(C(:,1),C(:,4), 'kx' )
    legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3','Cluster 4', 'Cluster Centroid' )
    cutoff = median([ZZ(end-2,3) ZZ(end-1,3)]);
    dendrogram(ZZ,'ColorThreshold',cutoff)
end;
if id==id3
     psr=squareform(pdist(X,'mahalanobis'));
    figure
    ZZ = linkage(psr, 'single');
    T = cluster(ZZ,'maxclust',cl);
    C=[];
    for i=1:cl
       Cc=[];
       Cc=[Cc;X(T==i,:)];
       if i==1
           if length(sum(Cc))==1
           C=[C; Cc]
           continue;
           end;
       C=[C sum(Cc)/length(Cc)];
       else
        if length(sum(Cc))==1
           C=[C; Cc]
           continue;
        end;
        C=[C;sum(Cc)/length(Cc)];   
       end;
    end;
    X=X1;
    gscatter(X(:,1),X(:,4),T, 'rkbgm' )
    hold on
    plot(C(:,1),C(:,4), 'kx' )
    legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3','Cluster 4', 'Cluster Centroid' )
    cutoff = median([ZZ(end-2,3) ZZ(end-1,3)]);
    dendrogram(ZZ,'ColorThreshold',cutoff)
end;
if id==id4
    psr=squareform(pdist(X,'mahalanobis'));
    figure
    ZZ = linkage(psr, 'complete');
    T = cluster(ZZ,'maxclust',cl);
     C=[];
    for i=1:cl
       Cc=[];
       Cc=[Cc;X(T==i,:)];
       if i==1
       C=[C sum(Cc)/length(Cc)];
       else
        if length(sum(Cc))==1
           C=[C; Cc]
           continue;
        end;
        C=[C;sum(Cc)/length(Cc)];   
       end;
    end;
    X=X1;
    gscatter(X(:,1),X(:,4),T, 'rkbgm' )
    hold on
    plot(C(:,1),C(:,4), 'kx' )
    legend( 'Cluster 1' , 'Cluster 2' , 'Cluster 3','Cluster 4', 'Cluster Centroid' )
    cutoff = median([ZZ(end-2,3) ZZ(end-1,3)]);
    dendrogram(ZZ,'ColorThreshold',cutoff)
end;
XCV=[regions X];