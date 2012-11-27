%% Simulation for peeling_Least degree %%%%
%% Try  to speed up the algorithm
function Cont_Peeling
%  clc;
 close all;
 global p;
 p=0.8;
 global q;
 q=0.4;
 
 
 global m;
 m=200;
 global na;

%% duration of simulation considered
dur=5000;
 
SMpMq=[0 0 0 0];
Pi=1;%% Percent of idle servers
Pp=0; %% Percent of P servers
Pq=0; %% Percent of Q servers

  %% arrival rate.  
 lambda1=[0.01 0.05 0.1 0.2];%[q-0.1:0.05:p-0.1-0.025];
  lambda2=[p-0.1:0.01:p-0.025];
  lambda3=[0.6:0.05:0.75];
  lambda=[lambda1 lambda2 lambda3];



Nlambda=length(lambda);
QsizeP=zeros(1,Nlambda);
dQP=zeros(1,Nlambda);
MpP=zeros(1,Nlambda);
dMpP=zeros(1,Nlambda);
MqP=zeros(1,Nlambda);
dMqP=zeros(1,Nlambda);
 
lamI=1;
%% for each arrival 
 while(lamI<=Nlambda)
     %%totoal arrival for the system
     ar=m*lambda(lamI);
     
     LAM=lambda(lamI)
     
     tQsize=0;
     tdQ=0;
     tMp=0;
     tdMp=0;
     tMq=0;
     tdMq=0;
     
     maxdat=ceil(dur*ar);   

 inqt=0; %% initial number of queueing tasks
 
 Rnqt=[];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nqt=inqt; %% # of queueing tasks
 %% local data servers index for queueing tasks %%
 QueTL=ceil(rand(inqt,3)*m);
  
 %% initial server states %%%
 SS=zeros(m,1); %% server state: p-sever or q-sever
 SD=zeros(m,1); %% server degree
 SJ=zeros(m,1);%% remaining service time for job processing at server
 
 nps=floor(m*Pp);%% number of p-sever
 nqs=floor(m*Pq);
 nis=m-nps-nqs;
 
 Randm=randperm(m); %%
 
%% NumberUnassignedT = 0; %%initialize the number of fresh tasks left in 
 %% the last time interval
 
 %%% Indicate P/Q servers %%%
SP=Randm(1:nps);
SQ=Randm(nps+1:nps+nqs);
SS(SP)=-1; %% P-server
SS(SQ)=-2; %% Q-server

 %% each arrival task
for ic=1:maxdat
%           ic
         interval=exprnd(1/ar);
         intvl=interval;
         %%nnewt = 1;
         nnewt=1;
        %% NumberUnassignedT = 0;     
         %% departure of processing tasks at servers        
         while(intvl>0)
             BusyS=find(SS~=0);
             [sm smid]=min(SJ(BusyS));    
             if(sm>intvl)
                  %% no leave during interval
                 SJ(BusyS)=SJ(BusyS)-intvl;
                 intvl=0;
             else
                 %% departure during interval
                 DepSID=BusyS(smid);
                 SS(DepSID)=0;
                 SJ(DepSID)=0;
                 StBusyS=setdiff(BusyS,DepSID);
                 SJ(StBusyS)=SJ(StBusyS)-sm;
                 intvl=intvl-sm;
                 
                 if(nqt>0)  %% nqt: number of queueing task
                    %% tasks queueing, so no other idle servers available
                    [jr,jc,v]=find(QueTL==DepSID);
                    
                    if(isempty(jr)) 
                        %%% no local data task randomly selected job to this server
                        trandp=randperm(nqt);
                        AssTID=trandp(1);
                        SJ(SID)=exprnd(1/q);
                        SS(SID)=-2;
                    else
                        %% Task with local data is assigned
                        LocalT=unique(jr);
                        trandp=randperm(length(LocalT));
                        AssTIDs=trandp(1);
                        AssTID=LocalT(AssTIDs);
                        SJ(SID)=exprnd(1/p);
                        SS(SID)=-1;
                        
                    end
                    
                    SD(QueTL(AssTID,:))=SD(QueTL(AssTID,:))-1;%% update server degree.
                    
                    TempQT=[QueTL(1:AssTID-1,:)
                            QueTL(AssTID+1:nqt,:)];

                   
                        
                       %% NumberUnassignedT = length(TempQT); 
                        %% render the number of unassigned tasks to fresh nodes
                        
            
                       nqt=nqt-1;
                       QueTL=TempQT;
                 end
             end
         end
       
         Rnqt=[Rnqt nqt];%% record number of queueing tasks
         tQsize=tQsize+nqt;
         tdQ=tdQ+nqt*nqt;

        %% Statics of P/Q/I servers and queueing tasks
         %%% P servers %%%
         SP=find(SS==-1);
         fp=length(SP)/m;
         tMp=tMp+fp;
         tdMp=tdMp+fp*fp;
         %%% Q servers %%%
         SQ=find(SS==-2);
         fq=length(SQ)/m;
         tMq=tMq+fq;
         tdMq=tdMq+fq*fq;         
         
        %% Assignment of server for new arrival
         NewTL=zeros(nnewt,3);
         
         for newarrial=1:nnewt
             trandp=randperm(m);
            NewTL(newarrial,:)=trandp(1:3);
         end
         
        SI=find(SS==0);
        %%% update degree of servers 
        X=[0:1:m];  
        dds1=hist(NewTL(:,1),X);
        dds2=hist(NewTL(:,2),X);
        dds3=hist(NewTL(:,3),X);
        dds=dds1(2:m+1)+dds2(2:m+1)+dds3(2:m+1);%% degrees of all servers
        SD=SD+dds'; %% degrees of servers updates       
        
        
        nqt=nqt+nnewt;
        QueTL=[QueTL
                   NewTL];   
               
        if(isempty(SI)==0) %% Some servers are idle
            
            %% send freshment to QueTL, nqt has been incredmented with 'nnewt'
           QueTL=ceil(rand(nqt,3)*m);
 
             nIdleS=length(SI);
            
            if(nqt<nIdleS) 
                NassT=nqt;
            else
                NassT=nIdleS;
            end
            %% Improved Peeling: start from servers of least degree
%                 SelecS=randperm(nIdleS,NassT);
                         
             
                tnqt=nqt;
                for ass=1:NassT
                    SI=find(SS==0);
                    nIdleS=length(SI);
                    dois=SD(SI); %% degree of idle servers
                    Non0DOISind=find(dois); 
                    smdeg=min(dois(Non0DOISind)); %%find the smallest degree (non-zero)
                    
                    
                    if(isempty(smdeg))
                        trandp=randperm(nIdleS);
                        SelecS=trandp(1);                
                        SID=SI(SelecS);
                    else
                        smdegSind=find(dois==smdeg);
                        candidateSID=SI(smdegSind);
                        trandp=randperm(length(candidateSID));
                        SelecS=trandp(1);
                        SID=candidateSID(SelecS);
                    end
                    
                    [jr,jc,v]=find(QueTL==SID);
                    
                    if(isempty(jr)) 
                        %%% no local data task randomly selected job to this server
                        trandp=randperm(tnqt);
                        AssTID=trandp(1);
                        SJ(SID)=exprnd(1/q);
                        SS(SID)=-2;
                    else
                        %% Task with local data is assigned
                        LocalT=unique(jr);
                        trandp=randperm(length(LocalT));
                        AssTIDs=trandp(1);
                        AssTID=LocalT(AssTIDs);
                        SJ(SID)=exprnd(1/p);
                        SS(SID)=-1;
                        
                    end
                        SD(QueTL(AssTID,:))=SD(QueTL(AssTID,:))-1;
                        
                        %%take the assigned task out of the queue
                        TempQT=[QueTL(1:AssTID-1,:)
                            QueTL(AssTID+1:tnqt,:)];
                        
                        tnqt=tnqt-1; %% decrease the number of tempQ
                        
                    %%    NumberUnassignedT = length(TempQT); 
                        %% render the number of unassigned tasks to fresh tasks
                        
                     %%   tnqt = tnqt - 1; %%reset the number of tempQ
                        
                        
                        QueTL=TempQT; %% update the Q with the tempQ
                end
                
                nqt=nqt-NassT;       
         end
        
        end
        QsizeP(lamI)=tQsize/maxdat;
        dQP(lamI)=tdQ/maxdat-QsizeP(lamI)^2;%% variance
        MpP(lamI)=tMp/maxdat;
        dMpP(lamI)=tdMp/maxdat-MpP(lamI)^2;
        MqP(lamI)=tMq/maxdat;
        dMqP(lamI)=tdMq/maxdat-MqP(lamI)^2;
            
        lamI=lamI+1;
        
 save(['randservernewr',num2str(LAM),'m',num2str(m),'dur',num2str(dur),'p',num2str(p),'q',num2str(q),'.mat'],'lambda','QsizeP','dQP','MpP','dMpP','MqP','dMqP','Rnqt');
  
 end
 

end


%% d.d for idle servers %%
function D=RI(NL,S)
global Rm;
global m;
  SI=find(S==0);
  D=zeros(1,Rm+1);
  if(~isempty(NL))
      for i=1:length(SI)
          ii=find(NL==SI(i));
          if(length(ii)<Rm+1)
              D(length(ii)+1)=D(length(ii)+1)+1;
          end
      end
      D=D/m;
  else
    D(1)=length(SI)/m;
  end
  
    
end

%% d.d for P servers %%
function D=RP(NL,S)
global Rm;
global m;
  SP=find(S==-1);
  D=zeros(1,Rm+1);
   if(~isempty(NL))
       for i=1:length(SP) 
           ip=find(NL==SP(i));
           if(length(ip)<Rm+1)
               D(length(ip)+1)=D(length(ip)+1)+1;
           end
       end
   else
       D(1)=length(SP);
   end
  
  D=D/m;
  
end

%% d.d for Q servers %%
function D=RQ(NL,S)
 global Rm;
 global m;
  SQ=find(S==-2);
  D=zeros(1,Rm+1);
  if(~isempty(NL))
      for i=1:length(SQ)
          iq=find(NL==SQ(i));
          if(length(iq)<Rm+1)
              D(length(iq)+1)=D(length(iq)+1)+1;
          end
      end
   else
     D(1)=length(SQ);
  end
  
  D=D/m; 
  
end