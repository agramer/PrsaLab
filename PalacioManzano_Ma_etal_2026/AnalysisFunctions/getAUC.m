function [AUC, Px, Py] = getAUC(x,y,Nx,Ny)
    
    % y is the positive class

    My=mean(y,1);
    DVx=zeros(Nx,1);
    for t=1:Nx
        Mx=mean(x(setdiff(1:Nx,t),:),1);
        DVx(t,1)=dot(x(t,:),Mx-My);
    end


    Mx=mean(x,1);
    DVy=zeros(Ny,1);
    for t=1:Ny
        My=mean(y(setdiff(1:Ny,t),:),1);
        DVy(t,1)=dot(y(t,:),Mx-My);
    end
    
    Clim=[min([DVx; DVy]) max([DVx; DVy])];
    Res=(Clim(2)-Clim(1))/1000;
    Pcnt=0;
    crit=Clim(1)-Res:Res:Clim(2)+Res;
    
    crit=sort([DVx' DVy']);
    
    Px=zeros(size(crit));
    Py=zeros(size(crit));
    
    for k=crit
        Pcnt=Pcnt+1;
        Px(1,Pcnt)=1-sum(DVx>k)/Nx;
        Py(1,Pcnt)=1-sum(DVy>k)/Ny;
    end
    

    AUC=trapz(Px,Py);


end