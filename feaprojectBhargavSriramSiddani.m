%Reading mesh file and creating nodes and elements data
fid=fopen('Mesh.txt');
line='cc';
no=0;%To track the line number in text file
while(line ~= -1)
    line=fgetl(fid);
    no=no+1;
    if(strcmp(line,'$PhysicalNames'))
       line=fgetl(fid);
       n1=str2double(string(line));
       C1=textscan(fid,'%d %d %s',n1);
       fclose(fid);
        break
    end
end
phy=zeros(n1,3);
phyName=string(zeros(n1,1));
for(i=1:n1)
    for(j=1:n1)
        if(i==C1{2}(j,1))
            phy(i,1)=i; %Physical Name "Number"
            phyName(i,1)=string(C1{3}(j,1)); %Physical Name "Name"
        end
    end
end
clear C1;
fid=fopen('BoundaryConditions.txt');
line='cc';
n5=0; %To track
while(line~=-1)
    line=fgetl(fid);
    n5=n5+1;
    for(i=1:n1)
        line1=string(line);
        if(strcmp(line1,phyName(i,1)))
            line2=char(line);
            line2=string(line2(2:5));
            if(strcmp(line2,"Conv")|strcmp(line2,"Inte"))
                C4=textscan(fid,'%f %f',1);
                phy(i,2)=C4{1}(1,1);
                phy(i,3)=C4{2}(1,1);
                fclose(fid);
                break
            else
                C4=textscan(fid,'%f',1);
                phy(i,2)=C4{1}(1,1);
                fclose(fid);
                break
            end
        end
    end
    fid=fopen('BoundaryConditions.txt');
    for(j=1:n5+1)
        line=fgetl(fid);
    end
    n5=n5+1;
end
clear C4 n5 line1 line2;
no=no+n1+4;%To start from nodes data
fid=fopen('Mesh.txt');
line='cc';
mo=0;
while(line~=-1)
    mo=mo+1;
    line=fgetl(fid);
    if(mo==no)
       NN = str2double(string(line)); %No of NODES
       C2=textscan(fid,'%d %f %f %f',NN);
       fclose(fid);
        break  
    end
end
nodes=zeros(NN,2);
for(i=1:NN)
    nodes(i,1)=C2{2}(i,1);
    nodes(i,2)=C2{3}(i,1);
end
clear C2;
no=mo+NN+3;
mo=0;
fid=fopen('Mesh.txt');
line='cc';
while(line~=-1)
    mo=mo+1;
    line=fgetl(fid);
    if(mo==no)
       NEL = str2double(string(line)); %No of ELEMENTS
       C3=textscan(fid,'%d %d %d %d %d %d %d %d',NEL);
       fclose(fid);
        break  
    end
end
elements=zeros(NEL,5);
for(i=1:NEL)
    elements(i,1)=C3{2}(i,1);
    elements(i,2)=C3{4}(i,1);
    elements(i,3)=C3{6}(i,1);
    elements(i,4)=C3{7}(i,1);
    elements(i,5)=C3{8}(i,1);
end
clear C3;
%Solver part
%[a]{T}={b}
%Applying Gauss Seidel Iterative Scheme
converg=0.000001; %Desired level of convergence
a=zeros(NN,NN);
b=zeros(NN,1);
T=zeros(NN,1);
Tcheck=zeros(NN,1);
n_triangles=0; %To find total number of triangular elements
for(i=1:NEL)
    elm = elements(i,2); %No corresponding to physical name number of element
    name=phyName(elm);
    name=char(name);
    name=string(name(2:5));
    if(strcmp(name,"Inte"))
        n_triangles=n_triangles+1;
        nd1=elements(i,3); %No of node1
        nd2=elements(i,4); %No of node2
        nd3=elements(i,5); %No of node3
        x1=nodes(nd1,1);
        y1=nodes(nd1,2);
        x2=nodes(nd2,1);
        y2=nodes(nd2,2);
        x3=nodes(nd3,1);
        y3=nodes(nd3,2);
        f1=x2*y3-x3*y2;
        f2=x3*y1-x1*y3;
        f3=x1*y2-x2*y1;
        b1=y2-y3;
        b2=y3-y1;
        b3=y1-y2;
        c1=x3-x2;
        c2=x1-x3;
        c3=x2-x1;
        Area=0.5*(x2*y3-x3*y2-x1*y3+x1*y2+x3*y1-x2*y1);%Area of the Triangle
        k=phy(elm,2); %Thermal Diffusivity
        Qe=phy(elm,3); %Heat Generation
        keff = k/(4*Area);
        a(nd1,nd1)=a(nd1,nd1)+keff*(b1*b1+c1*c1);
        a(nd1,nd2)=a(nd1,nd2)+keff*(b1*b2+c1*c2);
        a(nd1,nd3)=a(nd1,nd3)+keff*(b1*b3+c1*c3);
        a(nd2,nd1)=a(nd2,nd1)+keff*(b1*b2+c1*c2);
        a(nd2,nd2)=a(nd2,nd2)+keff*(b2*b2+c2*c2);
        a(nd2,nd3)=a(nd2,nd3)+keff*(b2*b3+c2*c3);
        a(nd3,nd1)=a(nd3,nd1)+keff*(b1*b3+c1*c3);
        a(nd3,nd2)=a(nd3,nd2)+keff*(b2*b3+c2*c3);
        a(nd3,nd3)=a(nd3,nd3)+keff*(b3*b3+c3*c3);
        Qeff=Qe*(Area/3);
        b(nd1,1)=b(nd1,1)+Qeff;
        b(nd2,1)=b(nd2,1)+Qeff;
        b(nd3,1)=b(nd3,1)+Qeff;
    end
    if(strcmp(name,"Flux"))
        q0=phy(elm,2);
        nd1=elements(i,3);
        nd2=elements(i,4);
        x1=nodes(nd1,1);
        y1=nodes(nd1,2);
        x2=nodes(nd2,1);
        y2=nodes(nd2,2);
        Length=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
        Length=sqrt(Length);
        q0eff=0.5*q0*Length;
        b(nd1,1)=b(nd1,1)+q0eff;
        b(nd2,1)=b(nd2,1)+q0eff;
    end
    if(strcmp(name,"Conv"))
        h=phy(elm,2); % h of the element
        Tinf=phy(elm,3); %T ambient
        nd1=elements(i,3);
        nd2=elements(i,4);
        x1=nodes(nd1,1);
        y1=nodes(nd1,2);
        x2=nodes(nd2,1);
        y2=nodes(nd2,2);
        Length=(x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
        Length=sqrt(Length);
        %Contribution to LHS
        h1eff=(h/6)*Length;
        a(nd1,nd1)=a(nd1,nd1)+h1eff*2;
        a(nd1,nd2)=a(nd1,nd2)+h1eff;
        a(nd2,nd1)=a(nd2,nd1)+h1eff;
        a(nd2,nd2)=a(nd2,nd2)+h1eff*2;
        %Contribution to RHS
        h2eff=0.5*h*Tinf*Length;
        b(nd1,1)=b(nd1,1)+h2eff;
        b(nd2,1)=b(nd2,1)+h2eff;
    end
    if(strcmp(name,"Temp"))
        Tem=phy(elm,2);
        nd1=elements(i,3);
        nd2=elements(i,4);
        Tcheck(nd1,1)=1;
        Tcheck(nd2,1)=1;
        T(nd1,1)=Tem;
        T(nd2,1)=Tem;
    end
end
clear elm name nd1 nd2 nd3 Length Area h mo n1;
clear h1eff h2eff q0eff Qeff q0 Tem k keff Qe Tinf;
clear x1 y1 x2 y2 x3 y3;
clear b1 b2 b3 f1 f2 f3 c1 c2 c3;
%Gauss Seidel Iteration
Tcon=10;
while(Tcon>converg)
    Tcon_sum=0;
    for(i=1:NN)
        if(Tcheck(i,1)==0)
            aT=0; %a*T
            for(j=1:NN)
                if(i~=j)
                    aT=aT+a(i,j)*T(j,1);
                end
            end
            T_old=T(i,1); %Previous T(i,1)
            T(i,1)=(1/a(i,i))*(b(i,1)-aT);
            Tcon_sum=Tcon_sum+(T(i,1)-T_old)*(T(i,1)-T_old); 
        end     
    end
    Tcon=sqrt(Tcon_sum/NN);
end
clear Tcon_sum T-old;
x=zeros(3,n_triangles);
y=zeros(3,n_triangles);
c=zeros(3,n_triangles);%To plot Temperatures
Tfreq=zeros(NN,1); %To keep a track of triangular elements connected to a node
fluxTx=zeros(NN,1); %flux of T along x
fluxTy=zeros(NN,1); %flux of T along y
fx=zeros(3,n_triangles); %To plot flux along x
fy=zeros(3,n_triangles); %To plot flux along y
clear n_triangles;
na=0; %Dummy variable to track and update traingle number
for(i=1:NEL)
    elm = elements(i,2); %No corresponding to physical name number of element
    name=phyName(elm);
    name=char(name);
    name=string(name(2:5));
    if(strcmp(name,"Inte"))
        na=na+1;
        nd1=elements(i,3); %No of node1
        nd2=elements(i,4); %No of node2
        nd3=elements(i,5); %No of node3
        x1=nodes(nd1,1);
        y1=nodes(nd1,2);
        x2=nodes(nd2,1);
        y2=nodes(nd2,2);
        x3=nodes(nd3,1);
        y3=nodes(nd3,2);
        T1=T(nd1,1);
        T2=T(nd2,1);
        T3=T(nd3,1);
        x(1,na)=x1;
        x(2,na)=x2;
        x(3,na)=x3;
        y(1,na)=y1;
        y(2,na)=y2;
        y(3,na)=y3;
        c(1,na)=T1;
        c(2,na)=T2;
        c(3,na)=T3;
        b1=y2-y3;
        b2=y3-y1;
        b3=y1-y2;
        c1=x3-x2;
        c2=x1-x3;
        c3=x2-x1;
        Area=0.5*(x2*y3-x3*y2-x1*y3+x1*y2+x3*y1-x2*y1);%Area of the Triangle
        k=phy(elm,2); %Thermal Diffusivity
        Tfreq(nd1,1)=Tfreq(nd1,1)+1;
        Tfreq(nd2,1)=Tfreq(nd2,1)+1;
        Tfreq(nd3,1)=Tfreq(nd3,1)+1;
        fTx=(0.5/Area)*(b1*T1+b2*T2+b3*T3)*(-k);
        fTy=(0.5/Area)*(c1*T1+c2*T2+c3*T3)*(-k);
        %gradTx
        fluxTx(nd1,1)=(fluxTx(nd1,1)*(Tfreq(nd1,1)-1))+fTx;
        fluxTx(nd1,1)=fluxTx(nd1,1)/Tfreq(nd1,1);
        fluxTx(nd2,1)=(fluxTx(nd2,1)*(Tfreq(nd2,1)-1))+fTx;
        fluxTx(nd2,1)=fluxTx(nd2,1)/Tfreq(nd2,1);
        fluxTx(nd3,1)=(fluxTx(nd3,1)*(Tfreq(nd3,1)-1))+fTx;
        fluxTx(nd3,1)=fluxTx(nd3,1)/Tfreq(nd3,1);
        %gradTy
        fluxTy(nd1,1)=(fluxTy(nd1,1)*(Tfreq(nd1,1)-1))+fTy;
        fluxTy(nd1,1)=fluxTy(nd1,1)/Tfreq(nd1,1);
        fluxTy(nd2,1)=(fluxTy(nd2,1)*(Tfreq(nd2,1)-1))+fTy;
        fluxTy(nd2,1)=fluxTy(nd2,1)/Tfreq(nd2,1);
        fluxTy(nd3,1)=(fluxTy(nd3,1)*(Tfreq(nd3,1)-1))+fTy;
        fluxTy(nd3,1)=fluxTy(nd3,1)/Tfreq(nd3,1);
    end
end
clear x1 x2 x3 y1 y2 y3 T1 T2 T3;
clear na nd1 nd2 nd3;
na=0; %Dummy variable to track and update traingle number
for(i=1:NEL)
    elm = elements(i,2); %No corresponding to physical name number of element
    name=phyName(elm);
    name=char(name);
    name=string(name(2:5));
    if(strcmp(name,"Inte"))
        na=na+1;
        nd1=elements(i,3); %No of node1
        nd2=elements(i,4); %No of node2
        nd3=elements(i,5); %No of node3
        x1=nodes(nd1,1);
        y1=nodes(nd1,2);
        x2=nodes(nd2,1);
        y2=nodes(nd2,2);
        x3=nodes(nd3,1);
        y3=nodes(nd3,2);
        x(1,na)=x1;
        x(2,na)=x2;
        x(3,na)=x3;
        y(1,na)=y1;
        y(2,na)=y2;
        y(3,na)=y3;
        fx(1,na)=fluxTx(nd1,1);
        fx(2,na)=fluxTx(nd2,1);
        fx(3,na)=fluxTx(nd3,1);
        fy(1,na)=fluxTy(nd1,1);
        fy(2,na)=fluxTy(nd2,1);
        fy(3,na)=fluxTy(nd3,1);
    end
end
figure
title('Temperature Distribution')
patch(x,y,c,'LineStyle','none')
colorbar
figure
title('Heat Flux along x')
patch(x,y,fx,'LineStyle','none')
colorbar
figure
title('Heat Flux along y')
patch(x,y,fy,'LineStyle','none')
colorbar
clear Area x1 x2 x3 y1 y2 y3;
clear aT b1 b2 b3 c1 c2 c3;
clear nd1 nd2 nd3;
fileid=fopen('Result.txt','wt');
fprintf(fileid,'NodeNo\tx\ty\tT\tqx\tqy\n');
for(i=1:NN)
    fprintf(fileid,'%d\t%f\t%f\t%f\t%f\t%f\n',i,nodes(i,1),nodes(i,2),T(i,1),fluxTx(i,1),fluxTy(i,1));
end
fclose(fileid);