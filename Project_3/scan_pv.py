filenameA = 'pv-FS.d';
filenameB = 'pv-FD.d';

fidA  = fopen(filenameA,'r');
fidB  = fopen(filenameB,'r');

np = 1000;
dt = 1e-3;
sf = 1000;
for i=1:np
    sA = get_frame(fidA);
    sB = get_frame(fidB);
    
    v = sA.x(:) -sB.x(:);
    d(i) = sqrt(dot(v,v));
    a(i) = dot(sA.v(:),sB.v(:))/sqrt(dot(sA.v(:),sA.v(:))*dot(sB.v(:),sB.v(:)));
end
fclose(fidA);
fclose(fidB);

np =50;
a = acos(a)*180/pi;

figure(1), 
hd1 = plot((1:np)*sf*dt,log(d(1:np)));
set(hd1,'LineWidth',1.5);
ha1 = gca;
xlabel('Time (ps)','FontSize',24)
ylabel('Log[Phasespace seperation (Ang)]','FontSize',24)
set(ha1,'LineWidth',1)
set(ha1,'Box','On')
set(ha1,'FontSize',24)
%set(ha1,'ytick',[])


figure(2), 
hd2 = plot((1:np)*sf*dt,a(1:np));
set(hd2,'LineWidth',1.5);
ha2 = gca;
xlabel('Time (ps)','FontSize',24)
ylabel('Angle between trajectories (deg)','FontSize',24)
set(ha2,'LineWidth',1)
set(ha2,'Box','On')
set(ha2,'FontSize',24)




function round_data(infile)

addpath /nfs/mohr/u1/MATS588/Shared/Tools/Matlab-LMP-Functions
%addpath /Users/Alex/Class/Shared/Tools/Matlab-LMP-Functions

%infile  = 'data.Si-equilibriated-300K';
outfile = [infile,'-single-prescission'];

[h,b] = load_datafile(infile)
write_data_file(outfile,h,b)
