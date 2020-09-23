function self = bdfscan(filename,plotfig)
% Open example blade definition file and plots geometry
%   filename = name of file to read in
%   plotfig = set to plot blade and streamsurfaces
fid=fopen(filename,'r');
for i=1:3
    line=fgets(fid);
end
nstubes=fscanf(fid,'%g',[1]);
self.nblades=fscanf(fid,'%g',[1]); % number of blades (or periodicity)
flipfoil=fscanf(fid,'%g',[1]);
line =fgets(fid);
line=fgets(fid);
le_tube_index = zeros(nstubes,1);
te_tube_index = zeros(nstubes,1);
alpha1 = zeros(nstubes,1);
alpha2 = zeros(nstubes,1);
rle = zeros(nstubes,1);
rte = zeros(nstubes,1);
te_index = zeros(nstubes,1);
points = zeros(nstubes,1);
for i=1:nstubes
    tubepoints=fscanf(fid,'%g',[1]);
    something=fscanf(fid,'%g',[2]);
    le_tube_index(i) = something(1);
    te_tube_index(i) = something(2);
    if te_tube_index(i)>tubepoints
        le_tube_index(i) = 1;
        te_tube_index(i) = tubepoints;
    end
    points(i)=fscanf(fid,'%g',[1]);
    te_index(i)=fscanf(fid,'%g',[1]);
    un=fscanf(fid,'%g',[2]);
    rle(i)=fscanf(fid,'%g',[1]);
    rte(i)=fscanf(fid,'%g',[1]);
    alpha1(i)=fscanf(fid,'%g',[1]);
    alpha2(i)=fscanf(fid,'%g',[1]);
    ztube{i} = zeros(tubepoints,1);
    rtube{i} = zeros(tubepoints,1);
    ytube{i} = zeros(tubepoints,1);
    for j=1:tubepoints
        ztube{i}(j)=fscanf(fid,'%g',[1]);
        rtube{i}(j)=fscanf(fid,'%g',[1]);
        ytube{i}(j)=0;
    end
    for j=1:points(i)
        r{i}(j)=fscanf(fid,'%g',[1]);
        theta{i}(j)=fscanf(fid,'%g',[1]);
        z{i}(j)=fscanf(fid,'%g',[1]);
    end
    line=fgets(fid);
    line=fgets(fid);
end
for i=1:nstubes
    if flipfoil==1
        theta{i} = -theta{i};
    end
    xx{i}=r{i}.*cos(theta{i});
    yy{i}=-r{i}.*sin(theta{i});
end
if exist('plotfig','var')
    xplot = z;
    yplot = yy;
    zplot = xx;
    figure(1)
    clf
    hold on    
    for i=1:nstubes
        plot3(xplot{i},yplot{i},zplot{i},'k')
        plot3(ztube{i},ytube{i},rtube{i},'r')
    end
    view(0,0);
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
end
self.r = r; % blade radial coordinate (a vector for each section inside a cell format)
self.z = z; % blade axial coordinate
self.t = theta; % blade angle coordinate
self.ztube = ztube; % streamsurface axial coordinate
self.rtube = rtube; % streamsurface radial coordinate
self.te_index = te_index; % index for the trailing edge of the section
self.alpha_in = alpha1; % inlet metal angle
self.alpha_out = alpha2; % exit metal angle
self.rle = rle;
self.rte = rte;
self.le_tube_index = le_tube_index; % index for the leading edge in the streamsurface
                                    % vector for the section
self.te_tube_index = te_tube_index; % index for the trailing edge in the streamsurface
                                    % vector for the section
self.points = points; % number of points defining each section
self.hub = [self.ztube{1},self.rtube{1}];
self.casing = [self.ztube{nstubes},self.rtube{nstubes}];
end

