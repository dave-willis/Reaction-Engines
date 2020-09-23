function write_geomTurbo(filename,hub,casing,blades,varargin)
%write_geomTurbo Writes out input file for Autogrid meshing.
%   Detailed explanation goes here
nargin = length(varargin);
even=0;
sqte = 0;
for n=1:nargin
    if strcmpi(varargin{n},'sqte') %if the entry varargin{n} is a string 'sqte'
        sqte =1; % square trailing edge
    elseif iscell(varargin{n}) %if varargin{n} is a cell array
        even = 1;
    elseif isstruct(varargin{n})
        shrouds = varargin{n}.shrouds;
        scale = varargin{n}.scale;
    elseif length(varargin{n})==1 
        scale = varargin{n};
    else
        splitter = varargin{n};
    end
end
% ~ negates the logical 1
for n=1:length(blades)
    if ~isfield(blades{n},'type') % if 'type' is not a field in the structured array blades
        blades{n}.type = 'normal';
    end
    if ~isfield(blades{n},'rpm')
        blades{n}.rpm = 0;
    end
end

if ~exist('scale','var') % if 'scale' does not exist as a variable
    scale = 1;
end
fn = [filename '.geomTurbo']; % make the file a .geomTurbo
fid = fopen(fn,'w'); % write a new file
fprintf(fid,'GEOMETRY TURBO\n'); % \n starts a new line
fprintf(fid,'VERSION\t5.6\n'); % \t writes a horiziontal tab
fprintf(fid,'TOLERANCE	 1e-006\n');
fprintf(fid,'UNITS\t Meters\n');
fprintf(fid,'UNITS-FACTOR\t 1\n');
if exist('splitter','var')
    fprintf(fid,'byPass yes\n');
else
    fprintf(fid,'byPass no\n');
end
fprintf(fid,'NI_BEGIN CHANNEL\n');
fprintf(fid,' NI_BEGIN basic_curve\n');
fprintf(fid,'  NAME\thub\n');
fprintf(fid,'  DISCRETISATION 10\n');
fprintf(fid,'  DATA_REDUCTION 0\n');
fprintf(fid,'  NI_BEGIN zrcurve\n');
fprintf(fid,'  ZR\n');
zrhub = scale.*hub;
zrlen = length(zrhub(:,1)); % number of rows in zrhub
fprintf(fid,'   %i\n',zrlen); % %i denotes base 10 value
for i=1:zrlen
    fprintf(fid,'   %0#6.12f  %0#6.12f\n',... % specifying precision of values
        [zrhub(i,1),zrhub(i,2)]);
end 
fprintf(fid,'  NI_END zrcurve\n');
fprintf(fid,' NI_END   basic_curve\n');
fprintf(fid,' NI_BEGIN basic_curve\n');
fprintf(fid,'  NAME\tshroud\n');
fprintf(fid,'  DISCRETISATION 10\n');
fprintf(fid,'  DATA_REDUCTION 0\n');
fprintf(fid,'  NI_BEGIN zrcurve\n');
fprintf(fid,'   ZR\n');
zrcas = scale.*casing;
zrlen = length(zrcas(:,1));
fprintf(fid,'   %i\n',zrlen);
for i=1:zrlen
    fprintf(fid,'   %0#6.12f  %0#6.12f\n',...
        [zrcas(i,1),zrcas(i,2)]);
end
fprintf(fid,'  NI_END zrcurve\n');
fprintf(fid,' NI_END   basic_curve\n');
if exist('shrouds', 'var')
    for nrow=1:length(blades);
        fprintf(fid,' NI_BEGIN basic_curve\n');
        fprintf(fid,'  NAME\tsolid\n');
        fprintf(fid,'  DISCRETISATION 10\n');
        fprintf(fid,'  DATA_REDUCTION 0\n');
        fprintf(fid,'  NI_BEGIN zrcurve\n');
        fprintf(fid,'   ZR polyline\n');
        zrbladeshr = scale.*shrouds.shroud{nrow}.blade;
        zrlen = length(zrbladeshr(:,1));
        fprintf(fid,'   %i\n',zrlen);
        for i=1:zrlen
            fprintf(fid,'   %0#6.12f  %0#6.12f\n',...
                [zrbladeshr(i,1),zrbladeshr(i,2)]);
        end
        fprintf(fid,'  NI_END zrcurve\n');
        fprintf(fid,' NI_END   basic_curve\n');
        fprintf(fid,' NI_BEGIN basic_curve\n');
        fprintf(fid,'  NAME\tsolid\n');
        fprintf(fid,'  DISCRETISATION 10\n');
        fprintf(fid,'  DATA_REDUCTION 0\n');
        fprintf(fid,'  NI_BEGIN zrcurve\n');
        fprintf(fid,'   ZR polyline\n');
        zrbladeshr = scale.*shrouds.shroud{nrow}.cas;
        zrlen = length(zrbladeshr(:,1));
        fprintf(fid,'   %i\n',zrlen);
        for i=1:zrlen
            fprintf(fid,'   %0#6.12f  %0#6.12f\n',...
                [zrbladeshr(i,1),zrbladeshr(i,2)]);
        end
        fprintf(fid,'  NI_END zrcurve\n');
        fprintf(fid,' NI_END   basic_curve\n');
        fprintf(fid,' NI_BEGIN basic_curve\n');
        fprintf(fid,'  NAME\trotor/stator\n');
        fprintf(fid,'  DISCRETISATION 10\n');
        fprintf(fid,'  DATA_REDUCTION 0\n');
        fprintf(fid,'  NI_BEGIN zrcurve\n');
        fprintf(fid,'   ZR polyline\n');
        zrbladeshr = scale.*shrouds.shroud{nrow}.sep;
        zrlen = length(zrbladeshr(:,1));
        fprintf(fid,'   %i\n',zrlen);
        for i=1:zrlen
            fprintf(fid,'   %0#6.12f  %0#6.12f\n',...
                [zrbladeshr(i,1),zrbladeshr(i,2)]);
        end
        fprintf(fid,'  NI_END zrcurve\n');
        fprintf(fid,' NI_END   basic_curve\n');
    end
end 
if exist('splitter','var')
    fprintf(fid,' NI_BEGIN basic_curve\n');
    fprintf(fid,'  NAME\tnozzle\n');
    fprintf(fid,'  DISCRETISATION 10\n');
    fprintf(fid,'  DATA_REDUCTION 0\n');
    fprintf(fid,'  NI_BEGIN zrcurve\n');
    fprintf(fid,'   ZR\n');
    zrcas = scale.*splitter;
    zrlen = length(zrcas(:,1));
    fprintf(fid,'   %i\n',zrlen);
    for i=1:zrlen
        fprintf(fid,'   %0#6.12f  %0#6.12f\n',...
            [zrcas(i,1),zrcas(i,2)]);
    end
    fprintf(fid,'  NI_END zrcurve\n');
    fprintf(fid,' NI_END   basic_curve\n');
end
fprintf(fid,' NI_BEGIN channel_curve hub\n');
fprintf(fid,'  NAME		 hub\n');
fprintf(fid,'  VERTEX	 CURVE_P hub 0\n');
fprintf(fid,'  VERTEX	 CURVE_P hub 1\n');
fprintf(fid,' NI_END   channel_curve hub\n');
fprintf(fid,' NI_BEGIN channel_curve shroud\n');
fprintf(fid,'  NAME		 shroud\n');
fprintf(fid,'  VERTEX	 CURVE_P shroud 0\n');
fprintf(fid,'  VERTEX	 CURVE_P shroud 1\n');
fprintf(fid,' NI_END   channel_curve shroud\n');
if exist('splitter','var')
    fprintf(fid,' NI_BEGIN channel_curve nozzle\n');
    fprintf(fid,'  NAME		 nozzle\n');
    fprintf(fid,'  VERTEX	 CURVE_P nozzle 0\n');
    fprintf(fid,'  VERTEX	 CURVE_P nozzle 1\n');
    fprintf(fid,' NI_END   channel_curve nozzle\n');
end
fprintf(fid,'NI_END   CHANNEL\n');

% Blade row
for nrow=1:length(blades);
    nblades = blades{nrow}.nblades;
    if nblades==0
        periodicity = 30;
    else
        periodicity = nblades;
    end
    type = blades{nrow}.type;
    rotation = blades{nrow}.rpm;
    reduction = 0;
    if isfield(blades{nrow},'reduction')
        reduction = blades{nrow}.reduction;
    end
    fprintf(fid,'NI_BEGIN nirow\n');
    fprintf(fid,' NAME			%s\n',['row ' num2str(nrow)]);
    fprintf(fid,' TYPE			%s\n',type);
    fprintf(fid,'PERIODICITY			%i\n',periodicity);
    fprintf(fid,'ROTATION_SPEED			%f\n',rotation);
    fprintf(fid,'NI_BEGIN NIBlade\n');
    fprintf(fid,' NAME			Main Blade\n');
    if isfield(blades{nrow},'tgap')
        fprintf(fid,'NI_BEGIN NITipGap\n');
        if ischar(blades{nrow}.tgap(1))
            fprintf(fid,' DEFINE_SHAPE  %s\n',blades{nrow}.tgap);
        else
            fprintf(fid,' WIDTH_AT_LEADING_EDGE  %f\n',blades{nrow}.tgap(1));
            fprintf(fid,' WIDTH_AT_TRAILING_EDGE  %f\n',blades{nrow}.tgap(2));
        end
        fprintf(fid,'NI_END NITipGap\n');
    end
    if isfield(blades{nrow},'hgap')
        fprintf(fid,'NI_BEGIN NIHubGap\n');
        if ischar(blades{nrow}.hgap(1))
            fprintf(fid,' DEFINE_SHAPE  %s\n',blades{nrow}.hgap);
        else
            fprintf(fid,' WIDTH_AT_LEADING_EDGE  %f\n',blades{nrow}.hgap(1));
            fprintf(fid,' WIDTH_AT_TRAILING_EDGE  %f\n',blades{nrow}.hgap(2));
        end
        fprintf(fid,'NI_END NIHubGap\n');
    end
    fprintf(fid,'NI_BEGIN nibladegeometry\n');
    fprintf(fid,'TYPE	GEOMTURBO\n');
    fprintf(fid,'GEOMETRY_MODIFIED	0\n');
    fprintf(fid,'GEOMETRY TURBO VERSION 5\n');
    fprintf(fid,['blade_expansion_factor_hub		0.06\n',...
        'blade_expansion_factor_shroud		0.03\n',...
        'intersection_npts    10\n',...
        'intersection_control 1\n',...
        'data_reduction                        %i\n',...
        'data_reduction_spacing_tolerance      1e-06\n',...
        'data_reduction_angle_tolerance        90\n'],reduction);
    fprintf(fid,'units		                        1\nnumber_of_blades		        %i\n',nblades);
    if nblades>0
        fprintf(fid,'pressure\nSECTIONAL\n');
        z = blades{nrow}.z;
        r = blades{nrow}.r;
        t = blades{nrow}.t;
        nsec = length(z);
        fprintf(fid,'%i\n',nsec);
        for i=1:nsec
            npoints = length(z{i});
            if even==1
                npoints = npoints-1;
            end
            te_index = blades{nrow}.te_index(i);
            if sqte==1
                te_index=te_index+1;
            end
            fprintf(fid,'#   section %i\nZRTH\n%i\n',i,npoints-te_index+1);
            for j=1:npoints-te_index+1
                fprintf(fid,'%0#6.12f %0#6.12f %0#6.12f\n',(scale.*z{i}(end-j+1)),...
                    (scale.*r{i}(end-j+1)),(t{i}(end-j+1)));
            end
        end
        fprintf(fid,'suction\nSECTIONAL\n%i\n',nsec);
        for i=1:nsec
            te_index = blades{nrow}.te_index(i);
            fprintf(fid,'#   section %i\nZRTH\n%i\n',i,te_index);
            for j=1:te_index
                fprintf(fid,'%0#6.12f %0#6.12f %0#6.12f\n',scale.*z{i}(j),...
                    scale.*r{i}(j),t{i}(j));
            end
        end
    end
    fprintf(fid,'NI_END   nibladegeometry\n');
    fprintf(fid,'  SOLID_BODY_CONFIGURATION	0\n');
    fprintf(fid,'NI_END   NIBlade\n');
    fprintf(fid,'NI_END   nirow\n');
end
fprintf(fid,'NI_END GEOMTURBO\n\n');
fprintf(fid,'NI_BEGIN GEOMETRY\n');
fprintf(fid,'NI_END   GEOMETRY');

fclose(fid);

end