clear all;
output_name = '20_stage'; %Output file name
input_name = '20_stage.csv'; %Geometry file name

%Set a few constant parameters
t_TE = 0.0003; %Trailing edge thickness
g_tip = 0.0003; %Tip clearance
RPM = 6782; %Rotational speed in rpm
epsilon = 0.0001; %some small number to ensure overlap
%Unpacking the geometry data
geom = csvread(input_name); %Read from the csv
geom_size = size(geom); %Get the dimensions of the resulting array
n_stages = geom_size(2)/6; %Number of stages
n_blades = geom(1, 1:2*n_stages); %Number of blades in each row
X1 = geom(2, :); %Inlet angles for each blade. 3 angles per blade for hub, midspan and tip
X2 = geom(3, :); %Exit angles for each blade. 3 angles per blade for hub, midspan and tip
Cx = geom(4, 1:2*n_stages); %Blade axial chords
hub_x = geom(5, 1:2*n_stages+1); %Hub axial coordinates
hub_r = geom(6, 1:2*n_stages+1); %Hub radial coordinates
case_x = geom(7, 1:2*n_stages+1); %Case axial coordinates
case_r = geom(8, 1:2*n_stages+1); %Case radial coordinates

casing = [transpose(case_x) transpose(case_r)]; %Make the coordinates an nx2 double
hub = [transpose(hub_x) transpose(hub_r)];
z = 0; %Start the axial position
for i=1:n_stages*2
    
    blades{i}.nblades = n_blades(i);
    blades{i}.type = 'normal';
    blades{i}.reduction = 1; 
    rm = (hub_r(i)+case_r(i))/2; %Stage mean radius
    
    for j=1:3 %1-3 for hub, midspan, tip
        [Z, RTH, te_index] = Profile(X1((i-1)*3+j), X2((i-1)*3+j), t_TE, Cx(i)); %Get the profile
        blades{i}.te_index(j) = te_index; %Index of point in middle of TE
        blades{i}.z{j} = Z+0.25*Cx(i)+z; %Add a quater chord to the blade position
        R = zeros(1, length(Z)); %Initialise the radial coordinate
        TH = zeros(1, length(Z)); %Initialise the tangential coordinate
        if j == 1 %Hub
            for k=1:length(Z)
                %Interpolating to be just below the hub line
                R(k) = (1-epsilon)*(hub_r(i)+(hub_r(i+1)-hub_r(i))*(Z(k)-hub_x(i))/(hub_x(i+1)-hub_x(i)));
                TH(k) = RTH(k)/R(k); %Convert the distance to an angle
            end
        end
        if j == 2 %Midspan
            for k=1:length(Z)
                R(k) = rm; %Constant radius
                TH(k) = RTH(k)/R(k);
            end
        end
        if j == 3 %Casing
            for k=1:length(Z)
                %Interpolating to be just above the casing line
                R(k) = (1+epsilon)*(case_r(i)+(case_r(i+1)-case_r(i))*(Z(k)-case_x(i))/(case_x(i+1)-case_x(i)));
                TH(k) = RTH(k)/R(k); %Convert the distance to an angle
            end
        end
        blades{i}.r{j} = R; %Set the radial coordinates
        blades{i}.t{j} = TH; %Set the tangential coordinates
    end
    
    if mod(i,2)==1 %Stator
        blades{i}.rpm = 0; %No rotation
        r1_blade = (hub_r(i)+(hub_r(i+1)-hub_r(i))*(z+0.125*Cx(i)+g_tip/2-hub_x(i))/(hub_x(i+1)-hub_x(i)))...
            *(1+epsilon); %Radius of first point on shroud attached to blade
        r2_blade = r1_blade - 0.001; %1mm shroud thickness
        r5_blade = (hub_r(i)+(hub_r(i+1)-hub_r(i))*(z+1.375*Cx(i)-g_tip/2-hub_x(i))/(hub_x(i+1)-hub_x(i)))...
            *(1+epsilon); %Radius of last point on shroud attached to blade
        r4_blade = r5_blade - 0.001; %1mm shroud thickness
        r3_blade = (r2_blade+r4_blade)/2;
        r1_hub = (hub_r(i)+(hub_r(i+1)-hub_r(i))*(z+0.125*Cx(i)-g_tip/2-hub_x(i))/(hub_x(i+1)-hub_x(i)))...
            *(1+epsilon); %Radius of first point on shroud attached to blade
        r2_hub = r1_hub - 0.001 - g_tip; %1mm shroud thickness
        r5_hub = (hub_r(i)+(hub_r(i+1)-hub_r(i))*(z+1.375*Cx(i)+g_tip/2-hub_x(i))/(hub_x(i+1)-hub_x(i)))...
            *(1+epsilon); %Radius of last point on shroud attached to blade
        r4_hub = r5_hub - 0.001 - g_tip; %1mm shroud thickness
        r3_hub = (r2_hub+r4_hub)/2;
        z1_blade = z+0.125*Cx(i)+g_tip/2;
        z3_blade = z+1.375*Cx(i)-g_tip/2;
        z2_blade = (z1_blade+z3_blade)/2;
        z1_hub = z+0.125*Cx(i)-g_tip/2;
        z3_hub = z+1.375*Cx(i)+g_tip/2;
        z2_hub = z2_blade;
        z_sep = z2_blade;
        r1_sep = r3_blade;
        r2_sep = r3_hub;
        shrouds.shroud{i}.blade = [z1_blade r1_blade; z1_blade r2_blade; z2_blade r3_blade;...
            z3_blade r4_blade; z3_blade r5_blade];
        shrouds.shroud{i}.cas = [z1_hub r1_hub; z1_hub r2_hub; z2_hub r3_hub; z3_hub r4_hub; z3_hub r5_hub];
        shrouds.shroud{i}.sep = [z_sep r1_sep; z_sep r2_sep];
    else
        blades{i}.rpm = RPM; %Set rpm
        r1_blade = (case_r(i)+(case_r(i+1)-case_r(i))*(z+0.125*Cx(i)+g_tip/2-case_x(i))/(case_x(i+1)-case_x(i)))...
            *(1-epsilon); %Radius of first point on shroud attached to blade
        r2_blade = r1_blade + 0.001; %1mm shroud thickness
        r5_blade = (case_r(i)+(case_r(i+1)-case_r(i))*(z+1.375*Cx(i)-g_tip/2-case_x(i))/(case_x(i+1)-case_x(i)))...
            *(1-epsilon); %Radius of last point on shroud attached to blade
        r4_blade = r5_blade + 0.001; %1mm shroud thickness
        r3_blade = (r2_blade+r4_blade)/2;
        r1_cas = (case_r(i)+(case_r(i+1)-case_r(i))*(z+0.125*Cx(i)-g_tip/2-case_x(i))/(case_x(i+1)-case_x(i)))...
            *(1-epsilon); %Radius of first point on shroud attached to blade
        r2_cas = r1_cas + 0.001 + g_tip; %1mm shroud thickness
        r5_cas = (case_r(i)+(case_r(i+1)-case_r(i))*(z+1.375*Cx(i)+g_tip/2-case_x(i))/(case_x(i+1)-case_x(i)))...
            *(1-epsilon); %Radius of last point on shroud attached to blade
        r4_cas = r5_cas + 0.001 + g_tip; %1mm shroud thickness
        r3_cas = (r2_cas+r4_cas)/2;
        z1_blade = z+0.125*Cx(i)+g_tip/2;
        z3_blade = z+1.375*Cx(i)-g_tip/2;
        z2_blade = (z1_blade+z3_blade)/2;
        z1_cas = z+0.125*Cx(i)-g_tip/2;
        z3_cas = z+1.375*Cx(i)+g_tip/2;
        z2_cas = z2_blade;
        z_sep = z2_blade;
        r1_sep = r3_blade;
        r2_sep = r3_cas;
        shrouds.shroud{i}.blade = [z1_blade r1_blade; z1_blade r2_blade; z2_blade r3_blade;...
            z3_blade r4_blade; z3_blade r5_blade];
        shrouds.shroud{i}.cas = [z1_cas r1_cas; z1_cas r2_cas; z2_cas r3_cas; z3_cas r4_cas; z3_cas r5_cas];
        shrouds.shroud{i}.sep = [z_sep r1_sep; z_sep r2_sep];
    end
    
    z = z+1.5*Cx(i); %Increment for the next blade row
end

scale = 1;
varargin.scale = scale;
varargin.shrouds = shrouds;
write_geomTurbo(output_name,hub,casing,blades,varargin)