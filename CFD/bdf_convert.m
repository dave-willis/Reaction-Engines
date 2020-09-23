output_name = 'output_lean_m10_stagger_p20';

bdf_name = 'Blade_Definition';

bdf = bdfscan(bdf_name);

blades{1}.z = bdf.z;
blades{1}.r = bdf.r;

toffset = -0.1*3.14/11;
staggeroffset = 0.2*3.14/11;
figure()
for i =1:21
    
    plot(bdf.z{i},bdf.t{i},'-xb')
    hold on
    fraction = -0.5+(bdf.z{i}-min(bdf.z{i}))/(max(bdf.z{i})-min(bdf.z{i}))
    bdf.t{i}=bdf.t{i}+toffset*(i-1)/21 +staggeroffset*fraction*((21.-i)/20.-0.5)
    plot(bdf.z{i},bdf.t{i},'-r')
    %show()
    
end

blades{1}.t = bdf.t;

blades{1}.nblades = bdf.nblades;
blades{1}.te_index = bdf.te_index;
blades{1}.type = 'normal';
blades{1}.rpm = 0;
blades{1}.reduction = 1;

casing = bdf.casing;
hub = bdf.hub;

scale = 1;

write_geomTurbo(output_name,hub,casing,blades,scale)