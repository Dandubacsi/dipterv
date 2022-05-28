function out = Geometry_updater(model)

% clear former model
model.component.clear
model.study.clear
model.sol.clear
model.result.clear
model.result.dataset.clear

% extracting parameters from modell
%--------------------------------------------------------------------------
w = char(model.param.get('w'));
d = char(model.param.get('d'));
eta = double(string(model.param.get('eta')));
mesh_size = double(string(model.param.get('mesh_size')));
resolution = char(model.param.get('resolution'));

w = double(string(w(1:end-4)));
d = double(string(d(1:end-4)));
resolution = double(string(resolution(1:end-4)));

% calculating and setting the necessary simulational space
%--------------------------------------------------------------------------

R = d*sqrt(1/(1-eta)^2-1);
model.param.set('R',string(R)+'[um]','Maximum sample size based on point charge field');

% creating new component
%--------------------------------------------------------------------------

comp1 = model.component.create('comp1',true);

% creating new geometry
%--------------------------------------------------------------------------
geom = model.component('comp1').geom.create('geom1', 3);
geom.lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);

% number of stray strips
num_strips = floor(log((1+R/w)/5)/mesh_size)-1;

% create sample surface
wp = geom.create('wp1', 'WorkPlane');
wp.set('quickz', '0');
wp.label('Sample');

wp.geom.create('sq1','Square');
wp.geom.feature('sq1').set('size','1.5*w');
wp.geom.create('sq2','Square');
wp.geom.feature('sq2').set('size','2.5*w');
wp.geom.create('sq3', 'Square');
wp.geom.feature('sq3').set('size', '2.5*w+R');
wp.geom.create('fil1','Fillet');
wp.geom.feature('fil1').selection('point').set('sq3',3);
wp.geom.feature('fil1').set('radius','2.5*w+R');

for i = 0:num_strips
    s = 'ca'+string(i+1);
    wp.geom.create(s, 'CircularArc');
	wp.geom.feature(s).set('r', '5*w*exp('+string(i*mesh_size)+')');
end

wp.geom.create('pard1', 'PartitionDomains');
wp.geom.feature('pard1').set('partitionwith', 'edges');
wp.geom.feature('pard1').selection('domain').set('fil1', 1);
wp.geom.feature('pard1').selection('edge').set('sq1', [2,3]);
wp.geom.feature('pard1').selection('edge').set('sq2', [2,3]);

for i = 0:num_strips
    s = 'ca'+string(i+1);
    wp.geom.feature('pard1').selection('edge').set(s, 1);
end

% create electrode surface
wp = geom.create('wp2', 'WorkPlane');
wp.set('quickz', 'd');
wp.label('Electrode');

wp.geom.create('sq1','Square');
wp.geom.feature('sq1').set('size','w/2');

% create shield surface
wp = geom.create('wp3', 'WorkPlane');
wp.set('quickz', 'd_shield');
wp.label('Shield');

wp.geom.create('sq1','Square');
wp.geom.feature('sq1').set('size','2*w');
wp.geom.create('sq2','Square');
wp.geom.feature('sq2').set('size','w');
wp.geom.create('dif1', 'Difference');
wp.geom.feature('dif1').selection('input').set('sq1');
wp.geom.feature('dif1').selection('input2').set('sq2');

% selections for surfaces
geom.create('boxsel1', 'BoxSelection');
geom.feature('boxsel1').set('entitydim', 2);
geom.feature('boxsel1').label('Sample_surface');
geom.feature('boxsel1').set('xmin', 0);
geom.feature('boxsel1').set('ymin', 0);
geom.feature('boxsel1').set('ymax', 'inf');
geom.feature('boxsel1').set('zmin', 0);
geom.feature('boxsel1').set('zmax', 0);

geom.create('boxsel2', 'BoxSelection');
geom.feature('boxsel2').set('entitydim', 2);
geom.feature('boxsel2').label('Electrode_surface');
geom.feature('boxsel2').set('xmin', 0);
geom.feature('boxsel2').set('xmax', 'w/2');
geom.feature('boxsel2').set('ymin', 0);
geom.feature('boxsel2').set('ymax', 'w/2');
geom.feature('boxsel2').set('zmin', 0.1);
geom.feature('boxsel2').set('zmax', 'inf');

geom.create('boxsel3', 'BoxSelection');
geom.feature('boxsel3').set('entitydim', 2);
geom.feature('boxsel3').label('Shield_surface');
geom.feature('boxsel3').set('xmin', 'w');
geom.feature('boxsel3').set('xmax', '2*w');
geom.feature('boxsel3').set('ymin', 0);
geom.feature('boxsel3').set('ymax', 'w/2');
geom.feature('boxsel3').set('zmin', 0.1);
geom.feature('boxsel3').set('zmax', 'inf');

geom.create('boxsel4', 'BoxSelection');
geom.feature('boxsel4').set('entitydim', 2);
geom.feature('boxsel4').label('Charge high resolution');
geom.feature('boxsel4').set('xmin', 0);
geom.feature('boxsel4').set('xmax', '1.5*w');
geom.feature('boxsel4').set('ymin', 0);
geom.feature('boxsel4').set('ymax', '1.5*w');
geom.feature('boxsel4').set('zmin', 0);
geom.feature('boxsel4').set('zmax', 0.1);
geom.feature('boxsel4').set('condition', 'inside');

geom.run;
geom.run('fin');

% setting up the materials
%--------------------------------------------------------------------------
mat = model.component('comp1').material.create('mat1', 'Common');
mat.selection.geom('geom1', 3);
mat.selection.allVoids;
mat.propertyGroup('def').set('relpermittivity', {'1'});
mat.label('Air');

% setting up the physics
%--------------------------------------------------------------------------
physics = comp1.physics.create('esbe', 'ElectrostaticsBoundaryElements', 'geom1');
physics.prop('ShapeProperty').set('shapeorder', 'p21');
physics.create('gnd1', 'Ground', 2);
physics.feature('gnd1').label('Sample');
physics.feature('gnd1').selection.named('geom1_boxsel1');
physics.create('pot1', 'ElectricPotential', 2);
physics.feature('pot1').label('Electrode');
physics.feature('pot1').selection.named('geom1_boxsel2');
physics.feature('pot1').set('V0', '1[V]');
physics.create('pot2', 'ElectricPotential', 2);
physics.feature('pot2').label('Shield');
physics.feature('pot2').selection.named('geom1_boxsel3');
physics.feature('pot2').set('V0', '1[V]');
physics.prop('Symmetry').set('sym1', 'even');
physics.prop('Symmetry').set('sym2', 'even');

physics.feature('pot1').label('Electrode');
physics.feature('pot2').label('Shield');

% setting up the mesh
%--------------------------------------------------------------------------

mesh = comp1.mesh.create('mesh1');
mesh.feature('size').set('custom', 'on');
mesh.feature('size').set('hmax', 'w/20');
mesh.feature('size').set('hmin', 'w/30');
mesh.create('size2','Size');
mesh.feature('size2').selection.geom('geom1',2);
mesh.feature('size2').selection.set(1);
mesh.feature('size2').set('custom', 'on');
mesh.feature('size2').set('hmax', 'resolution/2');
mesh.feature('size2').set('hmin', 'resolution/4');
mesh.create('size3','Size');
mesh.feature('size3').selection.geom('geom1',2);
mesh.feature('size3').selection.set([4 5]);
mesh.feature('size3').set('custom', 'on');
mesh.feature('size3').set('hmax', 'w/8');
mesh.feature('size3').set('hmin', 'w/10');

for i = 0:num_strips
    s = 'size'+string(i+4);
    mesh.create(s, 'Size');
    mesh.feature(s).selection.geom('geom1', 2);
    selection = i+6;
	mesh.feature(s).selection.set(selection);
    mesh.feature(s).set('custom', 'on');
    max_element_size = 'w*exp('+string(i)+'*mesh_size)*(1-exp(-mesh_size))';
    min_element_size = '0.5*w*exp('+string(i)+'*mesh_size)*(1-exp(-mesh_size))';
    mesh.feature(s).set('hmax', max_element_size);
    mesh.feature(s).set('hmaxactive', true);
    mesh.feature(s).set('hmin', min_element_size);
    mesh.feature(s).set('hminactive', true);
end

mesh.create('fq1', 'FreeQuad');
mesh.feature('fq1').selection.remaining;  

% setting up the study
%--------------------------------------------------------------------------

std1 = model.study.create('std1');
std1.create('stat', 'Stationary');
std1.feature('stat').activate('esbe', true);
std1.label('Static capacitances');
std1.setGenPlots(false);

sol1 =  model.sol.create('sol1');
sol1.study('std1');
sol1.attach('std1');
sol1.create('st1', 'StudyStep');
sol1.feature('st1').set('study', 'std1');
sol1.feature('st1').set('studystep', 'stat');
sol1.create('v1', 'Variables');
sol1.feature('v1').set('control', 'stat');

sol1.create('s1', 'Stationary');
sol1.feature('s1').create('i1', 'Iterative');
sol1.feature('s1').create('d1', 'Direct');
sol1.feature('s1').feature('i1').create('dp1', 'DirectPreconditioner');
sol1.feature('s1').set('stol', '1e-12');
sol1.feature('s1').feature('fcDef').set('linsolver', 'i1');
sol1.feature('s1').feature('i1').label('Iterative Solver');
sol1.feature('s1').feature('i1').set('linsolver', 'cg');
sol1.feature('s1').feature('i1').set('prefuntype', 'right');
sol1.feature('s1').feature('i1').feature('dp1').set('linsolver', 'dense');
sol1.feature('s1').feature('d1').label('Direct Dense Solver');
sol1.feature('s1').feature('d1').set('linsolver', 'dense');
sol1.feature('s1').feature('aDef').set('matrixformat', 'filled');
sol1.feature('s1').feature('i1').active(true);

% setting the results
%--------------------------------------------------------------------------

% grid for the potential calculations
grid1 = model.result.dataset.create('grid1', 'Grid3D');
grid1.set('source', 'data');
grid1.set('data', 'dset1');
grid1.set('par1', 'x');
grid1.set('par2', 'y');
grid1.set('par3', 'z');
grid1.set('parmin1', 0);
grid1.set('parmax1', '2.5*w');
grid1.set('parmin2', 0);
grid1.set('parmax2', '2.5*w');
grid1.set('parmin3', 0);
grid1.set('parmax3', '3*d_shield');
grid1.set('res1', 30);
grid1.set('res2', 30);
grid1.set('res3', 30);

% cutlines for the charge distributions
cl1 = model.result.dataset.create('cl1','CutLine3D');
cl1.set('genpoints', {'0' '0' 'd'; 'w/2' '0' 'd'});
cl1.label('Electrode_CutLine');

cl2 = model.result.dataset.create('cl2','CutLine3D');
cl2.set('genpoints', {'0' '0' '0'; 'w+R' '0' '0'});
cl2.label('Sample_CutLine');

cl3 = model.result.dataset.create('cl3','CutLine3D');
cl3.set('genpoints', {'w' '0' 'd_shield'; '2*w' '0' 'd_shield'});
cl3.label('Shield_CutLine');

% Surface for charge distribution on sample
surf1 = model.result.dataset.create('surf1', 'Surface');
surf1.label('Sample surface');
surf1.selection.named('geom1_boxsel1');

% Surface for charge distribution on shield
surf2 = model.result.dataset.create('surf2', 'Surface');
surf2.label('Shield surface');
surf2.selection.named('geom1_boxsel3');

% Surface for sensitivity
surf3 = model.result.dataset.create('surf3', 'Surface');
surf3.label('Sensitivity');
surf3.selection.set(1);

% Surface for electric field
cpl1 = model.result.dataset.create('cpl1', 'CutPlane');
cpl1.set('data', 'grid1');
cpl1.label('Electric field CutPlane');
cpl1.set('planetype', 'general');
cpl1.set('genmethod', 'pointnormal');
cpl1.set('genpnvec', [0 -1 0]);

% Electric potential
pg1 = model.result.create('pg1', 'PlotGroup3D');
pg1.label('Electric Potential');
pg1.set('data', 'grid1');

pg1.create('slc1', 'Slice');
pg1.feature('slc1').set('expr', {'esbe.V'});
pg1.feature('slc1').set('planetype', 'quick');
pg1.feature('slc1').set('quickplane', 'yz');
pg1.feature('slc1').set('quickxmethod', 'coord');
pg1.feature('slc1').set('smooth', 'internal');
pg1.feature('slc1').set('colortable', 'RainbowLight');

pg1.create('iso1', 'Isosurface');
pg1.feature('iso1').set('smooth', 'internal');
pg1.feature('iso1').set('allowmaterialsmoothing', false);
pg1.feature('iso1').set('inheritplot', 'slc1');
pg1.feature('iso1').set('resolution', 'normal');
pg1.feature('iso1').set('number', 7);

pg1.active(false);

% Surface charge density
pg2 = model.result.create('pg2', 'PlotGroup2D');
pg2.label('Sample charge density');
pg2.create('con1', 'Contour');
pg2.feature('con1').set('expr', 'esbe.nD');
pg2.feature('con1').set('unit', 'nC/m^2');
pg2.feature('con1').set('number', 50);
pg2.feature('con1').set('colortablerev', true);
pg2.feature('con1').set('contourtype', 'filled');

% Charge distributions
pg3 = model.result.create('pg3', 'PlotGroup1D');
pg3.label('Charge distribution electrode');
pg3.set('data', 'cl1');
pg3.create('lngr1', 'LineGraph');
pg3.feature('lngr1').set('expr', 'comp1.esbe.nD');
pg3.feature('lngr1').set('unit', 'nC/m^2');

pg4 = model.result.create('pg4', 'PlotGroup1D');
pg4.label('Charge distribution sample');
pg4.set('data', 'cl2');
pg4.create('lngr1', 'LineGraph');
pg4.feature('lngr1').set('expr', 'comp1.esbe.nD');
pg4.feature('lngr1').set('unit', 'nC/m^2');

pg5 = model.result.create('pg5', 'PlotGroup1D');
pg5.label('Charge distribution shield');
pg5.set('data', 'cl3');
pg5.create('lngr1', 'LineGraph');
pg5.feature('lngr1').set('expr', 'comp1.esbe.nD');
pg5.feature('lngr1').set('unit', 'nC/m^2');

% Electric field 
pg6 = model.result.create('pg6', 'PlotGroup2D');
pg6.create('con1', 'Contour');
pg6.create('str1', 'Streamline');
pg6.label('Electric field');
pg6.set('data', 'cpl1');
pg6.feature('con1').set('number', 10);
pg6.feature('con1').set('levelrounding', true);
pg6.feature('con1').set('contourtype', 'filled');
pg6.feature('con1').set('colortable', 'RainbowLight');
pg6.feature('con1').set('smooth', 'internal');
pg6.feature('con1').set('allowmaterialsmoothing', false);
pg6.feature('con1').set('resolution', 'normal');
pg6.feature('str1').set('expr', {'esbe.Ex' 'esbe.Ez'});
pg6.feature('str1').set('posmethod', 'magnitude');
pg6.feature('str1').set('pointtype', 'arrow');
pg6.feature('str1').set('arrowcount', 165);
pg6.feature('str1').set('arrowtype', 'arrow');
pg6.feature('str1').set('arrowlength', 'logarithmic');
pg6.feature('str1').set('arrowscale', 3.75E-5);
pg6.feature('str1').set('smooth', 'internal');
pg6.feature('str1').set('allowmaterialsmoothing', false);
pg6.feature('str1').set('arrowcountactive', false);
pg6.feature('str1').set('arrowscaleactive', false);
pg6.feature('str1').set('resolution', 'normal');

% Sensitivity
pg7 = model.result.create('pg7', 'PlotGroup2D');
pg7.label('Sensitivity');
pg7.create('con1', 'Contour');
pg7.set('data', 'surf3');
pg7.feature('con1').set('expr', '-esbe.nD');
pg7.feature('con1').set('unit', 'nC/m^2');
pg7.feature('con1').set('number', 50);
pg7.feature('con1').set('colortablerev', false);
pg7.feature('con1').set('contourtype', 'filled');

% Sample charge distribution high resolution
eg1 = model.result.evaluationGroup.create('eg1', 'EvaluationGroup');
eg1.set('data', 'dset1');
eg1.label('Charge distribution high resolution');
eg1.set('transpose', true);

% creating a square tile on the sample's surface around the measuring
% electrode using it's simmetry

% the indexing of the tiles is the following:

% y
% |      9
% |    7 8
% |  4 5 6
% |0 1 2 3
%  -------- x

index = 0;
for i = 0:floor(1.5*w/resolution)-1
    for j = i:floor(1.5*w/resolution)-1
        s = 'int'+string(index);
        eg1.create(s, 'IntSurface');
        left = string(j*resolution)+' [um]';
        right = string((j+1)*resolution)+' [um]';
        down = string(i*resolution)+' [um]';
        up = string((i+1)*resolution)+' [um]';
        eg1.feature(s).set('expr', '-4*esbe.nD*(x>='+left+' && x<'+right+' && y>='+down+' && y<'+up+')');
        eg1.feature(s).set('unit', {'fC'});
        eg1.feature(s).set('descr', 'Charge in box '+string(index));
        eg1.feature(s).selection.named('geom1_boxsel4');
        index = index + 1;
    end
end

% Capacitance coefficients
eg2 = model.result.evaluationGroup.create('eg2', 'EvaluationGroup');
eg2.set('data', 'dset1');
eg2.label('Capacitance coefficients 1');
eg2.create('int1', 'IntSurface');
eg2.feature('int1').label('C_electrode');
eg2.feature('int1').set('expr', {'4*esbe.nD/1[V]'});
eg2.feature('int1').set('unit', {'fF'});
eg2.feature('int1').set('descr', {'C10'});
eg2.feature('int1').selection.named('geom1_boxsel2');
eg2.create('int2', 'IntSurface');
eg2.feature('int2').label('C_shield');
eg2.feature('int2').set('expr', {'4*esbe.nD/1[V]'});
eg2.feature('int2').set('unit', {'fF'});
eg2.feature('int2').set('descr', {'C20'});
eg2.feature('int2').selection.named('geom1_boxsel3');

eg3 = model.result.evaluationGroup.create('eg3', 'EvaluationGroup');
eg3.set('data', 'dset1');
eg3.label('Capacitance coefficients 2');
eg3.create('int1', 'IntSurface');
eg3.feature('int1').label('C_mutal');
eg3.feature('int1').set('expr', {'-4*esbe.nD/1[V]'});
eg3.feature('int1').set('unit', {'fF'});
eg3.feature('int1').set('descr', {'C12'});
eg3.feature('int1').selection.named('geom1_boxsel3');

out = model;
end

