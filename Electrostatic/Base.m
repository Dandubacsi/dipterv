function out = model
%
% Base.m
%
% Model exported on Apr 28 2022, 21:54 by COMSOL 5.5.0.292.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath(['G:\OneDrive - Sch' native2unicode(hex2dec({'00' 'f6'}), 'unicode') 'nherz Zolt' native2unicode(hex2dec({'00' 'e1'}), 'unicode') 'n Koll' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gium\BME\_MSc\5\Dipterv 2\Electrostatic']);

model.label('base.mph');

model.param.set('resolution', '5[um]', 'Capacitance export''s resolution on sample''s surface');
model.param.set('eta', '0.95', 'Capacitance ratio');
model.param.set('mesh_size', '0.25', 'Mesh base size');
model.param.set('w', '25[um]', 'Electrode width');
model.param.set('d', '10[um]', 'Electrode distance from sample');
model.param.set('d_shield', '10[um]', 'Shielding electrode distance from sample');
model.param.set('R', '199.7498[um]', 'Maximum sample size based on point charge field');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 3);

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').lengthUnit([native2unicode(hex2dec({'00' 'b5'}), 'unicode') 'm']);
model.component('comp1').geom('geom1').create('wp1', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp1').label('Sample');
model.component('comp1').geom('geom1').feature('wp1').geom.create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq1').set('size', '1.5*w');
model.component('comp1').geom('geom1').feature('wp1').geom.create('sq2', 'Square');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq2').set('size', '2.5*w');
model.component('comp1').geom('geom1').feature('wp1').geom.create('sq3', 'Square');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('sq3').set('size', '2.5*w+R');
model.component('comp1').geom('geom1').feature('wp1').geom.create('fil1', 'Fillet');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').set('radius', '2.5*w+R');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('fil1').selection('point').set('sq3(1)', 3);
model.component('comp1').geom('geom1').feature('wp1').geom.create('ca1', 'CircularArc');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ca1').set('r', '5*w*exp(0)');
model.component('comp1').geom('geom1').feature('wp1').geom.create('ca2', 'CircularArc');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('ca2').set('r', '5*w*exp(0.25)');
model.component('comp1').geom('geom1').feature('wp1').geom.create('pard1', 'PartitionDomains');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').set('partitionwith', 'edges');
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').selection('domain').set('fil1(1)', 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').selection('edge').set('sq1(1)', [2 3]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').selection('edge').set('sq2(1)', [2 3]);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').selection('edge').set('ca1(1)', 1);
model.component('comp1').geom('geom1').feature('wp1').geom.feature('pard1').selection('edge').set('ca2(1)', 1);
model.component('comp1').geom('geom1').create('wp2', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp2').label('Electrode');
model.component('comp1').geom('geom1').feature('wp2').set('quickz', 'd');
model.component('comp1').geom('geom1').feature('wp2').geom.create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('wp2').geom.feature('sq1').set('size', 'w/2');
model.component('comp1').geom('geom1').create('wp3', 'WorkPlane');
model.component('comp1').geom('geom1').feature('wp3').label('Shield');
model.component('comp1').geom('geom1').feature('wp3').set('quickz', 'd_shield');
model.component('comp1').geom('geom1').feature('wp3').geom.create('sq1', 'Square');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('sq1').set('size', '2*w');
model.component('comp1').geom('geom1').feature('wp3').geom.create('sq2', 'Square');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('sq2').set('size', 'w');
model.component('comp1').geom('geom1').feature('wp3').geom.create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('wp3').geom.feature('dif1').selection('input').set({'sq1'});
model.component('comp1').geom('geom1').feature('wp3').geom.feature('dif1').selection('input2').set({'sq2'});
model.component('comp1').geom('geom1').run('fin');
model.component('comp1').geom('geom1').create('boxsel1', 'BoxSelection');
model.component('comp1').geom('geom1').feature('boxsel1').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('boxsel1').label('Sample_surface');
model.component('comp1').geom('geom1').feature('boxsel1').set('xmin', 0);
model.component('comp1').geom('geom1').feature('boxsel1').set('ymin', 0);
model.component('comp1').geom('geom1').feature('boxsel1').set('ymax', 'inf');
model.component('comp1').geom('geom1').feature('boxsel1').set('zmin', 0);
model.component('comp1').geom('geom1').feature('boxsel1').set('zmax', 0);
model.component('comp1').geom('geom1').create('boxsel2', 'BoxSelection');
model.component('comp1').geom('geom1').feature('boxsel2').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('boxsel2').label('Electrode_surface');
model.component('comp1').geom('geom1').feature('boxsel2').set('xmin', 0);
model.component('comp1').geom('geom1').feature('boxsel2').set('xmax', 'w/2');
model.component('comp1').geom('geom1').feature('boxsel2').set('ymin', 0);
model.component('comp1').geom('geom1').feature('boxsel2').set('ymax', 'w/2');
model.component('comp1').geom('geom1').feature('boxsel2').set('zmin', 0.1);
model.component('comp1').geom('geom1').feature('boxsel2').set('zmax', 'inf');
model.component('comp1').geom('geom1').create('boxsel3', 'BoxSelection');
model.component('comp1').geom('geom1').feature('boxsel3').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('boxsel3').label('Shield_surface');
model.component('comp1').geom('geom1').feature('boxsel3').set('xmin', 'w');
model.component('comp1').geom('geom1').feature('boxsel3').set('xmax', '2*w');
model.component('comp1').geom('geom1').feature('boxsel3').set('ymin', 0);
model.component('comp1').geom('geom1').feature('boxsel3').set('ymax', 'w/2');
model.component('comp1').geom('geom1').feature('boxsel3').set('zmin', 0.1);
model.component('comp1').geom('geom1').feature('boxsel3').set('zmax', 'inf');
model.component('comp1').geom('geom1').create('boxsel4', 'BoxSelection');
model.component('comp1').geom('geom1').feature('boxsel4').set('entitydim', 2);
model.component('comp1').geom('geom1').feature('boxsel4').label('Charge high resolution');
model.component('comp1').geom('geom1').feature('boxsel4').set('xmin', 0);
model.component('comp1').geom('geom1').feature('boxsel4').set('xmax', '1.5*w');
model.component('comp1').geom('geom1').feature('boxsel4').set('ymin', 0);
model.component('comp1').geom('geom1').feature('boxsel4').set('ymax', '1.5*w');
model.component('comp1').geom('geom1').feature('boxsel4').set('zmin', 0);
model.component('comp1').geom('geom1').feature('boxsel4').set('zmax', 0.1);
model.component('comp1').geom('geom1').feature('boxsel4').set('condition', 'inside');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').material.create('mat1', 'Common');
model.component('comp1').material('mat1').selection.geom('geom1', 3);
model.component('comp1').material('mat1').selection.allVoids;

model.component('comp1').physics.create('esbe', 'ElectrostaticsBoundaryElements', 'geom1');
model.component('comp1').physics('esbe').create('gnd1', 'Ground', 2);
model.component('comp1').physics('esbe').feature('gnd1').selection.named('geom1_boxsel1');
model.component('comp1').physics('esbe').create('pot1', 'ElectricPotential', 2);
model.component('comp1').physics('esbe').feature('pot1').selection.named('geom1_boxsel2');
model.component('comp1').physics('esbe').create('pot2', 'ElectricPotential', 2);
model.component('comp1').physics('esbe').feature('pot2').selection.named('geom1_boxsel3');

model.component('comp1').mesh('mesh1').create('size2', 'Size');
model.component('comp1').mesh('mesh1').create('size3', 'Size');
model.component('comp1').mesh('mesh1').create('size4', 'Size');
model.component('comp1').mesh('mesh1').create('fq1', 'FreeQuad');
model.component('comp1').mesh('mesh1').feature('size2').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size2').selection.set([3 4 5]);
model.component('comp1').mesh('mesh1').feature('size3').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size3').selection.set([6]);
model.component('comp1').mesh('mesh1').feature('size4').selection.geom('geom1', 2);
model.component('comp1').mesh('mesh1').feature('size4').selection.set([7]);
model.component('comp1').mesh('mesh1').feature('fq1').selection.remaining;

model.component('comp1').view('view2').axis.set('xmin', -81.19003295898438);
model.component('comp1').view('view2').axis.set('xmax', 343.4398193359375);
model.component('comp1').view('view2').axis.set('ymin', -13.11248779296875);
model.component('comp1').view('view2').axis.set('ymax', 275.3622741699219);

model.component('comp1').material('mat1').label('Air');
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.component('comp1').material('mat1').propertyGroup('def').set('relpermittivity_symmetry', '0');

model.component('comp1').physics('esbe').prop('Symmetry').set('sym1', 'even');
model.component('comp1').physics('esbe').prop('Symmetry').set('sym2', 'even');
model.component('comp1').physics('esbe').feature('gnd1').label('Sample');
model.component('comp1').physics('esbe').feature('pot1').set('V0', '1[V]');
model.component('comp1').physics('esbe').feature('pot1').label('Electrode');
model.component('comp1').physics('esbe').feature('pot2').set('V0', '1[V]');
model.component('comp1').physics('esbe').feature('pot2').label('Shield');

model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 'w/10');
model.component('comp1').mesh('mesh1').feature('size').set('hmin', 'w/20');
model.component('comp1').mesh('mesh1').feature('size2').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size2').set('hmax', 'w/6');
model.component('comp1').mesh('mesh1').feature('size2').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size2').set('hmin', 'w/10');
model.component('comp1').mesh('mesh1').feature('size2').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('size3').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size3').set('hmax', 'w*exp(0*mesh_size)*(1-exp(-mesh_size))');
model.component('comp1').mesh('mesh1').feature('size3').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size3').set('hmin', '0.5*w*exp(0*mesh_size)*(1-exp(-mesh_size))');
model.component('comp1').mesh('mesh1').feature('size3').set('hminactive', true);
model.component('comp1').mesh('mesh1').feature('size4').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size4').set('hmax', 'w*exp(1*mesh_size)*(1-exp(-mesh_size))');
model.component('comp1').mesh('mesh1').feature('size4').set('hmaxactive', true);
model.component('comp1').mesh('mesh1').feature('size4').set('hmin', '0.5*w*exp(1*mesh_size)*(1-exp(-mesh_size))');
model.component('comp1').mesh('mesh1').feature('size4').set('hminactive', true);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('i1', 'Iterative');
model.sol('sol1').feature('s1').create('d1', 'Direct');
model.sol('sol1').feature('s1').feature('i1').create('dp1', 'DirectPreconditioner');

model.study('std1').label('Static capacitances');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('s1').set('stol', '1e-12');
model.sol('sol1').feature('s1').feature('fcDef').set('linsolver', 'i1');
model.sol('sol1').feature('s1').feature('aDef').set('matrixformat', 'filled');
model.sol('sol1').feature('s1').feature('i1').label('Iterative Solver');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
model.sol('sol1').feature('s1').feature('i1').set('prefuntype', 'right');
model.sol('sol1').feature('s1').feature('i1').feature('dp1').set('linsolver', 'dense');
model.sol('sol1').feature('s1').feature('d1').label('Direct Dense Solver');
model.sol('sol1').feature('s1').feature('d1').set('linsolver', 'dense');
model.sol('sol1').runAll;

out = model;
