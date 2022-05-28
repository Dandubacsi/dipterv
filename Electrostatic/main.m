% main code
close all
clear all
clc

% Importing libraries
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.showProgress(true);

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% configuration space description
d_list = 10:0.5:50; % electrode sample distances
w_list = [25, 50, 75]; % electrode widths 
amp = 0; % excitation amplitude
resolution = 5; % sample discretization resolution for charge calculation
eta = 0.95; % defines the model space size
mesh_size = 0.25; % defines the mesh fineness (smaller better)

% path for model files
path = 'G:\Batch_models\';

% create model file
model = ModelUtil.create('Model');
model.modelPath('G:\OneDrive - Schönherz Zoltán Kollégium\BME\_MSc\5\Dipterv 2\Electrostatic\');
model.label('base.mph');

% set resolution for charge distribution extraction
model.param.set('resolution',string(resolution)+'[um]',"Capacitance export's resolution on sample's surface");
       
% eta defines the model space size
model.param.set('eta', string(eta), 'Capacitance ratio');

% defining the mesh sizeing
model.param.set('mesh_size', string(mesh_size), 'Mesh base size');

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% create log file
fid = fopen('log.txt','w');
fprintf(fid, '%s: %s\r\n', datestr(now, 0), 'Capacitance batch calculations.');
if fid == -1
    error('Cannot open log file.');
end
fclose(fid);

% create table for results
number_of_rows = size(w_list,2)*size(d_list,2)*(2*amp+1);
number_of_columns = 7+floor(1.5*max(w_list)/resolution)*floor(1.5*max(w_list)/resolution+1)/2;

vartypes = cell([1,number_of_columns]);
varnames = cell([1,number_of_columns]);

vartypes{1} = 'string';

for i = 2:number_of_columns
    vartypes{i} = 'doublenan';
end

varnames{1} = 'Solution time (hh:mm:ss)';
varnames{2} = 'Electrode width (um)';
varnames{3} = 'Electrode separation (um)';
varnames{4} = 'Shield separation (um)';
varnames{5} = 'Electrode capacitance (fF)';
varnames{6} = 'Shield capacitance (fF)';
varnames{7} = 'Mutual capacitance (fF)';

for i = 8:number_of_columns
    varnames{i} = 'Charges within box '+string(i-8)+' (fC)';
end

results = table('Size',[number_of_rows,number_of_columns],'VariableTypes',vartypes,'VariableNames',string(varnames));

% processing

run_count = 1;
for i = 1:length(w_list)
    for j = 1:length(d_list)        
        model.param.set('w', string(w_list(i))+'[um]','Electrode width');
        model.param.set('d', string(d_list(j))+'[um]','Electrode distance from sample');
        model.param.set('d_shield', string(d_list(j))+'[um]', 'Shielding electrode distance from sample');
        
        % update geometry
        model = Geometry_updater(model);
        
        for k = 1:2*amp+1
            try  
                tic
                
                % update geometry parameters
                dist = d_list(j)-amp + k - 1;
                shield = dist;
                model.param.set('d', string(dist)+'[um]');
                model.param.set('d_shield', string(shield)+'[um]');
                                
                % log state
                disp('Calculating w = '+string(w_list(i))+' um, d = '+string(dist)+' um.');

                fid = fopen('log.txt','a');
                Msg = 'Calculating w = '+string(w_list(i))+' um, d = '+string(dist)+' um.';
                fprintf(fid, '%s: %s\r\n', datestr(now, 0), Msg);
                fclose(fid);

                % update and solve model
                model.physics('esbe').feature('pot2').set('V0', '1[V]');
                model.sol('sol1').runAll;

                % calculate results
                model.result.evaluationGroup('eg2').run;
                
                % extract results
                cap = mphtable(model,'eg2');
                c10 = cap.data(1);
                c20 = cap.data(2);
                
                % update and solve model
                model.physics('esbe').feature('pot2').set('V0', '0[V]');
                model.sol('sol1').runAll;
                
                % calculate results
                model.result.evaluationGroup('eg1').run;
                model.result.evaluationGroup('eg3').run;
                
                % extract result
                charges = mphtable(model,'eg1');
                charges = charges.data';
                cap = mphtable(model,'eg3');
                c12 = cap.data(1);
                
                % save model
                filename = 'Capacitances_w='+string(w_list(i))+'um_d_='+string(dist)+'um_d_shield='+string(shield)+'um.mph';
                mphsave(model,path+filename);

                % stop timer and format
                time = toc;
                time = seconds(time);
                time.Format = 'hh:mm:ss'; 

                disp('w = '+string(w_list(i))+' um, d = '+string(dist)+' um, d_shield = '+string(shield)+' um is complete.');
                disp('solution took: '+string(time));

                % log progress
                fid = fopen('log.txt','a');
                Msg = 'w = '+string(w_list(i))+' um, d = '+string(dist)+' um, d_shield = '+string(shield)+' um is complete.';
                fprintf(fid, '%s: %s\r\n', datestr(now, 0), Msg);
                Msg = 'Solution took: '+string(time);
                fprintf(fid, '%s: %s\r\n', datestr(now, 0), Msg);
                fclose(fid);
                
                % save results to variable
                results{run_count,1} = string(time);
                results{run_count,2} = w_list(i);
                results{run_count,3} = dist;
                results{run_count,4} = shield;
                results{run_count,5} = c10;
                results{run_count,6} = c20;
                results{run_count,7} = c12;
                results{run_count,8:7+size(charges,2)} = charges;
                
                % increase run count
                run_count = run_count+1;
                
            catch ME
                yourMsg = getReport(ME);
                fid = fopen('log.txt','a');
                fprintf(fid, '\r\n%s: %s\r\n', datestr(now, 0), yourMsg);
                fclose(fid);

                % save model
                filename = 'Capacitances_w='+string(w_list(i))+'um_d_='+string(dist)+'um_d_shield='+string(shield)+'um.mph';
                mphsave(model,path+filename);

                % display error
                disp('Error during w = '+string(w_list(i))+' um, d = '+string(dist)+' um, d_shield = '+string(shield)+' um!');
            end
        end
    end
end

% save results to drive
writetable(results,'sim_res_w25_75.csv');  

% shut down pc after simulation is done
eval(['!shutdown -s -f -t ' num2str(60)])
