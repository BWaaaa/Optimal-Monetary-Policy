doc = struct;
doc.title = 'Dynare Prctice on New Keysian Simple Models: Rottemberg(1982) and Calvo(1983)';
doc.author = 'Brian Wang';
doc.date = date;
doc.mod = 'corresponding model';
doc.notes = 'correspodong explainative files';

addpath D:\...\4.5.7\matlab 
% Where you instal Dynare\4.5.7(or any latest version)\matlab
dynare NK_CalvoFairy_DeviatTaylor.mod
% Here run the corresponding_model_file.mod
% By default, dynare solve questions and store results in its \work folders. 
% Therefore, please save both this_implement_file.m and corresponding_model_file.mod in dynare\work, and keep current path at dynare\work
