function [] = calculate_LARS_PCE_parcel(parcel_idx, name, XED,YED,p,q, output_folder)
% parcel_idx, name : only for the nale of the files
% input_param_indices : input parameters kept after Morris (out of 145)
% p,q : truncation scheme
cd(output_folder);

MetaOpts.ExpDesign.X = XED; 
MetaOpts.ExpDesign.Y = YED ;

% Degree adaptive PCE is autom. enabled if MetaOpts.Degree is an array
MetaOpts.Bootstrap.Replications = 1000;
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'LARS';
MetaOpts.Degree = 1:p;
MetaOpts.TruncOptions.qNorm = q;

% Calculate and save the metamodel
myPCE_LARS = uq_createModel(MetaOpts);
uq_print(myPCE_LARS)
save(strcat('myPCE_LARS_',name,'_',int2str(parcel_idx)),'myPCE_LARS')

% Calculate and save Sobol indices
%SobolOpts.Type = 'Sensitivity';
%SobolOpts.Method = 'Sobol';
%mySobolAnalysisPCE = uq_createAnalysis(SobolOpts);
%mySobolResultsPCE = mySobolAnalysisPCE.Results;
%save(strcat('mySobolResultsPCE_LARS_',name,'_',...
%     int2str(parcel_idx)),'mySobolResultsPCE')

end

