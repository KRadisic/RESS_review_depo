function [] = calculate_OLS_PCE_parcel(parcel_idx, name, XED,YED,p,q,A, output_folder)
cd(output_folder);

MetaOpts.ExpDesign.X = XED; 
MetaOpts.ExpDesign.Y = YED ;

% Degree adaptive PCE is autom. enabled if MetaOpts.Degree is an array
MetaOpts.Bootstrap.Replications = 100;
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'PCE';
MetaOpts.Method = 'OLS';
MetaOpts.Degree = 1:p;
MetaOpts.TruncOptions.qNorm = q;
MetaOpts.TruncOptions.Custom = A;

% Calculate and save the metamodel
myPCE_OLS = uq_createModel(MetaOpts);
uq_print(myPCE_OLS)
save(strcat('myPCE_OLS_',name,'_',int2str(parcel_idx)),'myPCE_OLS')

% Calculate and save Sobol indices
%SobolOpts.Type = 'Sensitivity';
%SobolOpts.Method = 'Sobol';
%mySobolAnalysisPCE = uq_createAnalysis(SobolOpts);
%mySobolResultsPCE = mySobolAnalysisPCE.Results;
%save(strcat('mySobolResultsPCE_OLS_',name,'_',...
%     int2str(parcel_idx)),'mySobolResultsPCE')

end
