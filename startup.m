addpath([pwd,'/demo'])
addpath([pwd,'/pop_spike'])
addpath([pwd,'/temp_files'])
addpath([pwd,'/maxent'])
addpath([pwd,'/maxent_MCMC'])
addpath([pwd,'/dich_gauss'])
addpath([pwd,'/flat_models'])
addpath([pwd,'/util'])
addpath([pwd,'/util/minFunc'])


%compile C_Code  if neceessary
if 1
    cd ./maxent_MCMC/C_Code
    ! make all
    cd ../..
end

