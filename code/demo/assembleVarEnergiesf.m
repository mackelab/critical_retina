function assembleVarEnergiesf(ns, idxRepet, idxTs)

cd /home/marcel/criticalityIterScaling/results/energy/

if length(ns) == 1
 MoEs =   cell(length(idxRepet), length(idxTs)); 
 varE  = zeros(length(idxRepet), length(idxTs));
 for i = idxRepet
   names = strsplit(ls(['VarEnf', num2str(10*ns), 'idxRepet', num2str(i), 'T*']));
   if strcmp(names{end}, '')
       names = names(1:end-1); % last one at leat on the on the cluster always is the empty string ''
   elseif strcmp(names{1}, '')
       names = names(2:end);
   end
       
   idxT = zeros(length(idxTs),1);
   for j = 1:length(names) 
     load(names{j})     % gives MoE, T, Ts, lambdaT, x0, burnIn, nSamplesBatch
     MoEs{i,j}=MoE(:,1:end-1); 
     varE(i,j) = mean(MoE(2,3:end-1)) - mean(MoE(1,3:end-1)).^2; 
     idxT(j) = T;
     clear MoE
   end
   [~, idxT] = sort(idxT); 
   MoEs(i,:) = MoEs(i,idxT);
   varE(i,:)  = varE(i,idxT);
   
 end 
 save(['VarEnf', num2str(10*ns), 'HeatTraces.mat'], 'MoEs', 'varE', 'ns', 'idxRepet', 'Ts')
else % save the results for the full sweep over all network sizes

 MoEs =   cell(length(ns), length(idxRepet), length(idxTs)); 
 varE  = zeros(length(ns), length(idxRepet), length(idxTs));
  for n = ns
   for i = idxRepet
     names = strsplit(ls(['VarEnf', num2str(10*n), 'idxRepet', num2str(i), 'T*']));
     if strcmp(names{end}, '')
       names = names(1:end-1); % last one at leat on the on the cluster always is the empty string ''
     elseif strcmp(names{1}, '')
       names = names(2:end);
     end
     idxT = zeros(length(idxTs),1);
     for j = 1:length(names) % last one at leat on the on the cluster always is the empty string ''
       load(names{j})     % gives MoE, T, Ts, lambdaT, x0, burnIn, nSamplesBatch
       MoEs{n,i,j}=MoE(:,1:end-1); 
       varE(n,i,j) = mean(MoE(2,3:end-1)) - mean(MoE(1,3:end-1)).^2; 
       idxT(j) = T;
       clear MoE
     end
     [~, idxT] = sort(idxT); 
     MoEs(n,i,:) = MoEs(n, i,idxT);
     varE(n,i,:)  = varE(n,i,idxT);
   end 
  end
 save(['VarEFFFHeatTraces.mat'], 'MoEs', 'varE', 'ns', 'idxRepet', 'Ts')
end

end % end function

