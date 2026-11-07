function data = getQucsVariable(dataSet,name)
% Get the specified variable from Qucs data set
  data = [];
  if isempty(dataSet)
    error('Data set is empty');
  end
  for i=1:numel(dataSet)
    if strcmp(dataSet(i).name,name)
      data = dataSet(i).data;
      return;
    end
  end
  error(['Variable ' name ' not found in data set']);
end
