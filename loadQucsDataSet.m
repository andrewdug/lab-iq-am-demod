function dataSet = loadQucsDataSet(dataSetFile)
% Load a dataset produced by the Qucs-S / Qucsator circuit simulator.
    dataSet = [];
    fid = fopen(dataSetFile,'r');
    if fid < 0
        fprintf(1,'Data set %s not found\n',dataSetFile);
        return;
    end

    err = 0;
    idata = 0;
    idx = 0;
    numdatapoints = 0;
    while ~feof(fid) && ~err
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        if isempty(line)
            break;
        end
        ll = length(line);
        if line(1) == '<'
            if strcmp(line(1:min(7,ll)),'<indep ')
                out = textscan(line(8:end-1),'%s %d');
                if numel(out) ~= 2
                    err = 1;
                    break;
                end
                name = out{1}{1};
                numdatapoints = out{2};
                idx = idx + 1;
                dataSet(idx).name = name;
                dataSet(idx).nameDep = '-';
                dataSet(idx).dep = 0;
                dataSet(idx).data = zeros(1,numdatapoints);
                idata = 1;
            end
            if strcmp(line(1:min(5,ll)),'<dep ')
                out = textscan(line(5:end-1),'%s %s');
                if numel(out) ~= 2
                    err = 1;
                    break;
                end
                name = out{1}{1};
                nameDep = out{2}{1};
                idx = idx + 1;
                dataSet(idx).name = name;
                dataSet(idx).nameDep = nameDep;
                dataSet(idx).dep = 1;
                dataSet(idx).data = zeros(1,numdatapoints);
                idata = 1;
            end
            if strcmp(line(1:min(8,ll)),'</indep>') || strcmp(line(1:min(6,ll)),'</dep>')
                dataSet(idx).data = dataSet(idx).data(1:idata-1);
                dataSet(idx).len = idata - 1;
                idata = 0;
            end
        else
            if idata > 0
                jloc = strfind(line, 'j');
                if isempty(jloc)
                    [val,cnt] = sscanf(line,'%f');
                    if cnt ~= 1
                        error('Invalid data in data set %s\n', dataSetFile);
                    end
                else
                    line(jloc) = line(jloc-1);
                    line(jloc-1) = ' ';
                    [val,cnt] = sscanf(line,'%f');
                    if cnt ~= 2
                        error('Invalid data in data set %s\n', dataSetFile);
                    end
                    val = complex(val(1), val(2));
                end
                dataSet(idx).data(idata) = val;
                idata = idata + 1;
            end
        end
    end
    fclose(fid);
end
