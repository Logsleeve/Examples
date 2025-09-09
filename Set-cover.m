function results = product_setcover(folderPath, varargin)
% PRODUCT_SETCOVER Minimal substring set-cover across product files.
% results = product_setcover(folderPath, 'minLen',4,'maxLen',12,'nonOverlap',true)
%
% Each .txt in folderPath is treated as a product file. Each line is a part
% number string. Substrings are contiguous. No wildcards. Substrings used
% for product A are removed if they appear in any other product's strings.
% Greedy set cover selects substrings that cover most uncovered strings.
%
% Outputs:
% results(i).productName
% results(i).parts         - cellstr of parts (unmodified)
% results(i).selected     - cellstr of selected substrings
% results(i).coveredIndex - logical array showing coverage of parts
% results(i).notes        - brief stats

% Parse inputs
p = inputParser;
addRequired(p,'folderPath',@(s)ischar(s)||isstring(s));
addParameter(p,'minLen',3,@(x)isnumeric(x) && x>=1);
addParameter(p,'maxLen',12,@(x)isnumeric(x) && x>=1);
addParameter(p,'nonOverlap',true,@islogical);
addParameter(p,'verbose',false,@islogical);
parse(p,folderPath,varargin{:});
minLen = p.Results.minLen;
maxLen = p.Results.maxLen;
nonOverlap = p.Results.nonOverlap;
verbose = p.Results.verbose;

% Read files
files = dir(fullfile(folderPath,'*.txt'));
if isempty(files)
    error('No .txt files found in %s',folderPath);
end

% Load product parts
nProducts = numel(files);
products = struct('name',[],'parts',[]);
for i=1:nProducts
    fname = fullfile(folderPath,files(i).name);
    txt = readlines(fname,'EmptyLineRule','skip');
    txt = cellstr(strtrim(txt));
    % remove empty lines
    txt(cellfun(@isempty,txt)) = [];
    products(i).name = files(i).name;
    products(i).parts = txt;
end

% Build global index: map substring -> set of (product,index) occurrences
% We will instead compute candidate substrings for each product then filter
% those that appear in any other product.
allStringsByProduct = cell(nProducts,1);
for i=1:nProducts
    allStringsByProduct{i} = products(i).parts(:);
end

% Precompute for quick cross-product membership testing:
% For each other product, create a single concatenated big string with separators.
% But safer: just check substring membership across products using contains.
% Build candidates per product
candidates = cell(nProducts,1); % each cell: struct array with .substr and .covers (indices)
for pIdx = 1:nProducts
    parts = allStringsByProduct{pIdx};
    nParts = numel(parts);
    candMap = containers.Map('KeyType','char','ValueType','any'); % substr -> logical cover vector
    for sIdx = 1:nParts
        s = parts{sIdx};
        L = strlength(s);
        % scan start positions from 1..L
        for startPos = 1:double(L)
            maxPossible = double(L) - startPos + 1;
            for len = minLen:min(maxLen, maxPossible)
                sub = char(extractBetween(s,startPos,startPos+len-1));
                if isempty(sub), continue; end
                if isKey(candMap,sub)
                    vec = candMap(sub);
                    vec(sIdx) = true;
                    candMap(sub) = vec;
                else
                    vec = false(1,nParts);
                    vec(sIdx) = true;
                    candMap(sub) = vec;
                end
            end
        end
    end
    % convert map to struct array
    keysList = candMap.keys;
    kN = numel(keysList);
    arr = struct('substr',cell(kN,1),'covers',cell(kN,1),'score',[]);
    for k=1:kN
        arr(k).substr = keysList{k};
        arr(k).covers = candMap(keysList{k});
        arr(k).score = sum(arr(k).covers);
    end
    candidates{pIdx} = arr;
    if verbose, fprintf('Product %d=%s candidates=%d\n',pIdx,products(pIdx).name,kN); end
end

% Filter candidates that appear in other products (cross-product uniqueness)
for pIdx = 1:nProducts
    otherIdx = setdiff(1:nProducts,pIdx);
    arr = candidates{pIdx};
    keepMask = true(numel(arr),1);
    for k=1:numel(arr)
        sub = arr(k).substr;
        appearsElsewhere = false;
        for j = otherIdx
            % check if any string in product j contains sub
            if any(contains(allStringsByProduct{j}, sub))
                appearsElsewhere = true;
                break
            end
        end
        if appearsElsewhere
            keepMask(k) = false;
        end
    end
    candidates{pIdx} = arr(keepMask);
    if verbose
        fprintf('Product %d=%s filtered candidates -> %d\n',pIdx,products(pIdx).name,sum(keepMask));
    end
end

% Greedy set cover per product
results = repmat(struct('productName','','parts',{{}},'selected',{{}},'coveredIndex',[],'notes',''), nProducts,1);
for pIdx = 1:nProducts
    parts = allStringsByProduct{pIdx};
    nParts = numel(parts);
    arr = candidates{pIdx};
    if isempty(arr)
        % no unique candidates; fallback: leave empty solution
        results(pIdx).productName = products(pIdx).name;
        results(pIdx).parts = parts;
        results(pIdx).selected = {};
        results(pIdx).coveredIndex = false(nParts,1);
        results(pIdx).notes = 'No unique candidates after cross-product filtering.';
        continue
    end
    uncovered = true(1,nParts);
    selectedSubs = {};
    % Precompute covers as logical arrays
    coversMat = false(numel(arr), nParts);
    startsPos = zeros(numel(arr),1); % for tie-breaker: earliest start position among occurrences
    for k=1:numel(arr)
        coversMat(k,:) = arr(k).covers;
        % compute earliest start position among parts where it occurs
        idxs = find(coversMat(k,:));
        posVals = inf;
        for ii = idxs
            % find first occurrence position in that part for tie-breaker
            pos = strfind(parts{ii}, arr(k).substr);
            if ~isempty(pos), posVals = min(posVals, pos(1)); end
        end
        startsPos(k) = posVals;
    end

    % Greedy loop
    while any(uncovered)
        % compute coverage gain for each candidate relative to uncovered
        gain = coversMat * uncovered(:); % number of newly covered
        [bestGain, bestIdx] = max(gain);
        if bestGain==0
            break; % can't cover remaining parts with available candidates
        end
        % tie-breaking: prefer shorter substring, then earlier start pos
        bests = find(gain==bestGain);
        if numel(bests)>1
            lens = cellfun(@strlength, {arr(bests).substr});
            [~,iShort] = min(lens);
            tied = bests(lens==lens(iShort));
            if numel(tied)>1
                [~,idxMinPos] = min(startsPos(tied));
                bestIdx = tied(idxMinPos);
            else
                bestIdx = tied(1);
            end
        end

        pick = arr(bestIdx).substr;
        selectedSubs{end+1,1} = pick; %#ok<AGROW>
        newly = coversMat(bestIdx,:) & uncovered;
        uncovered(newly) = false;

        % If nonOverlap constraint: remove any candidate that is substring/ superstring of chosen
        if nonOverlap
            toRemove = false(numel(arr),1);
            for k=1:numel(arr)
                if k==bestIdx, toRemove(k)=true; continue; end
                a = arr(k).substr; b = pick;
                if contains(a,b) || contains(b,a)
                    toRemove(k)=true;
                end
            end
            % remove and update arrays
            arr(toRemove) = [];
            coversMat(toRemove,:) = [];
            startsPos(toRemove) = [];
            if isempty(arr), break; end
        else
            % remove only chosen candidate (avoid picking again)
            arr(bestIdx).score = -inf;
            coversMat(bestIdx,:) = false;
            startsPos(bestIdx) = inf;
        end
    end

    % final coverage check
    coveredIndex = ~any(uncovered,1) ? true(1,nParts) : ~uncovered; % boolean row
    % If some parts remain uncovered, report which ones
    results(pIdx).productName = products(pIdx).name;
    results(pIdx).parts = parts;
    results(pIdx).selected = selectedSubs;
    results(pIdx).coveredIndex = coveredIndex(:);
    numUncov = sum(~coveredIndex);
    notes = sprintf('Selected %d substrings. Uncovered parts: %d of %d.', numel(selectedSubs), numUncov, nParts);
    results(pIdx).notes = notes;
    if verbose
        fprintf('Product %s: %s\n', products(pIdx).name, notes);
    end
end
end
