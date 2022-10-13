% script_regulondb
%
% We used to try to convert protein names into gene names. Now we keep a
% separate look-up table.
%
% This script is an attempt to use the regulondb to search for motifs in e
% coli.

clear
close all


% =========================================================================
% Load-in the table look-up between gene(s) and TF complex(es)
% =========================================================================
%{

%
% First we load in the file that cross-references gene names with their
% protein products. This is essentially a look-up table that we'll use
% every time.
%
[~,textdata] = xlsread('RegulonDB/TFset.xlsx');
proteinnames = textdata(2:end,2);
genenames = textdata(2:end,3);

%
% Replace hyphens in protein names with underscores
%
hyphen = strfind(proteinnames,'-');
for i = 1:length(proteinnames)
	proteinnames{i}(hyphen{i}) = '_';
end

%
% Create a second set of genenames and proteinnames that have each gene in
% a list expanded out (some proteins are actually complexes that come from
% multiple genes)
%
proteinnames2 = proteinnames;
genenames2 = genenames;
semicolon = strfind(genenames,';');
k = find(strfindDU(genenames,';'));
for i = k'
	c = str2cell(genenames2{i},'; ');
	c = repeat_remove_cell(c);
	genenames2(i) = c(1)';
	genenames2(end+1:end+length(c)-1) = c(2:end)';
	proteinnames2(end+1:end+length(c)-1) = proteinnames2(i);
end
G = [proteinnames2,genenames2]; % concatenate


%
% Create two structures: one with protein names for fields (for looking up
% what genes make up the proteins in the TF complex), and the other with
% gene names for fields (reverse look-up). For the gene names
%
CELL_pre = {''};
CELL_pre2 = [proteinnames repmat({CELL_pre},length(proteinnames),1)]';
p2g = struct(CELL_pre2{:});
for i = 1:length(proteinnames)
	c = str2cell(genenames{i},'; ');
	c = repeat_remove_cell(c);
	p2g.(proteinnames{i}) = c;
end

genenames2_rr = repeat_remove_cell(genenames2);
CELL_pre = {''};
CELL_pre2 = [genenames2_rr repmat({CELL_pre},length(genenames2_rr),1)]';
g2p = struct(CELL_pre2{:});
for i = 1:length(genenames2_rr)
	g = genenames2_rr{i};
	v = strfindDU(genenames2,g);
	p = proteinnames2(v);
	g2p.(g) = p';
end


%}

% =========================================================================
% Now the GRN connections and making the matrix of interactions
% =========================================================================
%{

% 
% Load-in the GRN connections
% 
[numdata,txtdata] = xlsread('RegulonDB/network_tf_gene.xlsx');

regulators = txtdata(2:end,1);
regulatees = txtdata(2:end,2);
effect = txtdata(2:end,3);

%
% Even though sometimes the effect column is '+-', sometimes there is a
% double entry, one with a '+' and the other with a '-'. This is likely
% because the confidence is different, or the evidence comes from a
% separate place. In any case, my interest is in combining those two
% nearly-duplicate rows into a '+-' effect.
%
N = length(regulators);
c = cellfun(@strcat,regulators,regulatees,'UniformOutput',false);
V = repeatcheckcell(c,true); % see which indices are repeated
for i = 1:length(V)
	E = effect(V{i});
	if isequal(E,{'+';'-'}) || isequal(E,{'+';'+-'}) || ...
			isequal(E,{'+-';'-'}) || isequal(E,{'+';'+-';'-'})
		effect{V{i}(1)} = '+-';
		V{i}(1) = [];
		
	elseif isequal(E,{'+';'?'}) || isequal(E,{'-';'?'})
		V{i}(1) = [];
	end
end
k = [V{:}];
regulators(k) = [];
regulatees(k) = [];
effect(k) = [];
N = length(regulators);


% v = false(N,1);
% for i = 1:N-1
% 	if strcmp(regulators{i},regulators{i+1}) &&...
% 			strcmp(regulatees{i},regulatees{i+1})
% 		v(i) = true;
% 	end
% end
% k = find(v); % the indices of all the repeats (first in pair)
% V = conseccheck(k); % 
% for i = 1:length(k)-length(V)
% 	
% 	
% end


%
% Replace hyphens in protein names with regulators
%
hyphen = strfind(regulators,'-');
for i = 1:N
	regulators{i}(hyphen{i}) = '_';
end

%
% Make a list of all TFs that appear in "regulators", and all genes that
% appear in "regulatees".
%
TFlist = repeat_remove_cell(regulators);
genelist = repeat_remove_cell(regulatees);
nTFs = length(TFlist);
ngenes = length(genelist);
% Note: proteinnames has 8 more entries than TFlist:
% 'AlpA', 'DecR', 'FrmR', 'RacR', 'XynR', 'YdfH', 'YebK', 'YgiV'
% No big deal. I guess these don't have known regulatory connections?
% xlswrite('TFlist.xlsx',TFlist)

%
% Convert effect to a one or negative one. Additionally: a "2" means both,
% and a "3" means question mark.
%
effects = zeros(N,1);
v = strcmp(effect,'+');
effects(v) = 1;
v = strcmp(effect,'-');
effects(v) = -1;
v = strcmp(effect,'+-');
effects(v) = 2;
v = strcmp(effect,'?');
effects(v) = 3;

%
% Create the matrix, which is nTFs number of rows, ngenes number of cols.
%
M = sparse(nTFs,ngenes);
for i = 1:nTFs
	idx_regulators = find(strcmp(regulators,TFlist{i}));
	r = effects(idx_regulators);
	n = length(r);
	
	idx_regulatees = zeros(1,n);
	for j = 1:n
		
		k2 = find(strcmp(genelist,regulatees{idx_regulators(j)}));
		idx_regulatees(j) = k2;
	end
	
	% 1. Create separate index and value arrays.
	% 2. Call sparse to assemble the index and value arrays.
	M(i,idx_regulatees) = r;
end


%}


%% ========================================================================
% Now that we've created the big matrix M, which has nTFs rows, and ngenes
% columns. There are way more genes than TFs, though, so we will re-write
% this matrix to be nTFs-by-nTFs (plus one, where the final column is a
% generic one that shows TF "i" regulates some non-TF gene).
% =========================================================================
%{

I_ROW = cell(nTFs,1);
idx_tot = 1:nTFs;
M_cell = cell(nTFs,1);
endnode_cols = cell(1,2*ngenes); % store a list of which end-node genes are 
% in which column of M_cell1.
for i = 1:nTFs	
	
	if isnan(idx_tot(i))
		continue
	end
	
	%
	% Go forward and backward in a while loop. We continue until we are
	% finished. 
	%
	count = 1;
	I_row = i;
	done = false;
	I_row1 = cell(nTFs,1);
	M_cell1 = sparse(nTFs,2*ngenes);
	while ~done
		
		i1 = I_row(count);
		if isnan(idx_tot(i1))
			
			%
			% If we have more rows to investigate, then increment count and
			% try the next row. If not, we're done.
			%
			if length(I_row) > count
				I_row1{count} = [];
				count = count + 1;
			else
				done = true;
			end
			continue
		end
		
		%
		% Extract row and record that we've done this row (take this row
		% out of commission)
		%
		M1 = M(i1,:);
		idx_tot(i1) = NaN;
		
		
		
		% -----------------------------------------------------------------
		% Going forward: find the genes that are regulated
		% -----------------------------------------------------------------
		idx_regulatees = find(M1 ~= 0);
		effect_col = full(M1(idx_regulatees));
		downstream_genes = genelist(idx_regulatees);
		n_dg = length(downstream_genes);
		
		%
		% Connect the genes back to the proteins that they code for
		%
		I_col = zeros(1,n_dg); % prep to store ea. interaction
		for j = 1:n_dg
			g = downstream_genes{j};
			
			%
			% Check to see if the gene codes for a transcription factor. If
			% so, it will generate a new search within its own row. If not,
			% then it is an endpoint node.
			%
			if isfield(g2p,g)
				p = g2p.(g);
				
				%
				% The gene may act in may TF complexes, so we will loop
				% through them all.
				%
				for k = 1:length(p)
					
					%
					% If the TF complex is in our TFlist (don't forget that
					% 8 of them are not), then we will record that TF as
					% being one of the new rows to check out in this cycle
					% (ie, add to the end of "I_row"), and also add it to
					% the columns that we will mark as "true" for the
					% current TF in question (ie, add it to the jth elt of
					% "I_col").
					%
					if sum(strcmp(TFlist,p(k))) > 0
						idx = find(strcmp(TFlist,p(k)));
						
						if ~ismember(idx,I_row)
							I_row(end+1) = idx;
						end
						I_col(j) = idx;
					end
				end
			else
				
				v_g = strcmp(endnode_cols,g);
				if any(v_g)
					I_col(j) = find(v_g);
					
				else
					k_col = find(cellfun(@isempty,endnode_cols(nTFs+1:end))) + nTFs;
					k_col = k_col(1);
					endnode_cols{k_col} = g;
					I_col(j) = k_col;
				end
			end
		end
		
		%
		% Clean up I_col and effect_col, then add them to the matrix
		% M_cell1.
		%
% 		[I_col,isort] = sort(I_col);
% 		effect_col = effect_col(isort);
		[I_col,effect_col] = repeat_remove(I_col,effect_col);
		v_zero = I_col == 0;
		I_col(v_zero) = [];
		effect_col(v_zero) = [];
		M_cell1(i1,I_col) = effect_col; %#ok<SPRIX>


		% -----------------------------------------------------------------
		% Going backward: See which rows correspond to the TFs that
		% regulate our current TF. Our current TF may be a protein complex
		% composed of several different gene products.
		% -----------------------------------------------------------------
		g = p2g.(TFlist{i1});
		n_g = length(g);
		for j = 1:n_g
			i2 = find(strcmp(genelist,g(j)));
			M2 = M(:,i2);
			idx_regulators = find(M2 ~= 0);
			
			v = ismember(idx_regulators',I_row); % don't add old members
			I_row = [I_row idx_regulators(~v)'];
		end
		
		
		
		% -----------------------------------------------------------------
		% Check to see if we are done, or, if not, what row to continue
		% with next.
		% -----------------------------------------------------------------
		if length(I_row) > count
			n_rows = 0;
			for j = 1:count
				n_rows = n_rows + length(I_row1{j});
			end
			I_row1{count} = I_row(n_rows+1:end);
			
			count = count + 1;
		else
			done = true;
		end
	end
	
	M_cell{i} = M_cell1;
	I_ROW{i} = I_row;
end




%
% Removing the empty cells
%
v = cellfun(@isempty,I_ROW);
M_cell(v) = [];
I_ROW(v) = [];
n_cell = length(M_cell);

%
% Downsizing each matrix in M_cell
%
n_cols = 1;
for i = 1:n_cell
	s = find(sum(~~M_cell{i}) > 0);
	s = s(end);
	n_cols = max(n_cols,s);
end
for i = 1:n_cell
	A = M_cell{i};
	A(:,n_cols+1:end) = [];
	M_cell{i} = A;
end

%
% Creating one big matrix
%
M_out = sparse(nTFs,n_cols);
for i = 1:n_cell
	k = find(M_cell{i});
	M_out(k) = M_cell{i}(k);
end


%
% Create a new set of submatrices that have no fully zero columns or rows.
% Since these are just for show, we will collapse all endnode genes into
% one column.
%
M_cell1 = cell(n_cell,1);
for i = 1:n_cell
	M = M_cell{i};
	M1 = M(:,1:nTFs);
	M_end = M(:,nTFs+1:end);
	M_end = sum(~~M_end,2);
	M_end = ~~M_end + 0;
	
	sumM1 = sum([~~M1 M_end]);
	sumM2 = sum([~~M1 M_end],2);
	v = sumM1(1:end-1)' == 0 & sumM2 == 0;
	M1(:,v) = [];
	M1(v,:) = [];
	M_end(v) = [];
	
	M_cell1{i} = [M1 M_end];
end

save Mat/M_cell



%}

%% ========================================================================
% Export the individual matrices as jpgs
% =========================================================================
%{

%
% Full matrices
%
for i = 1:n_cell
	figure('visib','off')
	imagesc(M_cell{i})
	print(gcf,['Figs_matrix/',num2strDU(i,3),'.jpg'],'-r300','-djpeg')
	close(gcf)
end

%
% Simplified matrices: only export the "interesting" ones
%
for i = 1:length(M_cell1)
	
	M = M_cell1{i};
	if min(M(:)) < max(M(:)) && ~isequal(M,[0 1]) && ~isequal(M,[1 0])
		figure('visib','off')
		imagesc(M)
		print(gcf,['Figs_matrix2/',num2strDU(i,3),'.jpg'],'-r300','-djpeg')
		close(gcf)
	end
end

%}

%% =========================================================================
% Now that we've gone through the data and isolated the parts that are
% potential loops, from visual inspection, only the first matrix of M_cell
% can possibly have loops. We will focus on that matrix and search it for
% loops/motifs.
% =========================================================================
%{
load Mat/M_cell

%
% Extract the matrix in question
%
M = M_cell{1};

%
% First, we look for loops of size 2. Such a loop will have a one for both
% M(i,j) and M(j,i).
%
M1 = M(:,1:nTFs);
M2 = M1';
Diagprime = ~(eye(size(M2,1))); % eliminate diagonal elements

M33 = M1 & M2 & Diagprime;
M3 = triu(M33); % combine all, and take only upper triangle
% imagesc(M3)
[I,J] = find(M3);


%
% Now that we have a bunch of loops of size 2 (17 of them, to be exact), we
% have to double-check to make sure they aren't artifacts of the fact that
% some of these TFs are complexes, and they may regulate their own genes.
%
for i = 1:length(I)
	TF1 = TFlist{I(i)};
	TF2 = TFlist{J(i)};
	
	g1 = p2g.(TF1);
	g2 = p2g.(TF2);
	
	C = [{TF1} g1 num2str(M1(I(i),J(i))) {TF2} g2 num2str(M1(J(i),I(i)))];
	for j = 1:length(C)
		fprintf('%s,',C{j})
	end
	fprintf('\n')
end

%
% OK, so it turns out that none of those are clearly artifacts, and on top
% of that, several are either definitely negative feedback (signs of
% effects are correct) or may be negative feedback (one of the interactions
% is mixed). Only 5 of the 17 are definitely positive feedback, while three
% have mixed interactions, which means 9 of 17 are definitely negative
% feedback.
%
% The remaining question is, of these 9 (actually, we'll search all 17),
% which are the final node of a FFL? If any, which are I1-FFLs?
%
for i = 1:length(I)
	
	TF1 = TFlist{I(i)};
	TF2 = TFlist{J(i)};
	effect1 = M1(I(i),J(i));
	effect2 = M1(J(i),I(i));
	
	%
	% Try several different scenarios. First filter ensures neg fbk, then
	% the scenarios are nested underneath
	%
	if effect1*effect2 < 0 || effect1 == 2 || effect2 == 2
		if effect1 == 1 || effect2 == -1
			% Then TF1 is Z, TF2 is W
			iz = I(i);
			iw = J(i);
			
		elseif effect1 == -1 || effect2 == 1
			% Then TF1 is W, TF2 is Z
			iz = J(i);
			iw = I(i);
			
		else
			% Then we try both
			iz = [I(i) J(i)];
			iw = [J(i) I(i)];
			
		end
	else
		continue
	end
				
	%
	% Test our scenario
	%
	for j = 1:length(iz)
		
		col1 = M1(:,iz(j));
		col1([I(i),J(i)]) = 0;
% 		col1 = col1(setdiff(1:nTFs,[I(i),J(i)]));
		n_upstream = full(sum(~~col1));
		K = find(col1);
		
		if n_upstream < 2
% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' is not the output of an FFL'])
			continue
% 		else
% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' might be the output of an FFL'])
		end
		
		%
		% Now we have to recursively check the possibilities of upstream
		% genes that regulate Z
		%
		%
		% The interaction matrix should look like the following:
		%		 X	 Y	 Z	 W
		%	X: [ 0	 1	 1	 0]
		%	Y: [ 0	 0	-1	 0]
		%	Z: [ 0	 0	 0	 1]
		%	W: [ 0	 0	-1	 0]
		%
		% Can we find any of the submatrices that look like that?
		%
		for o = 1:n_upstream-1
			for p = o:n_upstream
				
				%
				% Construct the regulatory matrix, try #1
				%
				idx = [K(o) K(p) iz(j) iw(j)];
				M2 = full(M1(idx,idx));
				
				if all(M2(2:end,1) == 0) && all(M2(3:end,2) == 0) && ...
						all(M2(1,2:3) > 0) && all(M2(1:2,4) == 0)
					disp(['TF #',num2str(j),' of pair ',num2str(i),' is the output of an FFL'])
					disp(['X: ',TFlist{o},', Y: ',TFlist{p},', Z: ',TFlist{iz(j)},...
						', W: ',TFlist{iw(j)}])
					disp(M2)
				end
				
				%
				% Construct the regulatory matrix, try #2 (permute 1st and
				% 2nd rows and also 1st and 2nd cols)
				%
				idx = [K(p) K(o) iz(j) iw(j)];
				M2 = full(M1(idx,idx));
				
				if all(M2(2:end,1) == 0) && all(M2(3:end,2) == 0) && ...
						all(M2(1,2:3) > 0) && all(M2(1:2,4) == 0)
					disp(['TF #',num2str(j),' of pair ',num2str(i),' is the output of an FFL'])
					disp(['X: ',TFlist{o},', Y: ',TFlist{p},', Z: ',TFlist{iz(j)},...
						', W: ',TFlist{iw(j)}])
					disp(M2)
				end
			end
		end
		
		
		
		%
		% Same thing as the double-for-loop above, but this time, more
		% relaxed constraints.
		%
		%
		% The interaction matrix should look like the following:
		%		 X	 Y	 Z	 W
		%	X: [ 0	 1	 1	 0]
		%	Y: [ 0	 0	-1	 0]
		%	Z: [ 0	 0	 0	 1]
		%	W: [ 0	 0	-1	 0]
		%
		% Can we find any of the submatrices that look like that?
		%
		for o = 1:n_upstream-1
			for p = o:n_upstream
				
				%
				% Construct the regulatory matrix, try #1
				%
				idx = [K(o) K(p) iz(j) iw(j)];
				M2 = full(M1(idx,idx));
				
				if all(M2(2:3,1) == 0) && M2(3,2) == 0 && ...
						all(M2(1,2:3) > 0)
					
					fprintf('\n\nSecond try:\n')
					disp(['TF #',num2str(j),' of pair ',num2str(i),' is the output of an FFL'])
					disp(['X: ',TFlist{o},', Y: ',TFlist{p},', Z: ',TFlist{iz(j)},...
						', W: ',TFlist{iw(j)}])
					disp(M2)
				end
				
				%
				% Construct the regulatory matrix, try #2 (permute 1st and
				% 2nd rows and also 1st and 2nd cols)
				%
				idx = [K(p) K(o) iz(j) iw(j)];
				M2 = full(M1(idx,idx));
				
				if all(M2(2:end,1) == 0) && all(M2(3:end,2) == 0) && ...
						all(M2(1,2:3) > 0) && all(M2(1:2,4) == 0)
					fprintf('\n\nSecond try:\n')
					disp(['TF #',num2str(j),' of pair ',num2str(i),' is the output of an FFL'])
					disp(['X: ',TFlist{o},', Y: ',TFlist{p},', Z: ',TFlist{iz(j)},...
						', W: ',TFlist{iw(j)}])
					disp(M2)
				end
			end
		end
	end	
end



%}



%% =========================================================================
% Now that we've created the matrices and the lists of TFs and endpoint
% genes, we can export them into an excel file.
% =========================================================================
%{

load Mat/M_cell

%
% Trim the "endnode_cols" char cell array. This array contains the list of
% end-node gene names, but also contains a bunch of empty elements. We will
% remove the empty elements, check for duplicates, and sort. (BTW, when we
% sort, we must also keep track of how the sorting happened, so we can
% permute the cols of M_out).
%
v = cellfun(@isempty,endnode_cols);
endnode_cols(v) = [];
[endnode_sort,isort] = sort(endnode_cols);

n_endnode = length(endnode_sort);
endnode_rr = repeat_remove_cell(endnode_sort);

if length(endnode_rr) ~= n_endnode
	error('You had repeats in your endnode list')
end

M_targets = M_out(:,nTFs+1:end);
M_targets = M_targets(:,isort);
M_out(:,nTFs+1:end) = M_targets;

M_TF = M_out(:,1:nTFs);

A = [[{''} TFlist',endnode_sort];[TFlist num2cell(M_out)]];
xlswrite('TableS1.xlsx',A,'M_total')

A = [[{''} TFlist'];[TFlist num2cell(M_TF)]];
xlswrite('TableS1.xlsx',A,'M_TF')

A = [[{''} endnode_sort];[TFlist num2cell(M_targets)]];
xlswrite('TableS1.xlsx',A,'M_targets')







%}
















%% =========================================================================
% There are two more things we could check regarding the combined FF/FB
% mechanism. The first one is easy: are there any I1-FFLs in which Z turns
% around and activates Y? We can use the same off-diagonal symmetry
% arguments as before (because twis is also a cycle of length 2). The
% second one, if Z negatively regulates itself, will be in the cell below.
% =========================================================================
%{

load Mat/M_cell

%
% Extract the matrix in question
%
M = M_cell{1};

%
% First, we look for loops of size 2. Such a loop will have a one for both
% M(i,j) and M(j,i).
%
M1 = M(:,1:nTFs);
M2 = M1';
Diagprime = ~(eye(size(M2,1))); % eliminate diagonal elements

M33 = M1 & M2 & Diagprime;
M3 = triu(M33); % combine all, and take only upper triangle
% imagesc(M3)
[I,J] = find(M3);


%
% Now that we have a bunch of loops of size 2 (17 of them, to be exact), we
% have to double-check to make sure they aren't artifacts of the fact that
% some of these TFs are complexes, and they may regulate their own genes.
%
for i = 1:length(I)
	TF1 = TFlist{I(i)};
	TF2 = TFlist{J(i)};
	
	g1 = p2g.(TF1);
	g2 = p2g.(TF2);
	
	C = [{TF1} g1 num2str(M1(I(i),J(i))) {TF2} g2 num2str(M1(J(i),I(i)))];
	for j = 1:length(C)
		fprintf('%s,',C{j})
	end
	fprintf('\n')
end

%
% OK, so it turns out that none of those are clearly artifacts, and on top
% of that, several are either definitely negative feedback (signs of
% effects are correct) or may be negative feedback (one of the interactions
% is mixed). Only 5 of the 17 are definitely positive feedback, while three
% have mixed interactions, which means 9 of 17 are definitely negative
% feedback.
%
% The remaining question is, of these 9 (actually, we'll search all 17),
% which are the final node of a FFL? If any, which are I1-FFLs?
%
for i = 1:length(I)
	
	TF1 = TFlist{I(i)};
	TF2 = TFlist{J(i)};
	effect1 = M1(I(i),J(i));
	effect2 = M1(J(i),I(i));
	
	%
	% Try several different scenarios. First filter ensures neg fbk, then
	% the scenarios are nested underneath
	%
	if effect1*effect2 < 0 || effect1 == 2 || effect2 == 2
		if effect1 == 1 || effect2 == -1
			% Then TF1 is Z, TF2 is Y
			iz = I(i);
			iy = J(i);
			
		elseif effect1 == -1 || effect2 == 1
			% Then TF1 is Y, TF2 is Z
			iz = J(i);
			iy = I(i);
			
		else
			% Then we try both
			iz = [I(i) J(i)];
			iy = [J(i) I(i)];
			
		end
	else
		continue
	end
				
	%
	% Test our scenario
	%
	for j = 1:length(iz)
		
		col1 = M1(:,iz(j));
		col1([I(i),J(i)]) = 0;
% 		col1 = col1(setdiff(1:nTFs,[I(i),J(i)]));
		n_upstream = full(sum(~~col1));
		K = find(col1);
		
		if n_upstream < 1
% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' is not the output of an FFL'])
			continue
% 		else
% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' might be the output of an FFL'])
		end
		
		%
		% Now we have to recursively check the possibilities of upstream
		% genes that regulate Z
		%
		%
		% The interaction matrix should look like the following:
		%		 X	 Y	 Z
		%	X: [ -	 1	 1]
		%	Y: [ 0	 -	-1]
		%	Z: [ 0	 1	 -]
		%
		% Can we find any of the submatrices that look like that?
		%
		for o = 1:n_upstream
			
			%
			% Construct the regulatory matrix
			%
			idx = [K(o) iy(j) iz(j)];
			M2 = full(M1(idx,idx));
			
			if all(M2(2:end,1) == 0) && all(M2(1,2:3) > 0)
				disp(['TF #',num2str(j),' of pair ',num2str(i),' is the output of an FFL'])
				disp(['X: ',TFlist{o},', Y: ',TFlist{iy(j)},', Z: ',TFlist{iz(j)}])
				disp(M2)
			end
			
		end
		
	end
		
			
end



%}




%% =========================================================================
% There are two more things we could check regarding the combined FF/FB
% mechanism. The first one was above. The second one is if Z negatively
% regulates itself. Here we will have way more than 17 things to check out,
% but the idea is still the same.
% =========================================================================
% {

load Mat/M_cell

%
% Extract the matrix in question
%
M1 = M_out(:,1:nTFs);

%
% First, we look for TFs that negatively autoregulate.
%
d = diag(M1);
M3 = d == -1 | d == 2 | d == 3;
I = find(M3);


%
% OK, so it turns out we have 94 TFs that negatively autoregulate
% (actually, some of those are mixed, and some are "?"). The remaining
% question is, of those 94, which are the final node of at least one
% I1-FFL?
%
for i = 1:length(I)
		
	%
	% Test our scenario
	%
	col1 = M1(:,I(i));
	col1(I(i)) = 0;
	
	n_upstream = full(sum(~~col1));
	K = find(col1);
	
	if n_upstream < 2
		% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' is not the output of an FFL'])
		continue
		% 		else
		% 			disp(['TF #',num2str(j),' of pair ',num2str(i),' might be the output of an FFL'])
	end
	
	%
	% Now we have to recursively check the possibilities of upstream
	% genes that regulate Z
	%
	%
	% The interaction matrix should look like the following:
	%		 X	 Y	 Z
	%	X: [ -	 1	 1]
	%	Y: [ 0	 -	-1]
	%	Z: [ 0	 0	-1]
	%
	% Can we find any of the submatrices that look like that?
	%
	for o = 1:n_upstream
		for p = 1:n_upstream
			
			if o == p
				continue
			end
			
			%
			% Construct the regulatory matrix
			%
			idx = [K(o) K(p) I(i)];
			M2 = full(M1(idx,idx));
			
			if all(M2(2:end,1) == 0) && all(M2(1,2:3) > 0) && ...
					M2(2,3) == -1 % && M2(3,2) == 0
				disp(['TF #',num2str(i),' is the output of an FFL'])
				disp(['X: ',TFlist{o},', Y: ',TFlist{p},', Z: ',TFlist{I(i)}])
				disp(M2)
			end
		end
	end
	
	
	
	
end



%}

