% script_regulondb
%
% This script is "old2" because the matrices in M_cell used to be nTFs by
% (nTFs+1). Now we are allowing there to be ngenes cols (ie, the whole gene
% set number of cols, even though most of those are just outputs).
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




% =========================================================================
% Now the GRN connections and making the matrix of interactions
% =========================================================================

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





% =========================================================================
% Now that we've created the big matrix M, which has nTFs rows, and ngenes
% columns. There are way more genes than TFs, though, so we will re-write
% this matrix to be nTFs-by-nTFs (plus one, where the final column is a
% generic one that shows TF "i" regulates some non-TF gene).
% =========================================================================

I_ROW = cell(nTFs,1);
idx_tot = 1:nTFs;
M_cell = cell(nTFs,1);

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
	M_cell1 = sparse(nTFs,nTFs+1);
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
				I_col(j) = nTFs+1;
			end
		end
		I_col = repeat_remove(I_col);
		I_col(I_col == 0) = [];
		M_cell1(i1,I_col) = 1; %#ok<SPRIX>


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
% Creating one big matrix
%
M_out = sparse(nTFs,nTFs+1);
for i = 1:nTFs
	k = find(M_cell{i});
	M_out(k) = 1;
end

%
% Removing the empty cells
%
v = cellfun(@isempty,I_ROW);
M_cell(v) = [];
I_ROW(v) = [];

%
% Create a new set of submatrices that have no fully zero columns or rows
%
M_cell1 = cell(length(M_cell),1);
for i = 1:length(M_cell)
	M = M_cell{i};
	sumM1 = sum(M);
	sumM2 = sum(M,2);
	v = sumM1(1:end-1)' == 0 & sumM2 == 0;
	M(:,[v; false]) = [];
	M(v,:) = [];
	
	M_cell1{i} = M;
end







%% ========================================================================
% Export the individual matrices as jpgs
% =========================================================================
%{

%
% Full matrices
%
for i = 1:length(M_cell)
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

%
% Extract the matrix in question
%
M = M_cell{1};

%
% First, we look for loops of size 2. Such a loop will have a one for both
% M(i,j) and M(j,i).
%
M1 = M(:,1:end-1);
M2 = M1';
Diagprime = ~(eye(size(M2,1))); % eliminate diagonal elements

M3 = triu(M1 & M2 & Diagprime); % combine all, and take only upper triangle
imagesc(M3)
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
	
	C = [{TF1} g1 {TF2} g2];
	for j = 1:length(C)
		fprintf('%s,',C{j})
	end
	fprintf('\n')
end













