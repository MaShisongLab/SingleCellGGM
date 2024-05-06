classdef SingleCellGGM
	properties
		gene_num = [];
		gene_name = [];
		pcor_all = [];
		pcor_sampling_num = [];
		PCC = [];
		coexpressed_cell_num = [];
		samples_num = [];
		RoundNumber = [];
		SigEdges = [];
		DatasetName = [];
	end

	methods (Static)
		function obj=SingleCellGGM (x, round_num, gene_name, dataset_name)
			
			cut_off_pcor = 0.03;
			cut_off_coex_cell = 10;
			selected_num = 2000;

			if size(x,2) ~= size(gene_name,1)
				error('Error encountered with the inputs!\nThe number of columns in the expression matrix should be equal the number of gene names.\nNumber of columns: %d\nNumber of gene names: %d\nPlease check and run again.\n',size(x,2),size(gene_name,1));
			end

			data_name = 'na';
			if nargin > 3
				data_name = dataset_name;
			end

			[n,p] = size(x);
			
			a = x > 0;
			if issparse(a)
				a = full(a);
			end
			cellnum = sum(a);
			cellnum = cellnum';
			a = double(a);
			coex = a' * a;
			clearvars a;

			if issparse(x)
				x = full(x);
			end
			x = x - mean(x);
			cov_all = x' * x / (n - 1);
			clearvars x;
			rho = corrcov( cov_all );

			gene_id = regexp(gene_name,'([a-zA-Z0-9_\.-]+)','once','match');
			
			obj.samples_num = n;
			obj.gene_num = p;
			obj.RoundNumber = round_num;
			obj.gene_name = gene_id;

			pcor_all = ones( p, p );
			pcor_sampling_num = zeros( p, p );

			time_trend = zeros(100,1);

			rng(98);
			
			fprintf('Calculating pcor in %d iterations.\n', round_num);	

			for i = 1 : round_num
				loop_start_t = clock;
				tic;
				j = datasample(1:p, selected_num, 'Replace', false);
				
				cov_x = cov_all(j,j);
				ix = inv( cov_x );
				d = diag(sqrt(diag(ix)));
				d = inv( d );
				pc = - d * ix * d;
				pc = eye( selected_num ) * 2 + pc;

				for m = 1 : selected_num
					for n = 1: selected_num
						r = j(m);
						s = j(n);
						if r > s
							pcor_sampling_num(r,s) = pcor_sampling_num(r,s) + 1;
							if abs(pc(m,n)) < abs(pcor_all(r,s))
								pcor_all(r,s) = pc(m,n);
							end
						end
					end
				end

				loop_time = etime(clock, loop_start_t);
				idx_time  = mod(i,100) + 1;
				time_trend(idx_time) = loop_time;
				average_loop_time = mean(time_trend);
				time_left = (round_num - i) * average_loop_time / 3600;

				if i == 100
					fprintf('Estimated to complete in %.2f hours.\n', time_left);
				end

				if mod(i,1000) == 0 & i < round_num
					fprintf('%d iterations done. %.2f sec/iteration in average. %.2f hours to go.\n', i, average_loop_time, time_left);
				end

				if i == round_num
					fprintf('%d iterations done.\n',i);
				end

			end

			obj.pcor_sampling_num = int16(pcor_sampling_num);

			idx = find(pcor_sampling_num == 0);
			pcor_all(idx) = 0;

			idx = find(pcor_all >= cut_off_pcor & pcor_all < 1 & coex >= cut_off_coex_cell);
			[e2,e1] = ind2sub(size(pcor_all), idx);
			e1n = cellnum(e1);
			e2n = cellnum(e2);
			e1 = gene_id(e1);
			e2 = gene_id(e2);
			e3 = pcor_all(idx);
			e4 = pcor_sampling_num(idx);
			e5 = rho(idx);
			e6 = coex(idx);
			e7 = e1;
			e7(:) = cellstr(data_name);
			colName = {'GeneA','GeneB','Pcor','SamplingTime','r','Cell_num_A','Cell_num_B','Cell_num_coexpressed','Dataset'};
			obj.SigEdges = table(e1,e2,e3,e4,e5,e1n,e2n,e6,e7,'VariableNames',colName);
			obj.pcor_all = single(pcor_all);
			obj.PCC = single(rho);
			obj.coexpressed_cell_num = int32(coex);
			obj.DatasetName = data_name;
		end
	end
end			
