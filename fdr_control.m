classdef fdr_control
	properties
		gene_num = [];
		gene_name = [];
		samples_num = [];
		RoundNumber = [];
		SigEdges = [];
		fdr = [];
		DatasetName = [];
	end

	methods (Static)
		function obj=fdr_control (x, ggm_ori, permutation_fraction)
			
			cut_off_pcor = 0.02;
			cut_off_coex_cell = 10;
			selected_num = 2000;

			if nargin < 3
				permutation_fraction = 1;
			end

			dataset_name = ggm_ori.DatasetName;
			obj.DatasetName = dataset_name;

			[n,p] = size(x);

			permutation_num = round(p * permutation_fraction);

			rng('shuffle');
			permutation_id = datasample(1:p, permutation_num,'Replace',false);
			gene_id = ggm_ori.gene_name;
			round_num = ggm_ori.RoundNumber;

			for i = 1:permutation_num
				j = permutation_id(i);
				x(:,j) = x(randperm(n),j);
				gene_id(j) = cellstr(strcat('P_',gene_id(j)));
			end
			
			aa = x > 0;
			cellnum = sum(aa);
			cellnum = cellnum';
			aa = double(aa);
			coex = aa' * aa;
			clear aa;

			cov_all = cov(x);
			clearvars x;
			rho = corrcov( cov_all );
			
			obj.samples_num = n;
			obj.gene_num = p;
			obj.RoundNumber = round_num;
			obj.gene_name = gene_id;

			pcor_all = ones( p, p );
			pcor_sampling_num = zeros( p, p );

			time_trend = zeros(100,1);

			rng('shuffle');
			
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
			e7(:) = cellstr(dataset_name);
			colName = {'GeneA','GeneB','Pcor','SamplingTime','r','Cell_num_A','Cell_num_B','Cell_num_coexpressed','Dataset'};
			obj.SigEdges = table(e1,e2,e3,e4,e5,e1n,e2n,e6,e7,'VariableNames',colName);

			fdr_stat = zeros(91,5);
			idx_np = ones(p,p);
			idx_np(:,permutation_id) = 0;
			idx_np(permutation_id,:) = 0;
			idx_np = idx_np == 1;
			idx_p = ~idx_np;

			num_np = ( p - permutation_num)  * ( p - permutation_num - 1) / 2;
			num_permutated = p * ( p - 1)/2 - num_np;

			for i = 10:100
				idx_ori = ggm_ori.pcor_all >= i/1000 & ggm_ori.coexpressed_cell_num >= cut_off_coex_cell;
				idx_permutated = pcor_all >= i/1000 & coex >= cut_off_coex_cell;
				num_ori_sig = sum(sum(idx_ori));
				num_permutated_sig = sum(sum(idx_permutated & idx_p));
				pct1 = num_permutated_sig / num_permutated;
				pct2 = p * (p - 1) / 2 * pct1 / num_ori_sig;
				if pct2 > 1
					pct2 = 1;
				end
				fdr_stat(i - 9 ,1) = i / 1000;
				fdr_stat(i - 9 ,2) = num_ori_sig;
				fdr_stat(i - 9 ,3) = pct2;
				fdr_stat(i - 9 ,4) = num_permutated_sig;
				fdr_stat(i - 9 ,5) = pct1;
			end
			fdr_str = cellstr(num2str(fdr_stat(:,3),3));
			minFDR = 0;
			for i = 1:91
				if fdr_stat(i,3) > 0
					minFDR = fdr_stat(i,3);
					fdr_str(i) = cellstr(strcat("  ",fdr_str(i)));
				else
					fdr_str(i) = cellstr(strcat("< ",num2str(minFDR,3)));
				end
			end
			colName = {'Pcor','SigEdgeNum','FDR','SigPermutatedEdgeNum','SigPermutatedEdgeProportion'};
			obj.fdr = table(fdr_stat(:,1), num2str(fdr_stat(:,2),9),fdr_str,num2str(fdr_stat(:,4),9),num2str(fdr_stat(:,5),3),'VariableNames',colName);
		end
	end
end			
