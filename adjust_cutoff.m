classdef adjust_cutoff
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
		function obj=adjust_cutoff (ggm_ori, cut_off_pcor, cut_off_coex_cell)
			
			if nargin < 2
				cut_off_pcor = 0.03;
			end

			if nargin < 3
				cut_off_coex_cell = 10;
			end

			dataset_name = ggm_ori.DatasetName;

			coex = ggm_ori.coexpressed_cell_num;
			cellnum = diag(coex);

			obj.samples_num = ggm_ori.samples_num;
			obj.gene_num = ggm_ori.gene_num;
			obj.RoundNumber = ggm_ori.RoundNumber;
			obj.gene_name = ggm_ori.gene_name;
			obj.DatasetName = ggm_ori.DatasetName;

			pcor_sampling_num = ggm_ori.pcor_sampling_num;
			pcor_all = ggm_ori.pcor_all;
			gene_id = ggm_ori.gene_name;
			rho = ggm_ori.PCC;

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
			obj.pcor_all = single(pcor_all);
			obj.pcor_sampling_num = pcor_sampling_num;
			obj.PCC = single(rho);
			obj.coexpressed_cell_num = int32(coex);
		end
	end
end			
