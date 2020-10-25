% For each subtype
% type = 'Proneural';
function batchOpt(tp,idx,rmd,lambda2,ngenes)
rng(896454);
sigma = 5;
l1=[0];
l2=[str2double];
gstart=1;
gend=str2num(ngenes);
	type = char(tp)
    wkdir = '/ahg/regevdata/projects/txnRegModeling/regression/';
    global geneList;
    global miList;
    global TFExpr;
    global miExpr;
    fid = fopen(strcat(wkdir,'/targetGeneList.txt'),'r');
    InputText = textscan(fid,'%s','delimiter','\n');
    geneList = InputText{1};
    len = size(geneList);
    lenm = 0;
	if str2num(rmd) == 1
    	TFExpr = dlmread(strcat(wkdir,'/',type,'/CV_',num2str(idx),'/TF_expr.txt'),'\t',1,1);
	elseif str2num(rmd) == 2
		TFExpr = dlmread(strcat(wkdir,'/',type,'/randShfl_',num2str(idx),'/TF_expr.txt'),'\t',1,1);
	else
		TFExpr = dlmread(strcat(wkdir,'/',type,'/TF_expr.txt'),'\t',1,1);
	end
    global g;
    outfn = strcat(type,'_progresslog_',num2str(tp),'_rmd',num2str(rmd),'_idx',num2str(idx),'.txt');
    outfid = fopen(outfn,'a');
    % Regression for the genes
	for i = gstart:gend
		fwrite(outfid,([num2str(i/(gend-gstart+1)*100) '%' 10]),'char');
        tgt = geneList{i};
		if str2num(rmd) == 1
			g = dlmread(strcat(wkdir,'/',type,'/CV_',num2str(idx),'/genes/',tgt,'/exprs.txt'))';
			tfs = dlmread(strcat(wkdir,'/',type,'/CV_',num2str(idx),'/genes/',tgt,'/TFs.txt'));
		elseif str2num(rmd) == 2
			g = dlmread(strcat(wkdir,'/',type,'/randShfl_',num2str(idx),'/genes/',tgt,'/exprs.txt'))';
			tfs = dlmread(strcat(wkdir,'/',type,'/randShfl_',num2str(idx),'/genes/',tgt,'/TFs.txt'));
		else
			g = dlmread(strcat(wkdir,'/',type,'/genes/',tgt,'/exprs.txt'))';
			tfs = dlmread(strcat(wkdir,'/',type,'/genes/',tgt,'/TFs.txt'));
		end
        miRNAs = [];
        stf = size(tfs);stf = stf(1);
        smi = size(miRNAs);smi = smi(1);
        lenx0 = stf + smi + 1;
        ub = ones(lenx0,1);
        ub(1:(stf+1),:) = Inf;
		for lm=[l1';l2']
            lambda = lm(2);
            lambda1 = lm(1);
            bestfval = Inf(1);
            for k=[1:1:1]
                x0 = normrnd(0,sigma,(1+stf),1);
                x0 = 2.^x0;
                [res,fval] = autoOpt00(x0,power(2,TFExpr(tfs,:)'),[],g',ub,lambda,lambda1);
                if fval < bestfval
                    resv = res;
                    bestfval = fval;
                end
            end
            if str2num(rmd) == 1
            	fn = strcat(wkdir,'/',type,'/CV_',num2str(idx),'/genes/',tgt,'/res_ri_nomi_',char(num2str(lambda)),'_',char(num2str(lambda1)),'.txt');
            elseif str2num(rmd) == 2
                fn = strcat(wkdir,'/',type,'/randShfl_',num2str(idx),'/genes/',tgt,'/res_ri_nomi_',char(num2str(lambda)),'_',char(num2str(lambda1)),'.txt');
            else
                fn = strcat(wkdir,'/',type,'/genes/',tgt,'/res_ri_nomi_',char(num2str(lambda)),'_',char(num2str(lambda1)),'.txt');
            end
            dlmwrite(fn,[resv;bestfval]);

		end
    end
    fprintf('\n');
    fclose(outfid);
   % if isdeployed
   %    exit;
   % else
   %     close all
   % end
end
