% For each subtype
% type = 'Proneural';
function b = batchOpt()
rng(19896454);
sigma = 5;
lambda=0.1;
lambda1=0;
gstart=1;
gend=33;
for tp = {'Mesenchymal'}
	type = char(tp)
    wkdir = '/Users/yunpengl/Data/multilayerNetwork/data/regrOutput/scExprs';
    global geneList;
    global miList;
    global TFExpr;
    global miExpr;
    fid = fopen(strcat(wkdir,'/sc_targetGeneList.txt'),'r');
    InputText = textscan(fid,'%s','delimiter','\n');
    geneList = InputText{1};
    %fid = fopen(strcat(wkdir,'/targetmiRNAList.txt'),'r');
    %InputText = textscan(fid,'%s','delimiter','\n');
    %miList = InputText{1};
    len = size(geneList);
    %lenm = size(miList);
    TFExpr = dlmread(strcat(wkdir,'/',type,'/TF_expr.txt'),'\t',1,1);
    %miExpr = dlmread(strcat(wkdir,'/',type,'/mi_expr.txt'),'\t',1,1);
    global g;
    % Regression for the genes
	for i = gstart:gend
		fprintf(strcat(num2str(i/(gend-gstart+1)*100),'%%\r'));
        tgt = geneList{i};
        g = dlmread(strcat(wkdir,'/',type,'/genes/',tgt,'/exprs.txt'))';
        tfs = dlmread(strcat(wkdir,'/',type,'/genes/',tgt,'/TFs.txt'));
        %mtest = dir(strcat(wkdir,'/',type,'/genes/',tgt,'/miRNAs.txt'));
        %if mtest.bytes ~= 0
        %    miRNAs = dlmread(strcat(wkdir,'/',type,'/genes/',tgt,'/miRNAs.txt'));
        %else
        %    miRNAs = [];
        %end
		miRNAs = [];
        stf = size(tfs);stf = stf(1);
        smi = size(miRNAs);smi = smi(1);
        %if mtest.bytes ~= 0
        %    miRNAs(:,1) = int16(miRNAs(:,1));
        %end
        lenx0 = stf + smi + 1;
        %x0 = ones(lenx0,1);
        %x0((stf+2):end,1) = miRNAs(:,2);
        ub = ones(lenx0,1);
        ub(1:(stf+1),:) = Inf;
        if 0
            x0 = [normrnd(0,sigma,(1+stf),1);-abs(normrnd(0,sigma,smi,1))];
            x0 = 2.^x0;
            [res,fval] = autoOpt00(x0,power(2,TFExpr(tfs,:)'),power(2,miExpr(miRNAs(:,1),:)'),g,ub,lambda,lambda1);
        else
            x0 = normrnd(0,sigma,(1+stf),1);
            x0 = 2.^x0;
            [res,fval] = autoOpt00(x0,power(2,TFExpr(tfs,:)'),[],g',ub,lambda,lambda1);
        end
        fn = strcat(wkdir,'/',type,'/genes/',tgt,'/res_ri_nomi_',char(num2str(lambda)),'_',char(num2str(lambda1)),'.txt');
        dlmwrite(fn,[res;fval]);
    end
    fprintf('\n');
    % Regression for the miRNAs
%     for i = 1:lenm(1)
%         tgt = miList{i};
%         g = dlmread(strcat(wkdir,'/',type,'/miRNAs/',tgt,'/exprs.txt'))';
%         tfs = dlmread(strcat(wkdir,'/',type,'/miRNAs/',tgt,'/TFs.txt'));
%         stf = size(tfs);stf = stf(1);
%         lenx0 = stf + 1;
%         x0 = ones(lenx0,1);
%         ub = ones(lenx0,1);
%         ub(1:(stf+1),:) = Inf;
%         [res,fval] = autoOpt00(x0,power(2,TFExpr(tfs,:)'),[],g,ub,lambda,lambda1);
%         fn = strcat(wkdir,'/',type,'/miRNAs/',tgt,'/res_ri_nomi_',char(num2str(lambda)),'_',char(num2str(lambda1)),'.txt');
%         dlmwrite(fn,[res;fval]);
%     end
end
