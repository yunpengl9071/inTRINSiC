% X should be a column vector
% TFs should be #samples-by-#TF
% gene should be a column vector
function f = lsqL1L2(X,TFs,miRNAs,gene,lambda,lambda1)
    Fvec = X(2:end);
    dimTF = size(TFs);
    dimmiRNA = size(miRNAs);
    Fvec_TF = Fvec(1:dimTF(2),1);
    Fvec_miRNA = Fvec((dimTF(2)+1):end,1);
    f_g_TF = times(repmat(Fvec_TF',dimTF(1),1),TFs);
    logs_TF = log2((1+f_g_TF)./(1+TFs));
    if dimmiRNA(2) > 0
        f_g_miRNA = times(repmat(Fvec_miRNA',dimmiRNA(1),1),miRNAs);
        logs_miRNA = log2((1+f_g_miRNA)./(1+miRNAs));
        ev = log2(X(1)) + 0.9 .* sum(logs_TF,2) + 0.1 .* sum(logs_miRNA,2);
    else
        ev = log2(X(1)) + sum(logs_TF,2);
        %disp(size(ev));
    end
    %disp(size(gene));
    diff = gene - ev;
    err = sum(diff.^2);
    if lambda > 0
    	l2norm = lambda * sum(log2(X).^2);
    else
        l2norm = 0;
    end
	if lambda1 > 0
    	l1norm = lambda1 * sum(abs(log2(Fvec+1E-20)));
    else
        l1norm = 0;
    end
%     grad_0 = (2 * diff * (-1)) / (X(1) * log(2));
%     if dimmiRNA(2) > 0
%         grad_1_mat_TF = -2 * log(2) * repmat(diff,1,dimTF(2)) .* TFs ./ (1+(repmat(Fvec_TF',dimTF(1),1)).*TFs);
%         grad_1_mat_mi = -2 * log(2) * repmat(diff,1,dimmiRNA(2)) .* miRNAs ./ (1+(repmat(Fvec_miRNA',dimmiRNA(1),1)).*miRNAs);
%         grad_1_mat = [0.9 * grad_1_mat_TF 0.1 * grad_1_mat_mi];
% 	else
%         grad_1_mat = -2 * log(2) * repmat(diff,1,dimTF(2)) .* TFs ./ (1+(repmat(Fvec',dimTF(1),1)).*TFs);
%     end
%     
%     grad_1 = sum(grad_1_mat,1);
    f = err + l1norm + l2norm;
%     grad = [grad_0;grad_1'];
end
