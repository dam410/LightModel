% Copyright ?2013. The Regents of the University of California (Regents).
% All Rights Reserved. Permission to use, copy, modify, and distribute
% this software and its documentation for educational, research, and
% not-for-profit purposes, without fee and without a signed licensing
% agreement, is hereby granted, provided that the above copyright notice,
% this paragraph and the following two paragraphs appear in all copies,
% modifications, and distributions. Contact The Office of Technology
% Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA
% 94720-1620, (510) 643-7201, for commercial licensing opportunities.
%
% Created by Jonathan T Barron and Jitendra Malik, Electrical Engineering
% and Computer Science, University of California, Berkeley.
%
% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


function Zsmooth = inpaintZ(Z, lambda_grad, lambda_curve, SMOOTH_ALL)

if nargin < 4
  SMOOTH_ALL = 0;
end

if size(Z,3)>1
  
  Zsmooth = {};
  for c = 1:size(Z,3)
    Zsmooth{c} = inpaintZ(Z(:,:,c), lambda_grad, lambda_curve);
  end
  Zsmooth = cat(3, Zsmooth{:});
  
else

  lambda_constrain = 10^3;
  
  Zvalid = ~isnan(Z);
  
  addpath(genpath('./bmorph'));
  fidx = find(conv2(double(~Zvalid), ones(3,3), 'same') > 0);
  
  [fi1, fj1] = ind2sub(size(Z), fidx);
  fi0 = fi1 - 1;
  fi2 = fi1 + 1;
  fj0 = fj1 - 1;
  fj2 = fj1 + 1;
  
  fis = [fi0, fi1, fi2];
  fjs = [fj0, fj1, fj2];
  
  keepi = all(fis <= size(Z,1) & fis >= 1,2);
  keepj = all(fjs <= size(Z,2) & fjs >= 1,2);
  
  idx1 = [sub2ind(size(Z), fi0(keepi), fj1(keepi)), sub2ind(size(Z), fi1(keepi), fj1(keepi)), sub2ind(size(Z), fi2(keepi), fj1(keepi))];
  idx2 = [sub2ind(size(Z), fi1(keepj), fj0(keepj)), sub2ind(size(Z), fi1(keepj), fj1(keepj)), sub2ind(size(Z), fi1(keepj), fj2(keepj))];
  
  idx_curve = [idx1; idx2];
  idx_curve = idx_curve(sum(~isnan(Z(idx_curve)),2) < 3,:);
  
  idx_grad = [idx_curve(:,1:2); idx_curve(:,2:3)];
  idx_grad = idx_grad(sum(~isnan(Z(idx_grad)),2) < 2,:);
  
  if SMOOTH_ALL
    
    Acurve = [conv2mat(size(Z), [-1,2,-1]); conv2mat(size(Z), [-1;2;-1])];
    Agrad = [conv2mat(size(Z), [-1,1]); conv2mat(size(Z), [-1;1])];
    
    bcurve = sparse(size(Acurve,1),1);
    bgrad = sparse(size(Agrad,1),1);
    
    n = numel(Z);
    Aeq = speye([n, n]);
    valid = ~isnan(Z);
    Aeq = Aeq(valid,:);
    beq = Z(valid);
    
    A = [lambda_curve*Acurve; lambda_grad*Agrad; lambda_constrain*Aeq];
    b = [lambda_curve*bcurve; lambda_grad*bgrad; lambda_constrain*beq];
    
    X = A \ b;
    Zsmooth = reshape(full(X), size(Z));
    
  else
    
    keep = find(conv2(double(~Zvalid), ones(5,5), 'same') > 0);
    Zkeep = Z(keep);
    Zkeep_valid = ~isnan(Zkeep);
    
    [junk, idx_curve] = ismember(idx_curve, keep);
    [junk, idx_grad] = ismember(idx_grad, keep);
    
    n = size(Zkeep,1);
    
    Acurve = 2*sparse(1:size(idx_curve,1), idx_curve(:,2), 1, size(idx_curve,1), n) - sparse(1:size(idx_curve,1), idx_curve(:,1), 1, size(idx_curve,1), n) - sparse(1:size(idx_curve,1), idx_curve(:,3), 1, size(idx_curve,1), n);
    bcurve = sparse(size(Acurve,1),1);
    
    Agrad = sparse(1:size(idx_grad,1), idx_grad(:,1), 1, size(idx_grad,1), n) - sparse(1:size(idx_grad,1), idx_grad(:,2), 1, size(idx_grad,1), n);
    bgrad = sparse(size(Agrad,1),1);
    
    Aeq = speye([n, n]);
    Aeq = Aeq(Zkeep_valid,:);
    beq = Zkeep(Zkeep_valid);
    
    A = [lambda_curve*Acurve; lambda_grad*Agrad; lambda_constrain*Aeq];
    b = [lambda_curve*bcurve; lambda_grad*bgrad; lambda_constrain*beq];
    
    %   tic;
    X = A \ b;
    %   toc;
    
    Zsmooth = Z;
    Zsmooth(keep) = X;
    
  end
  
  
end
% figure(fig); imagesc([visualizeZ([Z, Zsmooth]); visualizeNormals([Z, Zsmooth])]); imtight;




% function Zsmooth = inpaintZ(Z)
% 
% addpath(genpath('./bmorph'));
% 
% fig = figure; visualizeZ(Z); drawnow;
% 
% % Z(bmorph(isnan(Z), true([3,3]))) = nan;
% 
% Zvalid = ~isnan(Z);
% 
% fidx = find(bmorph(~Zvalid, true([3,3])));
% 
% [fi1, fj1] = ind2sub(size(Z), fidx);
% fi0 = fi1 - 1;
% fi2 = fi1 + 1;
% fj0 = fj1 - 1;
% fj2 = fj1 + 1;
% 
% fis = [fi0, fi1, fi2];
% fjs = [fj0, fj1, fj2];
% 
% keepi = all(fis <= size(Z,1) & fis >= 1,2);
% keepj = all(fjs <= size(Z,2) & fjs >= 1,2);
% 
% idx1 = [sub2ind(size(Z), fi0(keepi), fj1(keepi)), sub2ind(size(Z), fi1(keepi), fj1(keepi)), sub2ind(size(Z), fi2(keepi), fj1(keepi))];
% idx2 = [sub2ind(size(Z), fi1(keepj), fj0(keepj)), sub2ind(size(Z), fi1(keepj), fj1(keepj)), sub2ind(size(Z), fi1(keepj), fj2(keepj))];
% 
% idx = [idx1; idx2];
% idx = idx(sum(~isnan(Z(idx)),2) < 3,:);
% 
% keep = find(bmorph(~Zvalid, true([5,5])));
% 
% Zkeep = Z(keep);
% [junk, idx] = ismember(idx, keep);
% 
% Zkeep_valid = ~isnan(Zkeep);
% 
% n = size(Zkeep,1);
% 
% Asmooth = 2*sparse(1:size(idx,1), idx(:,2), 1, size(idx,1), n) - sparse(1:size(idx,1), idx(:,1), 1, size(idx,1), n) - sparse(1:size(idx,1), idx(:,3), 1, size(idx,1), n);
% bsmooth = sparse(size(Asmooth,1),1);
% 
% Aeq = speye([n, n]);
% Aeq = Aeq(Zkeep_valid,:);
% beq = Zkeep(Zkeep_valid);
% 
% lambda = 10^8;
% 
% A = [Asmooth; lambda*Aeq];
% b = [bsmooth; lambda*beq];
% 
% tic; X = A \ b; toc
% 
% Zsmooth = Z;
% Zsmooth(keep) = X;
% 
% figure(fig); imagesc([visualizeZ([Z, Zsmooth]); visualizeNormals([Z, Zsmooth])]); imtight;



% function Zsmooth = inpaintZ(Z)
% 
% if nargin < 2
%   MODE = 1;
% end
% 
% addpath(genpath('./bmorph'));
% 
% fig = figure; visualizeZ(Z); drawnow;
% 
% Zvalid = ~isnan(Z);
% 
% Zvalid_cropped = ~bmorph(~Zvalid, true([3,3]));
% Zvalid_cropped(:,[1,end]) = true;
% Zvalid_cropped([1,end],:) = true;
% fidx = find(~Zvalid_cropped);
% 
% [fi1, fj1] = ind2sub(size(Z), fidx);
% fi0 = fi1 - 1;
% fi2 = fi1 + 1;
% fj0 = fj1 - 1;
% fj2 = fj1 + 1;
% 
% idx1 = [sub2ind(size(Z), fi0, fj1), sub2ind(size(Z), fi1, fj1), sub2ind(size(Z), fi2, fj1)];
% idx2 = [sub2ind(size(Z), fi1, fj0), sub2ind(size(Z), fi1, fj1), sub2ind(size(Z), fi1, fj2)];
% 
% idx = [idx1; idx2];
% 
% keep = find(bmorph(~Zvalid, true([5,5])));
% % keep = [1:numel(Zvalid)]';
% 
% Zkeep = Z(keep);
% [junk, idx] = ismember(idx, keep);
% 
% Zkeep_valid = ~isnan(Zkeep);
% 
% n = size(Zkeep,1);
% 
% Asmooth = 2*sparse(1:size(idx,1), idx(:,2), 1, size(idx,1), n) - sparse(1:size(idx,1), idx(:,1), 1, size(idx,1), n) - sparse(1:size(idx,1), idx(:,3), 1, size(idx,1), n);
% bsmooth = sparse(size(Asmooth,1),1);
% 
% Aeq = speye([n, n]);
% Aeq = Aeq(Zkeep_valid,:);
% beq = Zkeep(Zkeep_valid);
% 
% lambda = 10^8;
% 
% A = [Asmooth; lambda*Aeq];
% b = [bsmooth; lambda*beq];
% 
% % if MODE == 1
% %   
%   tic; X = A \ b; toc
% % 
% % elseif MODE == 2
% %   
% %   tic; X = leastAbsolute(A,b, [], 5); toc
% %   
% % elseif MODE == 3
% %   
% %   cvx_begin
% %     variable X(size(Asmooth,2),1)
% %     minimize( norm(Asmooth*X - bsmooth, 2) )
% %     subject to
% %       Aeq * X == beq;
% %   cvx_end
% % 
% % elseif MODE == 4
% %   
% %   cvx_begin
% %     variable X(size(Asmooth,2),1)
% %     minimize( norm(Asmooth*X - bsmooth, 1) )
% %     subject to
% %       Aeq * X == beq;
% %   cvx_end
% %   
% % end
% 
% Zsmooth = Z;
% Zsmooth(keep) = X;
% 
% % [fi,fj] = find(Zsmooth == 0);
% % for ii = 1:size(fi,1)
% %   z = Zsmooth(unique(min(max(1, fi(ii) + [-1:1]), size(Z,1))), unique(min(max(1, fj(ii) + [-1:1]), size(Z,2))));
% %   Zsmooth(fi(ii),fj(ii)) = median(z(z~=0));
% % end
% 
% Zsmooth = medfilt2(Zsmooth, [3,3], 'symmetric');
% 
% figure(fig); imagesc([visualizeZ([Z, Zsmooth]); visualizeNormals([Z, Zsmooth])]); imtight;

