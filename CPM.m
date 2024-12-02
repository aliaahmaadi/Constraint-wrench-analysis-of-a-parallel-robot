function [AA] = CPM(A)
%CPM Summary of this function goes here
%   cross-product-matrix
AA=[0 -A(3) A(2);A(3) 0 -A(1);-A(2) A(1) 0];
end

