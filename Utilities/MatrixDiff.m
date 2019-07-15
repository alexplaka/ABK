function[Diff]=MatrixDiff(Matrix1,Matrix2)
% This function compares 2 matrices (whose second dimension must be the
% same) and returns the 'difference' matrix.  It is preferred that Matrix2
% is larger than Matrix1.  Matched rows are deleted from Matrix2, and
% what remains of it is returned by the function.

% Author: Alex Plakantonakis,   Copyright (c) 2019
% License: GNU GPLv3

if size(Matrix1,2) ~= size(Matrix2,2)
    error('Wrong input - Both matrices must have the same number of columns');
    usage('MatrixDiff(A,B) where B has more than or equal to rows to A');
end

a = size(Matrix1,1);

for x = 1:a
    for y = size(Matrix2,1):-1:1
        if isequal(Matrix1(x,:),Matrix2(y,:))
            Matrix2(y,:)=[];
        end
    end
end

Diff = Matrix2;