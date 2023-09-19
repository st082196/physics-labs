function B = remove(A,i)
%Remove the elements from array
%
%   A [vector | matrix | multidimensional array]
%   Input array
%   If A is multidimensional, it is converted to vector.
%
%   i [vector | matrix | multidimensional array]
%   Linear or logical indices of elements to be removed

B = A;
B(i) = [];