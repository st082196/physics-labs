function unpackstruct(S)
%Export structure fields into workspace as variables

for i = convertCharsToStrings(fieldnames(S)).'
    assignin('base',i,S.(i));
end