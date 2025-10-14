function matrix1D_forC = convertIntMatrixForC(matrix3D)

if isempty(matrix3D) 
    
    matrix1D_forC = matrix3D;
    
else

    sizes = size(matrix3D);

    allElements = sizes(2) * sizes(1) * sizes(3);

    inputMatrix = permute(matrix3D,[2 1 3]);

    matrix1D_forC(1:1:allElements)  = single(inputMatrix);

end

