function array3D = convertMatrixFromCToMatlab(array1D_fromC, dimX,dimY,dimZ)
%It converts 1D float by float (re im re im...) array to 3D complex matrix

allElements = dimX * dimY * dimZ;

re = array1D_fromC(1:2:2*allElements);
im = array1D_fromC(2:2:2*allElements);

CmplxM = complex(re + 1j* im);

CmplxM3D = reshape(CmplxM,dimY,dimX,dimZ);

array3D = permute(CmplxM3D,[2 1 3]);
end

