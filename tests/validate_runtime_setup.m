function result = validate_runtime_setup(config)
%VALIDATE_RUNTIME_SETUP Exercise the actual sparse and Laplace MEX routines.
% Supply the same config struct accepted by setup_tutorial. Requires MONALISA
% already built for the current MATLAB and platform; no raw data are needed.
config=setup_tutorial(config);
matrix=sparse([1 1 2 3 4],[1 3 2 4 1],[1 2 3 4 5],4,4);
operator=bmSparseMat_matlabSparse2bmSparseMat(matrix,[4 1 1],[1 1 1],'gauss',3,[0.61 10]);
operator.cpp_prepare('one_block',[],[]);
input=complex(single(reshape(1:12,4,3)),single(reshape(12:-1:1,4,3)));
expected=matrix*double(input);
actual=bmSparseMat_vec(operator,input,'omp','complex',false);
relativeError=norm(double(actual(:))-expected(:))/norm(expected(:));
assert(relativeError<1e-6,'Parallel sparse multiplication failed.');
mask=false(12,12,12); mask(5:8,5:8,5:8)=true;
start=zeros(size(mask),'single'); start(mask)=1;
laplace=bmImLaplaceEquationSolver(start,mask,1000,1e-4,'omp');
assert(all(isfinite(laplace(:))) && all(laplace(mask)==1),'Laplace boundary values changed.');
assert(max(abs(laplace(:)-1))<0.02,'Harmonic extension of a constant boundary failed.');
before=dir(which('bmSparseMat_cC_oBlock_omp_mex'));
setup_tutorial(config);
after=dir(which('bmSparseMat_cC_oBlock_omp_mex'));
assert(before.datenum==after.datenum && before.bytes==after.bytes,'Setup rebuilt an unchanged cache.');
result=struct('status','passed','platform',computer('arch'),'matlab',version, ...
    'sparseRelativeError',relativeError,'laplaceMaximumError',double(max(abs(laplace(:)-1))), ...
    'sparseMex',which('bmSparseMat_cC_oBlock_omp_mex'), ...
    'laplaceMex',which('bmImLaplaceEquationSolver3_omp_mex'));
disp(result);
end
