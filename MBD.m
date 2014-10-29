function [sol,NP] = MBD(problem,param)
%Find a minimal branching decomposition given the network stoichiometry information.
%function input:
%problem is a structure with the following fields:
%   CbModel: Cobra model at least with the following field:
%               S:      stoichiometric matrix
%               rev:    reversibility of reactions
%   flux: flux distribution
%
%   (optional fields in 'problem':)
%   EFM: EFM matrix 
%       If not supplied, it will be calculated using EFMtool (freely available at http://www.csb.ethz.ch/tools/efmtool). 
%       This requires the CalculateFluxModes.m of EFMtool is in the Matlab
%       search path.
%   imb: metabolites not to be balanced when calculating the EFM matrix, 
%         logical vector: 1 for imbalance, 0 for balance
%         or index vector: index of rows in the stoichiometric matrix
%         (Note: quite impossible to compute all the EFM if all are balanced in
%         genome-scale metabolic network. And some cofactors sometimes
%         better not balanced for a meaningful decomposition to reveal
%         certain functions of an EFM.)
%   w:  weight vector for metabolites, default to be all 1   
%   (optional fields in 'problem' that can be useful when enumerating 
%       optimal or/and suboptimal solutions:)
%   x0: inital set of EFMs provided, index for EFMs in problem.EFM
%   sol0: already-found solution, which will be block in the computation.
%       It is an Ns0 by Nefm logical matrix, where Ns0 is the number of 
%       already-found solution and Nefm is the number of EFMs. 
%       '1' in the i-th row and j-th column means the j-th EFM is in the
%       i-th solution, i.e. EFM(:, sol0(k, :)) is the set of EFMs of the k-th solution.
%   fval (optional): lower bound for the objective function value 
%
%(optional input)
%param: structure with the following fields:
%   solver: currently supported solver, 'cobra' (default) or 'cplex'
%   eps1:   tolerance for the difference between the sum of fluxes in EFMs and the flux in the original flux distribution.
%   eps2:   tolerance for zeros of the weight (alpha) during optimization
%   other optimization parameter for solveCobraMILP or cplexmilp
%
%function output:
%sol is a structure with the following fields:
%   a: the weight (alpha) of each EFM.
%   Eindex: the index of EFMs in the MBD found with respect to the input
%           EFM matrix or the E0 matrix found
%   E0index: the index of contributable EFMs with respect to the input
%           EFM matrix. It appears only if the EFM matrix is given as input.
%   E0: the set of EFMs contributable to the flux distribution,
%           scaled with respect to problem.flux. It appears only no EFM
%           matrix is given as input.
%   fval: objective function value
%   info: the original solution returned by the solver
%NP: number-of-locally-different-path matrix. the i,j-th entry is the
%   number of locally different paths passing through the i-th metabolite in
%   the j-th EFM.

if exist('param','var')
    if isfield(param,'eps1')
        eps1=param.eps1;
    else
        eps1=10^(-6);
    end
    if isfield(param,'eps2')
        eps2=param.eps2;
    else
        eps2=10^(-8);
    end
    if isfield(param,'solver')
        if strcmpi(param.solver,'cobra') || strcmpi(param.solver,'cplex')
            solver=param.solver;
        else
            warning('Solver not supported. Defaulted Cobra')
            solver='cobra';
        end
    else
        solver='cobra';
    end
else
    eps1=10^(-7);
    eps2=10^(-8);
    solver='cobra';
end

if ~isfield(problem,'imb')
    problem.imb=[];
end
S=problem.CbModel.S; %extract the stoichiometrix matrix [m x Nir+Nre]
rev=logical(problem.CbModel.rev);%extract the reversibility vector
%transform S into irreversible stoichiometric matrix [m x Nir+2*Nre]
%Sir=[S | S_reverse of reversible reactions]
Sir=[S -S(:,rev)];
vir=[problem.flux; -problem.flux(rev)];%transform the flux vector [Nir+2*Nre x 1]
vir(vir<0)=0;%transform the flux vector
S4opti=Sir(:,vir~=0);%reduce the matrix to reactions with non-zero fluxes only for optimization [m x Nnz]
%reduce the matrix to exclude metabolites not to be balanced for calculation of EFMs [m-imb x Nnz]
S4efm=S4opti;
S4efm(problem.imb,:)=[];
x0good=true;
%Calculate the whole set of EFMs and filter out EFMs with a=0
if ~isfield(problem,'EFM')||isempty(problem.EFM)  %EFM not user-supplied
    if strcmpi(solver,'cobra')
        global CBT_MILP_SOLVER
        solverCobra=CBT_MILP_SOLVER;
        clear CBT_MILP_SOLVER
    end
    %get the path of efmtool
    efmtoolpath=which('CalculateFluxModes');
    currentpath=pwd;
    if isempty(efmtoolpath)
        error('Unable to locate efmtools in the search path. Please add the folder ''efmtool'' into Matlab search path');
    end
    efmtoolpath=strsplit(efmtoolpath,'CalculateFluxModes.m');
    efmtoolpath=efmtoolpath{1};
    %go to the folder of efmtool
    cd(efmtoolpath)
    %calculate EFMs that are contributable to problem.flux
    mnet = CalculateFluxModes(full(S4efm),zeros(sum(vir~=0),1));
    %go back to the orginal path
    cd(currentpath)
    EFMir=mnet.efms; %EFM matrix [Nnz x K]
    %reset the global variable for Cobra solver because efmtool appears to
    %disrupt it.
    if strcmpi(solver,'cobra')
        if isempty(solverCobra)
            initCobraToolbox;
        else
            changeCobraSolver(solverCobra,'MILP');
        end
    end
else    %user-supplied EFM matrix
    EFM=problem.EFM; %user supplied EFM matrix, assumed to be [Nir+Nre x K0]
    EFMir=[EFM;-EFM(rev,:)]; %transform it into irreversible [Nir+2*Nre x K0]
    EFMir(EFMir<0)=0;
    %find out EFMs with alpha=0:
    cover=false(size(EFMir,2),1);
    nz=vir>eps2;
    for j=1:size(EFMir,2)
        cover(j)=~any(~nz & EFMir(:,j)>eps2);
    end
    EFMir=EFMir(:,cover);%exclude those EFMs with alpha=0
    EFMir=EFMir(vir~=0,:); %retain rows with non-zero fluxes only [Nnz x K]
    if isfield(problem,'x0')
        covInd=find(cover);
        if islogical(problem.x0) && numel(problem.x0)==size(problem.EFM,2)
            x0Ind=find(problem.x0);
        else
            x0Ind=problem.x0;
        end
        [tf,x0]=ismember(x0Ind(:),covInd(:));
        if any(~tf) || isempty(tf)
            warning('probelm.x0 not valid index of problem.EFM')
            x0good=false;
        end
    end
    if isfield(problem,'sol0')
        if ~isempty(problem.sol0)
            sol0=problem.sol0(:,cover);
        end
    end
end
if ~isfield(problem,'w')||isempty(problem.w)
    w=ones(1,size(S,1));    %default weight vector
else
    w=problem.w(:)';
end

%compute the %number-of-locally-different-path matrix for optimization
E0=EFMir;
v=vir(vir~=0);
S=S4opti;
[~,n]=size(S);
K=size(E0,2);
%scale the EFM matrix E with respect to v
[E,scaleR] = scaleEFM(E0,v);
%compute the objective coefficient vector
[NP] = NPmatrix(S,E0);

%additional constraint from already-found solutions:
if isfield(problem,'sol0')
    if ~isempty(problem.sol0)
        A2=zeros(size(sol0,1),2*K);
        b2=zeros(size(sol0,1),1);
        for j=1:size(sol0,1)
            A2(j,K+1:2*K)=sol0(j,:);
            b2(j)=sum(sol0(j,:))-1;
        end
    else
        A2=[];
        b2=[];
    end
else
    A2=[];
    b2=[];
end
%constraint matrix
EYES=sparse((1:K)',(1:K)',ones(K,1),K,K);
A=[[sparse([E;-E]) sparse([],[],[],2*n,K)];[EYES -EYES];sparse(A2)];
%rhs
b=[([v;-v]+eps1);zeros(K,1);b2];
%objective coefficient vector
c=[zeros(K,1);(w*NP)'];
%lower bound
lb=zeros(2*K,1);
%upper bound
ub=ones(2*K,1);
%optimization sense
osense=1;
%constraint sense
if isfield(problem,'sol0')
    if ~isempty(problem.sol0)
        csense=char([ones(1,2*n)*'L' ones(1,K+size(sol0,1))*'L']);
    else
        csense=char([ones(1,2*n)*'L' ones(1,K)*'L']);
    end
else
    csense=char([ones(1,2*n)*'L' ones(1,K)*'L']);
end
%constraint on objective value
if isfield(problem,'fval')
    A = [ A; -sparse(ones(K,1),(K+1:2*K)',(w*NP)',1,2*K)];
    b = [ b ; -problem.fval];
    csense = [csense 'L'];
end
%initial solution provided
if isfield(problem,'x0') && isfield(problem,'EFM') && x0good
    [w1]=FeasibleDecomp(E(:,x0),v,solver);
    if ~isempty(w1)
        a0=zeros(K,1);
        a0(x0)=1;
        w0=zeros(K,1);
        w0(x0)=w1;
        x0=[w0;a0];
    else
        x0=[];
    end
else
    x0=[];
end
%variable type assignment
vartype=char([ones(1,K)*'C' ones(1,K)*'B']);

%Call solver to solve:
MILPproblem=struct();
switch solver
    case 'cobra'
        MILPproblem.A=A;
        MILPproblem.b=b;
        MILPproblem.c=c;
        MILPproblem.lb=lb;
        MILPproblem.ub=ub;
        MILPproblem.osense=osense;
        MILPproblem.csense=csense;
        MILPproblem.x0=x0;
        MILPproblem.vartype=vartype;
        %call solveCobraMILP.m of the COBRA toolbox, the current solveCobraMILP.m 
        %may have problem if the solver in use is Gurobi 5. You are welcome to
        %contact Joshua Chan (joshua.chan@connect.polyu.hk) for a version
        %specifically corrected for the solver Gurobi 5.
        if exist('param','var')
             info= solveCobraMILPedit(MILPproblem,param);
        else
             info= solveCobraMILPedit(MILPproblem);
        end
    case 'cplex'
        MILPproblem.f=c;
        MILPproblem.Aineq=A;
        MILPproblem.bineq=b;
        MILPproblem.lb=lb;
        MILPproblem.ub=ub;
        MILPproblem.x0=x0;
        MILPproblem.ctype=vartype;
        MILPproblem.options=cplexoptimset;
        MILPproblem.options.TolFun=10^(-8);
        MILPproblem.options.TolRLPFun=10^(-8);
        MILPproblem.options.TolXInteger=10^(-8);
        if exist('param','var')
            fn=fieldnames(param);
            fn=setdiff(fn,{'solver';'eps1';'eps2'});
            for j=1:numel(fn)
                eval(['MILPproblem.options.' fn{j} ...
                    ' = param.' fn{j} ';'])
            end
        end
        [x,fval,~,output] = cplexmilp(MILPproblem) ;
        info=struct('cont',[],'int',[],'full',[],'obj',[],...
            'solver','','stat',[],'origStat',[],'time',[]);
        info.solver='cplex_direct';
        info.time=output.time;
        info.origStat=output;
        if ~isempty(x)
            info.cont=x(1:K);
            info.int=x(K+1:2*K);
            info.full=x;
            info.obj=fval;
        else
            info.cont=[];
            info.int=[];
            info.full=[];
            info.obj=[];
        end
        if output.cplexstatus==101 || output.cplexstatus==102
            info.stat=1;
        elseif output.cplexstatus==118
            info.stat=2;
        elseif output.cplexstatus==119 || output.cplexstatus==103
            info.stat=0;
        else
            info.stat=3;
        end
end
        

a=info.cont;%weight (alpha)
b=a>eps2;%on/off of an EFM (beta)
%get the scaled EFM matrix in the original space:
EFMir2=zeros(size(Sir,2),size(E,2));%[Nir+2*Nre x K]
EFMir2(vir~=0,:)=E;
EFM=EFMir2(1:size(problem.CbModel.S,2),:); %[Nir+Nre x K]
EFM(rev,:)=EFM(rev,:)-EFMir2((size(problem.CbModel.S,2)+1):end,:);
%output solution
sol=struct();
if isfield(problem,'EFM')
    sol.a=a(b>eps1).*scaleR(b>eps1);
    cover=find(cover);
    sol.Eindex=cover(b>eps1);
    sol.E0index=cover;
else
    sol.a=a(b>eps1);
    sol.Eindex=find(b);
    sol.E0=sparse(EFM);
end
sol.fval=info.obj;
sol.info=info;
end

function [NP] = NPmatrix(S,E)
%number of locally-different-path matrix
[m,n] = size(S);
C=zeros(m,n);
P=zeros(m,n);
C(S<0)=1;%consumption matrix
P(S>0)=1;%production matrix
D=logical(E);%non-zero entries in E
NP=(C*D).*(P*D);%number-of-locally-different-path matrix [m x K]
end

function [P,r] = scaleEFM(Pm,fx)
%Scale the EFM matrix Pm up to a flux distribution fx
%i.e. P(:,j)*a <= fx for all j for a<=1 
%if a column is not scalable to fx, then keep the original
P=Pm;
r=zeros(size(Pm,2),1);
for k=1:size(Pm,2)
    l=Pm(:,k)~=0;
    r(k)=min(fx(l)./Pm(l,k));
    if r(k)~=0
        P(:,k)=P(:,k)*r(k);
    end
end
end

function [w]=FeasibleDecomp(P,flux,solver,eps1)
LP=struct();
switch solver
    case 'cobra'
        LP.A=[P; -P];
        LP.b=[flux; -flux] + eps1;
        LP.c=zeros(size(P,2),1);
        LP.lb=zeros(size(P,2),1);
        LP.ub=ones(size(P,2),1);
        LP.osense=1;
        LP.csense=char('L'*ones(2*size(P,1),1));
        sol=solveCobraLP(LP);
        if sol.stat~=1
            w=[];
        else
            w=sol.full;
        end
    case 'cplex'
        LP.f=zeros(size(P,2),1);
        LP.Aineq=[P -P];
        LP.bineq=[flux; -flux] + eps1;
        LP.lb=zeros(size(P,2),1);
        LP.ub=ones(size(P,2),1);
        [sol,~,~,out]=cplexlp(LP);
        if out.cplexstatus~=1
            w=[];
        else
            w=sol;
        end
end
end